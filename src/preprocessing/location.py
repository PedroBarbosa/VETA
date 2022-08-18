import hgvs.parser
import hgvs.enums
import hgvs.exceptions
import numpy as np
import pandas as pd


CONSEQUENCES = {"missense_variant": ['coding', 'missense'],
                "synonymous_variant": ['coding', 'synonymous'],
                "frameshift_variant": ['coding', 'frameshift'] , 
                "inframe_deletion": ['coding', 'other_coding'],
                "inframe_insertion": ['coding', 'other_coding'], 
                "incomplete_terminal_codon_variant": ['coding', 'other_coding'],
                "coding_sequence_variant": ['coding', 'other_coding'],
                "protein_altering_variant": ['coding', 'other_coding'] ,
                "start_lost": ['coding', 'other_coding'], 
                "stop_lost": ['coding', 'other_coding'],
                "stop_gained": ['coding', 'other_coding'],
                "start_retained_variant": ['coding', 'other_coding'],
                "stop_retained_variant": ['coding', 'other_coding'],
                "transcript_ablation": ['coding', 'other_coding'],
                
                "splice_acceptor_variant": ['splice_region', 'splice_site'],
                "splice_donor_variant": ['splice_region', 'splice_site'], 
                "splice_donor_region_variant": ['splice_region', 'splice_site'],
                "splice_polypyrimidine_tract_variant": ['splice_region', 'splice_region'],
                "splice_donor_5th_base_variant": ['splice_region', 'splice_region'],
                "splice_region_variant": ['splice_region', 'splice_region'],
                "intron_variant": ['intronic', 'intronic'],
                
                "5_prime_UTR_variant": ['5primeUTR', '5primeUTR'],
                "3_prime_UTR_variant": ['3primeUTR', '3primeUTR'],
                
                "regulatory_region_ablation": ['regulatory', 'regulatory'],
                "regulatory_region_variant": ['regulatory', 'regulatory'],
                "TF_binding_site_variant": ['regulatory', 'regulatory'],
                "TFBS_ablation": ['regulatory', 'regulatory'],
                "TFBS_amplification": ['regulatory', 'regulatory'],
                
                "non_coding_transcript_exon_variant": ['noncodingRNA', 'noncodingRNA'],
                "non_coding_transcript_variant": ['noncodingRNA', 'noncodingRNA'],
                
                "upstream_gene_variant": ['other_location', 'other_types'],
                "downstream_gene_variant": ['other_location', 'other_types'],
                "intergenic_variant": ['other_location', 'other_types'],
                "mature_miRNA_variant": ['other_location', 'other_types'],
                "NMD_transcript_variant": ['other_location', 'other_types']
}
    
def get_location(x, hp, aggregate: bool):
    """
    Get location of variant from HGVSc expression

    :param pd.Series x: Single row
    :param hp: HGVS parser
    :param aggregate: Whether variants in coding/splice site regions
    should be aggregated and analyzed together 
    
    :return str: location
    """
    import itertools

    try:
        v = hp.parse_hgvs_variant(x['HGVSc'].split(" ")[0])
        if v.type == "m":
            return 'mithocondrial'
        elif v.type == "n":
            return 'noncodingRNA'
        elif v.posedit.pos.start.base < 0:
            return '5primeUTR'
        elif v.posedit.pos.start.datum == hgvs.enums.Datum.CDS_END and v.posedit.pos.start.base > 0:
            return '3primeUTR'
        elif v.posedit.pos.start.base > 0 and v.posedit.pos.start.offset == 0:
            return 'coding' if aggregate else get_location_from_consequence(x.Consequence, aggregate) 
        elif v.posedit.pos.start.base > 0 and abs(v.posedit.pos.start.offset) > 15:
            return get_location_from_consequence(x.Consequence, aggregate) 
        elif v.posedit.pos.start.base > 0 and abs(v.posedit.pos.start.offset) > 100:
            return 'deep_intronic'
        elif v.posedit.pos.start.base > 0 and 3 <= abs(v.posedit.pos.start.offset) <= 15:
            return 'splice_region'
        elif v.posedit.pos.start.base > 0 and 1 <= abs(v.posedit.pos.start.offset) <= 2:
            return 'splice_region' if aggregate else 'splice_site'
        else:
            raise ValueError('Weird HGVSc expression: {}'.format(v))

    except hgvs.exceptions.HGVSParseError:
        return get_location_from_consequence(x.Consequence, aggregate)
   

def get_location_from_consequence(x: str,
                                  aggregate: bool):
    """
    Get location/class of input variant
    
    :param str x: Value of consequence field for a single variant
    :param bool aggregate: Level of detail to report 
    """
    idx = 0 if aggregate else 1
    conseqs = x.split("&")
    
    # splice_region variants are prioritized over synonymous
    if all(x in conseqs for x in ["splice_region_variant",
                                  "synonymous_variant"]):
        conseq = "synonymous_variant"
    else:
        conseq = conseqs[0]

    try:
        return CONSEQUENCES[conseq][idx]
    
    except KeyError:
        raise ValueError("Consequence not registered. {}".format(conseq))

ranges_split_at_ss = [(0, 2, '1-2'),
          (3, 10, '3-10'),
          (11, 40, '11-40'),
          (41, 200, '41-200'),
          (201, 500, '201-500'),
          (501, 1000, '501-1000'),
          (1001, 5000000, '1000+')]

ranges = [(0, 10, '1-10'),
          (11, 40, '11-40'),
          (41, 200, '41-200'),
          (201, 500, '201-500'),
          (501, 1000, '501-1000'),
          (1001, 5000000, '1000+')]


def assign_intronic_bins(hgvs_exp, hp, location, aggregate):
    if location in {"splice_site", "splice_region", "deep_intronic", "intronic",
                    "5primeUTR", "3primeUTR", "noncodingRNA"}:
        try:
            v = hp.parse_hgvs_variant(hgvs_exp.split(" ")[0])
            
        except hgvs.exceptions.HGVSParseError:
            return pd.Series([np.nan, np.nan, np.nan]) 
            
        _offset = v.posedit.pos.start.offset
        offset = abs(_offset)
        if offset <= 1000:
            which_ss = "donor" if _offset > 0 else "acceptor"
        else:
            which_ss = "unclear"
            
        _ranges = ranges if aggregate else ranges_split_at_ss
        for i in _ranges:
            if i[0] <= offset <= i[1]:
                return pd.Series([i[2], np.int(offset), which_ss])
    else:
        return pd.Series([np.nan, np.nan, np.nan])
