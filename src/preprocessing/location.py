import re
import hgvs.parser
import hgvs.enums
import hgvs.exceptions
import numpy as np
import pandas as pd


def get_location(x, hp):
    try:
        v = hp.parse_hgvs_variant(x.split(" ")[0])
        if v.type == "m":
            return 'mithocondrial'
        elif v.type == "n":
            return 'noncodingRNA'
        elif v.posedit.pos.start.base < 0:
            return '5primeUTR'
        elif v.posedit.pos.start.datum == hgvs.enums.Datum.CDS_END and v.posedit.pos.start.base > 0:
            return '3primeUTR'
        elif v.posedit.pos.start.base > 0 and v.posedit.pos.start.offset == 0:
            return 'coding'
        elif v.posedit.pos.start.base > 0 and abs(v.posedit.pos.start.offset) >= 20:
            return 'deepintronic'
        else:
            return 'splicesite'

    except hgvs.exceptions.HGVSParseError:
        return 'unknown'


def get_location_from_consequence(x):
    conseq = x.split("&")[0]
    coding = ["coding_sequence_variant", "missense_variant",
              "frameshift_variant", "inframe_deletion",
              "inframe_insertion", "protein_altering_variant",
              "start_lost", "start_retained_variant",
              "stop_gained", "stop_lost",
              "stop_retained_variant", "synonymous_variant",
              "transcript_ablation"]
    splice_site = ["splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"]
    deep_intronic = ["intron_variant"]
    five_prime_utr = ["5_prime_UTR_variant"]
    three_prime_utr = ["3_prime_UTR_variant"]
    regulatory_variant = ["regulatory_region_ablation", "regulatory_region_variant", "TF_binding_site_variant"]
    noncodingRNA= ["non_coding_transcript_exon_variant", "non_coding_transcript_variant"]
    unknown = ["upstream_gene_variant", "downstream_gene_variant", "intergenic_variant", "mature_miRNA_variant",
               "NMD_transcript_variant"]

    if conseq in coding:
        return 'coding'
    elif conseq in splice_site:
        return 'splicesite'
    elif conseq in deep_intronic:
        return 'deepintronic'
    elif conseq in five_prime_utr:
        return '5primeUTR'
    elif conseq in three_prime_utr:
        return '3primeUTR'
    elif conseq in noncodingRNA:
        return 'noncodingRNA'
    elif conseq in regulatory_variant:
        return 'regulatory_variant'
    elif conseq in unknown:
        return 'unknown'
    else:
        print("Consequence not registered. {}".format(conseq))
        return conseq


ranges = [(0, 10, '0-10'),
          (11, 30, '10-30'),
          (31, 100, '30-100'),
          (101, 200, '100-200'),
          (201, 500, '200-500'),
          (501, 5000000, '500+')]


def assign_intronic_bins(hgvs, hp, location):

    if location in {"splicesite", "deepintronic", "5primeUTR", "3primeUTR", "noncodingRNA"}:
        v = hp.parse_hgvs_variant(hgvs.split(" ")[0])
        offset = abs(v.posedit.pos.start.offset)
        for i in ranges:
            if i[0] <= offset <= i[1]:
                return pd.Series([i[2], np.int(offset)])
    else:
        return pd.Series([np.nan, np.nan])
