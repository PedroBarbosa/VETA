import logging
import multiprocessing
import os
import subprocess
import tempfile
from collections import OrderedDict, defaultdict
from functools import partial
import re
import pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm
from preprocessing import osutils


def process_vcf(vcf: str,
                tools_config: dict,
                thresholds: list,
                select_conseq: str,
                af_col: str,
                is_clinvar: bool = False):
    """
    Scans a VEP annotated VCF file and extracts all the
    tools scores

    :param str vcf: VCF file
    
    :param dict tools_config: Dict with the
        map between available tools and VCF
        annotations (already parsed, only
        existing fields should be passed here).
        
    :param list thresholds: List of tools
        to process based on the tools config file
        as well as the selected scope for the
        analysis
        
    :param str select_conseq: Strategy to select
    a single consequence block per variant 
    
    :param str af_col: VCF field that measures
            population frequency of variants.
    :param bool is_clinvar: Whether `vcf` is
        from clinvar database.

    :return dict: Dict with scores keyed by variant
    :return list: List with the tools provided in the
    config that were not found in the VCF. (Will be
    used to update the config_tools variable)

    """
    logging.info("Extracting predictions from VCF.")

    # List of required INFO fields when
    # dataset is from Clinvar
    CLINVAR_FIELDS = ['CLNREVSTAT', 'CLNSIG']
    clinvar_confirmed = False

    # List of tools belonging to the
    # provided scope that will be analysed
    # Requires to be present in the updated config
    # (some tools might be removed from config after processing first file)
    TOOLS_TO_ANALYSE = [t[0] for t in thresholds if t[0] in tools_config.keys()]

    # MAP with the list of VCF fields per tool
    TOOLS_CONFIG_MAP = {k: v[0] for k, v in tools_config.items() if k in TOOLS_TO_ANALYSE}

    # List with the VCF fields to check
    # Present fields in the header will be
    # popped out sequentially
    MISSING_VCF_FIELDS_FLAT = [i for sublist in [v for _k, v in TOOLS_CONFIG_MAP.items()]
                               for i in sublist]

    # Missing tools from VCF to return
    # and exclude from analysis
    MISSING_VCF_TOOLS = []

    # Indexes is an ordered dict with the index(s)
    # where information about a given tool is stored
    # in the VEP annotation field
    vep_indexes = OrderedDict()

    # List of VCF fields (not tool names) that are not
    # found in the VEP annotation field
    absent_from_VEP_annot = []

    # Dict with the list of tools and corresponding
    # VCF fields found in the INFO field as a single
    # annotation
    present_in_INFO = defaultdict(list)

    # List of tool names with all annotation fields
    # found in the VEP annotation field
    all_vep_annotations = []
    present_in_VEP_annot = []

    # Parse header and count number of variants
    vcf_data = VCF(vcf)
    if vcf_data.contains("ANN"):
        VEP_TAG = "ANN"
    elif vcf_data.contains("CSQ"):
        VEP_TAG = "CSQ"
    else:
        raise ValueError("VETA requires VEP annotations. ANN or CSQ field was not found in "
                         "the INFO field of the input VCF file. Exiting.\n")

    for field in vcf_data.header_iter():
      
        # If field is in VEP annotations
        if field["HeaderType"] == "INFO" and field["ID"] == VEP_TAG:
            all_vep_annotations = field["Description"].split("Format:")[1][:-1].strip().split("|")
            if af_col in all_vep_annotations:
                vep_indexes[af_col] = [all_vep_annotations.index(af_col)]

            # Looks only for scores belonging
            # to the specific scope set
            # Also look for the AF column set
            for _tool in TOOLS_TO_ANALYSE:

                _tool_field = TOOLS_CONFIG_MAP[_tool]
                _present, _absent = _check_if_field_exists(_tool_field, all_vep_annotations)

                if _present:
                    _fields_present = [v[0] for v in _present]
                    present_in_VEP_annot.extend(_fields_present)
                    vep_indexes[_tool] = [v[1] for v in _present]
                    [MISSING_VCF_FIELDS_FLAT.remove(i) for i in _fields_present]

                if _absent:
                    absent_from_VEP_annot.extend(_absent)

        elif field["HeaderType"] == "INFO":
            if field['ID'] in MISSING_VCF_FIELDS_FLAT:
                _tool_name = [_t for _t, v in TOOLS_CONFIG_MAP.items() if field['ID'] in v][0]
                present_in_INFO[_tool_name].append(field['ID'])
                MISSING_VCF_FIELDS_FLAT.remove(field['ID'])

            if field['ID'] in CLINVAR_FIELDS:
                clinvar_confirmed = True

            if field['ID'] == af_col and af_col not in vep_indexes.keys() and af_col not in present_in_INFO.keys():
                present_in_INFO[af_col].append(field['ID'])

    if is_clinvar and not clinvar_confirmed:
        raise ValueError('Input VCF is not from Clinvar ({}) fields were '
                         'not found in the VCF INFO annotations. They are '
                         'required for clinical significance extraction.\n'
                         'Additionally, be aware that when running in '
                         'benchmark mode, non-Clinvar files must be passed '
                         'via an input directory where target files '
                         'referring to each class are located.'.format(CLINVAR_FIELDS))

    # if there are VCF fields from
    # the config absent in VCF header
    # LOG absent fields
    if MISSING_VCF_FIELDS_FLAT:

        for tool, fields in TOOLS_CONFIG_MAP.items():

            for _f in fields:
                if _f in MISSING_VCF_FIELDS_FLAT:
                    logging.info("\'{}\' field does not exist "
                                 "in VCF. Corresponding tool "
                                 "(\'{}\') will be discarded from "
                                 "analysis.".format(_f, tool))

                    MISSING_VCF_TOOLS.append(tool)
                    # Delete tool from analysis
                    # if at least one subfield
                    # is not present
                    try:
                        del vep_indexes[tool]
                    except KeyError:
                        pass

                    try:
                        del present_in_INFO[tool]
                    except KeyError:
                        pass

    n_variants = 0
    for _ in vcf_data:
        n_variants += 1
    vcf_data.close()

    # Pack kwargs
    args = {'VEP_TAG': VEP_TAG,
            'all_vep_annotations': all_vep_annotations,
            'present_in_INFO': present_in_INFO,
            'vep_indexes': vep_indexes,
            'is_clinvar': is_clinvar,
            'select_conseq': select_conseq}
    
    # Variant processing
    if n_variants > 5000:
        header = subprocess.Popen(["bcftools", "view", "-h", vcf],
                                  stdout=subprocess.PIPE).stdout.readlines()

        body_variants = tempfile.NamedTemporaryFile()
        subprocess.call(["bcftools", "view", "-H", vcf], stdout=body_variants)
 
        list_files = osutils.split_file_in_chunks(body_variants.name, header, n_variants)

        with multiprocessing.Pool() as p:
            df_list = list(tqdm(p.imap(partial(_iter_variants, **args), list_files)))

        logging.info("Merging data from parallel VCF processing.")
        scores = pd.concat(df_list)
        [os.remove(f) for f in list_files]
        logging.info("Done.")
        p.close()

    else:
        scores = _iter_variants(vcf, **args)
    return scores, MISSING_VCF_TOOLS


def _iter_variants(vcf: str, **kwargs):
    """
    Iterates over a VCF file

    :param str vcf: vcf file
    :param dict **kwargs: Additional info obtained
        when processing the VCF header

    :return dict: Dictionary with scores for each variant
    """
    scores = []
    keys = []

    vcf_data = VCF(vcf)

    # Iterate over VCF
    for record in vcf_data:
        
        colnames = ['index', 'chr', 'pos', 'ref', 'alt', 'id', 'type', 'subtype']
        _single_out = []
        
        # Only keep variants with VEP annotations
        vep_annotation = record.INFO.get(kwargs['VEP_TAG'])
        if not vep_annotation:
            continue
        
        if record.ID != "." and record.ID is not None:
            key = record.CHROM + "_" + str(record.POS) + "_" + str(record.ID) + "_" + str(record.ALT[0])

        else:
            key = record.CHROM + "_" + str(record.POS) + "_" + str(record.REF) + "_" + str(record.ALT[0])
        
        if key in keys:
            logging.info("Variant {},{},{},{}, {} is repeated in the VCF. "
                         "Keeping only first occurrence".format(record.CHROM,
                                                                record.POS,
                                                                record.REF,
                                                                record.ALT[0],
                                                                record.ID))
            continue

        else:
            keys.append(key)

        _single_out.append(key)

        # MNPs
        if record.var_type == "indel" and len(record.REF) == len(record.ALT[0]):
            _single_out.extend([record.CHROM,
                                record.POS,
                                record.REF,
                                record.ALT[0],
                                record.ID,
                                record.var_type, 
                                "mnp"])

        else:
            _single_out.extend([record.CHROM,
                                record.POS,
                                record.REF,
                                record.ALT[0],
                                record.ID,
                                record.var_type,
                                record.var_subtype])

        if record.POS == 44874309:
            a=1
        vep_info = _select_consequence(vep_annotation, **kwargs)
        # Add some VEP fields, if they exist
        for _field in ['Existing_variation', 'HGVSc', 'HGVSg', 'SYMBOL', 'Consequence']:
            colnames.append(_field)
            try:
                _single_out.append(vep_info[kwargs['all_vep_annotations'].index(_field)])
            except ValueError:
                _single_out.append(None)

        # Look for externally annotated scores
        # that are single fields in the INFO field
        if len(kwargs['present_in_INFO']) > 0:

            for _tool, _fields in kwargs['present_in_INFO'].items():
                colnames.append(_tool)
                tool_out = []
                for _f in _fields:
                    tool_out.append(record.INFO.get(_f))
                _single_out.append(tool_out)

        # Append fields from VEP annotation
        for _tool, idx in kwargs['vep_indexes'].items():
            colnames.append(_tool)
            # If tool has multiple fields
            # and at least 1 is within VEP annot
            # and other in the INFO annot.
            if _tool in kwargs['present_in_INFO'].keys():
     
                scores_so_far = scores[key]
                _updated = []
                for _v in scores_so_far:
                    if isinstance(_v, tuple) and _v[0] == _tool:
                        val = _v[1] + [vep_info[i] for i in idx]
                        _updated.append((_v[0], val))
                    else:
                        _updated.append(_v)

                scores[key] = _updated

            else:
     
                _single_out.append([vep_info[i] for i in idx])

        if kwargs['is_clinvar']:
            colnames.extend(['CLNREVSTAT', 'CLNSIG', 'CLNDISDB'])
            _single_out.extend([record.INFO.get("CLNREVSTAT"), 
                                record.INFO.get("CLNSIG"),
                                record.INFO.get("CLNDISDB")])
        
        df = pd.DataFrame( _single_out).T
        df.columns = colnames
        scores.append(df)
    vcf_data.close()

    if len(scores) != 0:
        return pd.concat(scores).set_index('index')
    


def _check_if_field_exists(field: list, available: list):
    """
    Checks whether a given field
    provided in the tools config (2nd col)
    exists in the VCF header

    :param list field: VCF field (s) from
    which predictions of a given tool should
    be extracted. If len(field) > 1,  multiple
    fields were provided in the tools config file

    :param list available: List of available fields
    to search for

    :return list: List tuples where 1st element
        is the VCF field name found and the second the
        index where it is found in the list of available
         fields

    :return list: List of VCF field name (s) that weren't
        found in the list of available fields
    """
    present, absent = [], []

    for _f in field:
        if _f in available:
            present.append((_f, available.index(_f)))
        else:
            absent.append(_f)

    return present, absent


def _select_consequence(vep_annotations: str, **kwargs):
    """
    Selects top VEP consequence block based on a
    choice provided by the user.

    :param str vep_annotations: Str with VEP
     annotations. If multiple blocks, "," splits
     them

    :return str: Top block
    """

    if "," in vep_annotations:
        
        if kwargs['select_conseq'] == "first":
            return vep_annotations.split(",")[0].split("|")
        
        elif kwargs['select_conseq'] in ["gene_body", "smallest_offset"]:
       
            try:
                hgvsc_index = kwargs['all_vep_annotations'].index('HGVSc')
            
            except ValueError:
                raise ValueError('HGVSc not present in VEP annotations.'
                                 'Can\'t select top consequence block based '
                                 'on the \'gene_body\' flag')
            
         
            if kwargs['select_conseq'] == "gene_body":
            
                for block in vep_annotations.split(","):
                    _block = block.split("|")
                
                    if _block[hgvsc_index]:
                        return _block
        
                return vep_annotations.split(",")[0].split("|")
            
            else:
                i = 0
                offset = 5000000
                for _i, block in enumerate(vep_annotations.split(",")):
                    _block = block.split("|")
                    
                    # if gene body
                    if _block[hgvsc_index]:
                    
                        _aux = re.findall(r'c.+[+-](\d+)', _block[hgvsc_index].split(" ")[0])
             
                        if len(_aux) > 0:
                            _offset = int(_aux[0])
                         
                        else:
                            _offset = 0
                        
                        # keep the smallest offset
                        if _offset < offset:
                            offset = _offset
                            i = _i
                            
                return vep_annotations.split(",")[i].split("|")             
            
    else:
        return vep_annotations.split("|")
