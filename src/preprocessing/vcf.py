from collections import OrderedDict, defaultdict
from .utils import *
from cyvcf2 import VCF
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def process_vcf(vcf: str,
                tools_config: dict,
                thresholds: list,
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
    :param bool is_clinvar: Whether `vcf` is
        from clinvar database.
    :return dict: Dict with scores keyed by variant

    """
    logging.info("Extracting predictions from VCF.")

    # List of required INFO fields when
    # dataset is from Clinvar
    CLINVAR_FIELDS = ['CLNREVSTAT', 'CLNSIG']
    clinvar_confirmed = False

    # List of tools belonging to the
    # provided scope that will be analysed
    TOOLS_TO_ANALYSE = [t[0] for t in thresholds]

    # MAP with the list of VCF fields per tool
    TOOLS_CONFIG_MAP = {k: v[0] for k, v in tools_config.items() if k in TOOLS_TO_ANALYSE}

    # List with the VCF fields to check
    # Present fields in the header will be
    # popped out sequentially
    MISSING_VCF_FIELDS_FLAT = [i for sublist in [v for _k, v in TOOLS_CONFIG_MAP.items()]
                               for i in sublist]

    # Indexes is an ordered dict with the index(s)
    # where a information about a given tool is stored
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
    present_in_VEP_annot = []

    # Final dict where scores are stored
    scores = defaultdict(list)

    vcf_data = VCF(vcf)
    if vcf_data.contains("ANN") or vcf_data.contains("CSQ"):
        # Parse VCF header
        for field in vcf_data.header_iter():
            if field["HeaderType"] == "INFO" and (field["ID"] == "ANN" or field["ID"] == "CSQ"):
                tools_annotated_with_vep = field["Description"].split("Format:")[1][:-1].strip().split("|")

                # Looks only for scores belonging
                # to the specific scope set
                for _tool in TOOLS_TO_ANALYSE:

                    _tool_field = TOOLS_CONFIG_MAP[_tool]
                    _present, _absent = _check_if_field_exists(_tool_field, tools_annotated_with_vep)

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

        # Iterate over VCF
        keys = []
        for record in vcf_data:

            key = record.CHROM + "_" + str(record.POS) + "_" + str(record.ALT[0])
            if key in keys:
                logging.info("Variant {},{},{},{} is repeated in the VCF. "
                             "Keeping only first occurrence".format(record.CHROM,
                                                                    record.POS,
                                                                    record.REF,
                                                                    record.ALT[0]))
                continue

            else:
                keys.append(key)

            # MNPs
            if record.var_type == "indel" and len(record.REF) == len(record.ALT[0]):
                scores[key].extend([record.CHROM,
                                    record.POS,
                                    record.REF,
                                    record.ALT[0],
                                    record.ID,
                                    record.var_type, "mnp"])

            else:
                scores[key].extend([record.CHROM,
                                    record.POS,
                                    record.REF,
                                    record.ALT[0],
                                    record.ID,
                                    record.var_type,
                                    record.var_subtype])

            annotation = record.INFO.get("ANN")
            if annotation:
                info = annotation.split(",")[0].split("|")

                try:
                    scores[key].append(info[tools_annotated_with_vep.index('Existing_variation')])
                    scores[key].append(info[tools_annotated_with_vep.index('HGVSc')])
                    scores[key].append(info[tools_annotated_with_vep.index('SYMBOL')])
                    scores[key].append(info[tools_annotated_with_vep.index('Consequence')])
                    scores[key].append(info[tools_annotated_with_vep.index('gnomAD_AF')])
                    scores[key].append(("gnomAD_genomes", record.INFO.get("gnomADg_AF")))

                except ValueError:
                    pass

                # Look for externally annotated scores
                # that are single fields in the INFO field
                if len(present_in_INFO) > 0:

                    for _tool, _fields in present_in_INFO.items():
                        tool_out = []
                        for _f in _fields:
                            tool_out.append(record.INFO.get(_f))
                        scores[key].append((_tool, tool_out))

                # Append fields from VEP annotation
                for _tool, idx in vep_indexes.items():

                    # If tool has multiple fields
                    # and at least 1 is within VEP annot
                    # and other in the INFO annot.
                    if _tool in present_in_INFO.keys():
                        scores_so_far = scores[key]
                        _updated = []
                        for _v in scores_so_far:
                            if isinstance(_v, tuple) and _v[0] == _tool:
                                val = _v[1] + [info[i] for i in idx]
                                _updated.append((_v[0], val))
                            else:
                                _updated.append(_v)

                        scores[key] = _updated

                    else:
                        scores[key].append((_tool, [info[i] for i in idx]))

            if is_clinvar:
                scores[key].append(('CLNREVSTAT', record.INFO.get("CLNREVSTAT")))
                scores[key].append(('CLNSIG', record.INFO.get("CLNSIG")))

    else:
        raise ValueError("Program requires ANN/CSQ field to be present in the INFO field of the input VCF file.\n")

    vcf_data.close()
    return scores


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
