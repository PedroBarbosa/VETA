from collections import OrderedDict, defaultdict
from .utils import *
from .location import *
from cyvcf2 import VCF
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

prediction_tools = ["ReMM", "Gwava", "ExACpLI", "dpsi_zscore", "ada_score", "rf_score",
               'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref', 'SpliceRegion',
               'GERP', 'phyloP', '29way_logOdds', 'phastCons',
               'MutationAssessor_score', 'SIFT_score', 'Polyphen2_HDIV_score',
               'Polyphen2_HVAR_score', 'MutationTaster_score', 'MutationTaster_pred',
               'MutationTaster_model', 'MutationTaster_AAE', 'FATHMM_score', 'fathmmMKL',
               'PROVEAN_score', 'PROVEAN_pred', 'CADD_PHRED',
               'DANN', 'Eigen', 'Eigen-PC', 'Eigen-pred_coding', 'Eigen-PC-phred_coding', 'MetaSVM_score',
               'MetaSVM_pred', 'MetaLR_score',
               'MetaLR_pred', 'REVEL_score', 'fitcons', 'LRT_score',
               'LRT_pred', 'MutPred_score', 'Condel', 'CAROL',
               'VEST4_score', 'M-CAP_score', 'SpliceAI', 'linsight_g', 'traP',
               'funseq2', 'HAL_DIFF', 'SCAP', 'kipoisplice_4', 'kipoisplice_4cons',
               'maxentscan_3', 'maxentscan_5', 'mmsplice_deltaPSI', 'mmsplice_efficiency', 'mmsplice_patho'
                    ]

tools_map = {'CADD_phred': 'CADD_PHRED',
               'GERP++_RS': 'GERP',
               'phyloP20way_mammalian': 'phyloP',  # new_is_100way
               'phastCons20way_mammalian': 'phastCons',  # new_is_100way
               '29way_logOdds': 'SiPhy',
               'integrated_fitCons_score': 'fitcons',
               'VEST3_score': 'VEST4_score',
               'fathmm-MKL_coding_score': 'fathmmMKL',
               'Eigen-PC-phred': 'Eigen-PC-phred_coding',
               'Eigen-phred': 'Eigen-pred_coding',
               'GWAVA': 'Gwava',
               }


def process_vcf_scores(vcf_file, is_clinvar=False):
    indexes, scores, absent = OrderedDict(), defaultdict(list), []
    vcf_data = VCF(vcf_file)
    if vcf_data.contains("ANN") or vcf_data.contains("CSQ"):
        spliceai_scores_single_field = False
        for field in vcf_data.header_iter():
            if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
                tools_annotated_with_vep = field["Description"].split("Format:")[1][:-1].strip().split("|")
                for tool in prediction_tools:
                    if tool in tools_annotated_with_vep:
                        indexes[tool] = tools_annotated_with_vep.index(tool)
                    else:
                        absent.append(tool)
            elif field["HeaderType"] == "INFO" and field["ID"] == "SpliceAI":
                spliceai_scores_single_field = True

        keys=[]
        for record in vcf_data:

            key = record.CHROM + "_" + str(record.POS) + "_" + str(record.ALT[0])
            if key in keys:
                print("Variant {},{},{},{} is repeated in the VCF. Keeping only first ocurrence".
                      format(record.CHROM, record.POS, record.REF, record.ALT[0]))
                continue
            else:
                keys.append(key)
            if record.var_type == "indel" and len(record.REF) == len(record.ALT[0]):
                scores[key].extend([record.CHROM, record.POS, record.REF, record.ALT[0], record.ID,
                                    record.var_type, "mnp"])
            else:
                scores[key].extend([record.CHROM, record.POS, record.REF, record.ALT[0], record.ID,
                                    record.var_type, record.var_subtype])
            annotation = record.INFO.get("ANN")
            if annotation:
                info = annotation.split(",")[0].split("|")

                try:
                    scores[key].append(info[tools_annotated_with_vep.index('Existing_variation')])
                    scores[key].append(info[tools_annotated_with_vep.index('HGVSc')])
                    scores[key].append(info[tools_annotated_with_vep.index('SYMBOL')])
                    scores[key].append(info[tools_annotated_with_vep.index('Consequence')])
                    scores[key].append(info[tools_annotated_with_vep.index('gnomAD_AF')])
                    #scores[key].append(info[tools_annotated_with_vep.index('gnomADg_AF')])

                except ValueError:
                    pass
                if len(absent) > 0:
                    scores[key].append(("gnomAD_genomes", record.INFO.get("gnomADg_AF")))
                    scores[key].append(("GERP", record.INFO.get("GERP")))
                    scores[key].append(("phyloP", record.INFO.get("phyloP")))
                    scores[key].append(("phastCons", record.INFO.get("phastCons")))
                    scores[key].append(("SiPhy", record.INFO.get("SiPhy")))
                    scores[key].append(("fitcons", record.INFO.get("fitcons")))
                    scores[key].append(("LINSIGHT", record.INFO.get("linsight_g")))
                    scores[key].append(("ReMM", record.INFO.get("ReMM")))
                    scores[key].append(("DANN", record.INFO.get("DANN")))
                    scores[key].append(("GWAVA", record.INFO.get("GWAVA")))
                    scores[key].append(("FATHMM-MKL", record.INFO.get("fathmmMKL")))
                    scores[key].append(("Eigen", record.INFO.get("Eigen")))
                    scores[key].append(("Eigen-PC", record.INFO.get("Eigen-PC")))
                    scores[key].append(("FunSeq2", record.INFO.get("funseq2")))
                    scores[key].append(("dpsi_zscore", record.INFO.get("dpsi_zscore")))
                    scores[key].append(("TraP", record.INFO.get("TraP")))
                    scores[key].append(("HAL", record.INFO.get("HAL_DIFF")))
                    scores[key].append(("S-CAP", record.INFO.get("SCAP")))
                    scores[key].append(("kipoiSplice4", record.INFO.get("kipoisplice_4")))
                    scores[key].append(("kipoiSplice4_cons", record.INFO.get("kipoisplice_4cons")))
                    scores[key].append(("kipoi_maxentscan_3", record.INFO.get("maxentscan_3")))
                    scores[key].append(("kipoi_maxentscan_5", record.INFO.get("maxentscan_5")))
                    scores[key].append(("mmsplice_deltaLogitPSI", record.INFO.get("mmsplice_deltaPSI")))
                    scores[key].append(("mmsplice_efficiency", record.INFO.get("mmsplice_efficiency")))
                    scores[key].append(("mmsplice_pathogenicity", record.INFO.get("mmsplice_patho")))
                    if spliceai_scores_single_field:
                        scores[key].append(("SpliceAI", record.INFO.get("SpliceAI")))
                        scores[key].append(("SpliceAI_ind", record.INFO.get("SpliceAI_ind")))
                    else:
                        scores[key].append(("SpliceAI", format_spliceai_fields(record,
                                                                    info[tools_annotated_with_vep.index('SYMBOL')])))

                for pl, i in indexes.items():
                    scores[key].append((pl, info[i]))

            if is_clinvar:
                scores[key].append(('CLNREVSTAT', record.INFO.get("CLNREVSTAT")))
                scores[key].append(('CLNSIG', record.INFO.get("CLNSIG")))

    else:
        logging.error("Program requires ANN field to be present in the INFO field of the input VCF file.\n")
        exit(1)

    vcf_data.close()
    return scores


def get_df_ready(vcf, isclinvar, isbenign, loc, deeper_intronic_analysis):
    logging.info("Processing {} file".format(os.path.basename(vcf)))
    scores = process_vcf_scores(vcf, isclinvar)
    df = pd.DataFrame.from_dict(scores, orient='index')

    new_col_names = ['hg19.chr', 'hg19.pos', 'ref', 'alt', 'id', 'type', 'subtype', 'rsID', 'HGVSc', 'Gene',
                     'Consequence', 'gnomAD_exomes']
    for column in df:
        if isinstance(df[column].iloc[0], (tuple,)):
            if df[column].iloc[0][0] == "gnomADg_AF":
                new_col_names.append("gnomAD_genomes")
            else:
                new_col_names.append(df[column].iloc[0][0])
            df[column] = df[column].map(tuple2float)

    rename_dict = {i: j for i, j in zip(list(df), new_col_names)}
    df.rename(columns=rename_dict, inplace=True)
    df['gnomAD_genomes'] = df['gnomAD_genomes'].replace(r'^\s*$', "0", regex=True)
    df['gnomAD_exomes'] = df['gnomAD_exomes'].replace(r'^\s*$', "0", regex=True)

    if loc == "HGVSc":
        hp = hgvs.parser.Parser()
        df['location'] = df['HGVSc'].apply(get_location, hp=hp)
        if deeper_intronic_analysis:
            df[['intron_bin', 'intron_offset']] = df.apply(lambda x: assign_intronic_bins(x['HGVSc'], hp, x['location']), axis=1)

    elif loc == "Consequence" and deeper_intronic_analysis:
        logging.error("If --intronic is set, --location must be HGVSc because intronic bin analysis will be performed"
                      " based on the HGVS nomenclature.")

    elif loc == "Consequence":
        df['location'] = df['Consequence'].apply(get_location_from_consequence)

    df['class'] = False if isbenign else True
    return df.reset_index()
