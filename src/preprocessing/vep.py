from collections import OrderedDict, defaultdict
from .utils import *
from .location import *
from cyvcf2 import VCF


vep_plugins = ["ReMM", "Gwava", "ExACpLI", "dpsi_zscore", "ada_score", "rf_score",
               'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref', 'SpliceRegion',
               'GERP', 'phyloP', '29way_logOdds', 'phastCons',
               'MutationAssessor_score', 'SIFT_score', 'Polyphen2_HDIV_score',
               'Polyphen2_HVAR_score', 'MutationTaster_score', 'MutationTaster_pred',
               'MutationTaster_model', 'MutationTaster_AAE', 'FATHMM_score', 'fathmmMKL',
               'PROVEAN_score', 'PROVEAN_pred', 'CADD_PHRED',
               'DANN', 'Eigen', 'Eigen-PC', 'Eigen-pred_coding', 'Eigen-PC-phred_coding', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score',
               'MetaLR_pred', 'REVEL_score', 'fitcons', 'LRT_score',
               'LRT_pred', 'MutPred_score', 'Condel', 'CAROL',
               'VEST4_score','M-CAP_score','SpliceAI', 'linsight_g', 'traP',
               'funseq2', ]

plugins_map = {'CADD_phred':'CADD_PHRED',
               'GERP++_RS':'GERP',
               'phyloP20way_mammalian': 'phyloP', #new_is_100way
               'phastCons20way_mammalian':'phastCons', #new_is_100way
               '29way_logOdds': 'SiPhy',
               'integrated_fitCons_score' : 'fitcons',
                'VEST3_score': 'VEST4_score',
                'fathmm-MKL_coding_score': 'fathmmMKL',
               'Eigen-PC-phred':'Eigen-PC-phred_coding',
                'Eigen-phred': 'Eigen-pred_coding',
               'GWAVA': 'Gwava',
                }


def process_vep_scores(vcf_file, is_clinvar=False):
    indexes, scores, absent = OrderedDict(), defaultdict(list), []
    vcf_data = VCF(vcf_file)
    if vcf_data.contains("ANN") or vcf_data.contains("CSQ"):
        spliceai_scores_fromtool=False
        for field in vcf_data.header_iter():
            if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
                tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
                for plugin in vep_plugins:
                    if plugin in tools:
                        indexes[plugin] = tools.index(plugin)
                    else:
                        absent.append(plugin)
            elif field["HeaderType"] == "INFO" and field["ID"] == "SpliceAI":
                spliceai_scores_fromtool=True

        for record in vcf_data:

            key = record.ID + "_" + str(record.POS)
            if record.var_type == "indel" and len(record.REF) == len(record.ALT[0]):
                scores[key].extend([record.CHROM,record.POS, record.REF, record.ALT[0], record.ID,
                                    record.var_type, "mnp"])
            else:
                scores[key].extend([record.CHROM, record.POS, record.REF, record.ALT[0], record.ID,
                                    record.var_type, record.var_subtype])
            annotation = record.INFO.get("ANN")
            if annotation:
                info = annotation.split(",")[0].split("|")

                try:
                    scores[key].append(info[tools.index('Existing_variation')])
                    scores[key].append(info[tools.index('HGVSc')])
                    scores[key].append(info[tools.index('SYMBOL')])
                    scores[key].append(info[tools.index('Consequence')])
                except ValueError:
                    pass
                if len(absent) > 0:
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
                    scores[key].append(("funseq2", record.INFO.get("funseq2")))
                    scores[key].append(("dpsi_zscore", record.INFO.get("dpsi_zscore")))
                    scores[key].append(("traP", record.INFO.get("traP")))
                    if spliceai_scores_fromtool:
                        scores[key].append(("SpliceAI", record.INFO.get("SpliceAI")))
                    else:
                        scores[key].append(("SpliceAI", format_spliceai_fields(record, info[tools.index('SYMBOL')])))

                for pl, i in indexes.items():
                    scores[key].append((pl, info[i]))

            if is_clinvar:
                scores[key].append(('CLNREVSTAT', record.INFO.get("CLNREVSTAT")))
                scores[key].append(('CLNSIG', record.INFO.get("CLNSIG")))


    else:
        print("Program requires ANN field to be present in the INFO field of the input VCF file.\n")
        exit(1)

    vcf_data.close()
    return scores


def get_df_ready(vcf,isclinvar,isbenign,loc):
    scores = process_vep_scores(vcf,isclinvar)
    df = pd.DataFrame.from_dict(scores, orient='index')
    new_col_names = ['hg19.chr', 'hg19.pos', 'ref', 'alt', 'id','type', 'subtype', 'rsID','HGVSc', 'Gene','Consequence']
    for column in df:
        if isinstance(df[column].iloc[0], (tuple,)):
            new_col_names.append(df[column].iloc[0][0])
            df[column] = df[column].map(tuple2float)

    rename_dict = {i: j for i, j in zip(list(df), new_col_names)}
    df.rename(columns=rename_dict, inplace=True)

    if loc == "HGVSc":
        hp = hgvs.parser.Parser()
        df['location'] = df['HGVSc'].apply(get_location, hp=hp)

    else:
        df['location'] = df['Consequence'].apply(get_location_from_consequence)

    df['class'] = False if isbenign else True
    return df.reset_index()
