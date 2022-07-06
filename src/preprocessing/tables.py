import logging
import os
import sys
import pandas as pd
from typing import List

def generate_consequence_table(df: pd.DataFrame, out_dir:str):
    """
    Generates a tsv table with counts of 
    each variant consequence/location type
    
    :param pd.DataFrame df: Df with variants to be evaluated
    :param str out_dir: Output directory
    """
    out_dir = os.path.join(out_dir, "variant_counts")
    os.makedirs(out_dir, exist_ok=True)
    
    _counts = df[['location', 'outcome']].value_counts().reset_index().rename(columns={0: 'counts'})
    counts = _counts.pivot(index='location', columns='outcome', values='counts').fillna(0)

    if all(x in counts.columns for x in ['Pathogenic', 'Benign']):
        counts = counts[['Pathogenic', 'Benign']].astype(int).sort_values('Pathogenic', ascending=False)
    elif 'Benign' in counts.columns:
        counts['Pathogenic'] = 0
    else:
        counts['Benign'] = 0

    counts.to_csv(os.path.join(out_dir, "counts_per_consequence_type.tsv"), sep="\t")
    
    
def generate_clinvar_table(datasets: dict, filters_var_type: List, out_dir: str, clinvar_stars: str):
    """
    Generates a ready latex table with variant counts of
    each significance level

    :param dict datasets: Dictionary with one df for
        each clinvar significance level
    :param List filters_var_type: List with filters
        for each variant type
    :param str out_dir: Output directory to write
        the table
        
    :param str clinvar_stars: Level of clinvar stars used
    to filter the dataset. 
    """

    out_dir = os.path.join(out_dir, "variant_counts")
    os.makedirs(out_dir, exist_ok=True)
    datasets[clinvar_stars][['id', 'outcome']].to_csv(os.path.join(out_dir, 'clinvar_ids_used.tsv'), 
                                      header=False,
                                      index=False)
    
    for var_type, _func in filters_var_type:
        fname = os.path.join(out_dir, "{}_counts_per_clinvar_star_level.tex".format(var_type))
        fname_tsv = os.path.join(out_dir, "{}_counts_per_clinvar_star_level.tsv".format(var_type))
        
        if os.path.exists(fname):
            os.remove(fname)
        
        if os.path.exists(fname_tsv):
            os.remove(fname_tsv)

        stars = ["stars", "0*", "", "", "1*", "", "", "2*", "", "", "3*", "", "", "4*", "", ""]
        columns = ["counts", "B", "P", "T", "B", "P", "T", "B", "P", "T", "B", "P", "T", "B", "P", "T"]
        out_tsv = [stars, columns]

        df_per_var_type = {}

        for _name, df in datasets.items():
            df_per_var_type[_name] = _func(df).copy()

        locations = sorted(datasets['0s']['location'].unique())
        with open(fname, "w") as f:
            f.write("\\begin{tabular}{ | l | l | r r r | r r r | r r r | r r r | r r r | }\n")
            f.write("\\cline{3-14}")
            f.write("\\multicolumn{2}{c|}{} & \\multicolumn{3}{c|}{0*} & \\multicolumn{3}{c|}{1*} & \\multicolumn{3}{c|}{2*} & \\multicolumn{3}{"
                    "c|}{3*} & \\multicolumn{3}{c|}{4*} \\\\ \\cline{3-14}\n")
            cols = " & ".join(5 * ["B & P & T"])
            f.write(" \\multicolumn{2}{c|}{} & " + cols + " \\\\ \\hline\n")

            for i, likely in enumerate(["", "_l"]):

                if i == 0:
                    out_tsv.append(['Sure'])
                else:
                    out_tsv.append(['Sure and Likely'])

                for location in locations:
                    if location == locations[0]:
                        label = "\shortstack{Sure and \\\\ likely}" if likely == "_l" else "Sure"
                        f.write("\multirow{{{}}}{{*}}{{{}}}".format(len(locations), label))
                    else:
                        pass
                    f.write(" & {} ".format(location))

                    counts = [location]
                    for stars in [0, 1, 2, 3, 4]:
                        
                        dname = str(stars) + "s" + likely
                        df = df_per_var_type[dname]
 
                        sub_df = df[df['location'] == location]
                        count = sub_df.shape[0]
                        b = sub_df[sub_df['label'] == False].shape[0]
                        p = sub_df[sub_df['label'] == True].shape[0]
                        counts.extend([b, p, count])
                        f.write(" & {} & {} & {}".format(b, p, count))
                    if location == locations[-1]:
                        f.write("\\\\ \\hline \n")
                    else:
                        f.write("\\\\ \n")
                    out_tsv.append(counts)

            f.write("\\end{tabular}\n")

            with open(fname_tsv, "a") as out:
       
                for row in out_tsv:
                    out.write('\t'.join([str(x) for x in row]) + "\n")


def generate_proposed_thresholds_latex(thresholds: dict,
                                       location: str,
                                       out_dir: str):
    """
    Generates a ready latex table with new thresholds
    for each tool at each beta value

    :param dict thresholds: Dictionary with reference and
        proposed thresholds for all the tools
    :param str location: Location filter for the variants analysed

    :param str out_dir: Output directory to write
        the table
    """
    with open(os.path.join(out_dir, "proposed_thresholds_{}.tex".format(location)), 'w') as out:
        out.write("""\\begin{tabular}{ p{3cm} >{\\raggedleft\\arraybackslash}p{1.5cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm}}
    \\hline
    Tool       & Original & 1/1   & 1/2   & 1/3   & 1/5   & 1/10  \\\\
    \\hline
    """)
        for tool, new_t in thresholds.items():
            out.write("{}&\t{}\n".format(tool, '\t&'.join([str(t) for t in new_t])) + "\\\\\n")
        out.write("\\end{tabular}")
