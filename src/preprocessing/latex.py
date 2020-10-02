import os
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
from typing import List


def generate_clinvar_table(datasets: dict, filters_var_type: List, out_dir: str):
    """
    Generates a ready latex table with variant counts of
    each significance level

    :param dict datasets: Dictionary with one df for
        each clinvar significance level
    :param List filters_var_type: List with filters
        for each variant type
    :param str out_dir: Output directory to write
        the table
    """

    out_dir = os.path.join(out_dir, "clinvar_stats")
    os.makedirs(out_dir)

    for var_type, _func in filters_var_type:
        fname = os.path.join(out_dir, "{}_table_datasets.tex".format(var_type))

        df_per_var_type = {}
        for _name, df in datasets.items():
            df_per_var_type[_name] = _func(df).copy()

        locations = sorted(datasets['clinvar']['location'].unique())
        with open(fname, "w") as f:
            f.write("\\begin{tabular}{ | l | l | r r r | r r r | r r r | r  r  r | }\n")
            f.write("\\cline{3-14}")
            f.write(" \\multicolumn{2}{c|}{} & \\multicolumn{3}{c|}{1*} & \\multicolumn{3}{c|}{2*} & \\multicolumn{3}{c|}{3*} & \\multicolumn{3}{c|}{4*} \\\\ \\cline{3-14}\n")
            cols = " & ".join(4*["B & P & T"])
            f.write(" \\multicolumn{2}{c|}{} & " + cols + " \\\\ \\hline\n")

            for likely in ["", "_l"]:

                for location in locations:
                    if location == locations[0]:
                        label = "\shortstack{Sure and \\\\ likely}" if likely == "_l" else "Sure"
                        f.write("\multirow{{{}}}{{*}}{{{}}}".format(len(locations), label))
                    else:
                        # Do nothing
                        pass
                    f.write(" & {} ".format(location))
                    for stars in [1,2,3,4]:
                        dname = str(stars) + "s" + likely
                        df = df_per_var_type[dname]
                        sub_df = df[df['location'] == location]
                        count = sub_df.shape[0]
                        b = sub_df[sub_df['label'] == False].shape[0]
                        p = sub_df[sub_df['label'] == True].shape[0]
                        f.write(" & {} & {} & {}".format(b, p, count))
                    if location == locations[-1]:
                        f.write("\\\\ \\hline \n")
                    else:
                        f.write("\\\\ \n")

            f.write("\\end{tabular}\n")
