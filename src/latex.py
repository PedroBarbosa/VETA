import os
from osutils import ensure_folder_exists
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

def generate_datasets_table(datasets, filters_var_type, OUT_DIR):

    if "3s_l" in datasets.keys(): #if clinvar data added
        logging.info("Generating clinvar dataset stats in latex..")
        for vartype, vartypefunction in filters_var_type:
            fname = os.path.join(OUT_DIR, "preprocessing", "generated_table_datasets_{}.tex".format(vartype))
            ensure_folder_exists(os.path.join(OUT_DIR, "preprocessing"))

            df_per_var_type = {}
            for dfname,df in datasets.items():
                df_per_var_type[dfname] = vartypefunction(df).copy()

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
                            b = sub_df[sub_df['class'] == False].shape[0]
                            p = sub_df[sub_df['class'] == True].shape[0]
                            f.write(" & {} & {} & {}".format(b, p, count))
                        if location == locations[-1]:
                            f.write("\\\\ \\hline \n")
                        else:
                            f.write("\\\\ \n")

                f.write("\\end{tabular}\n")
