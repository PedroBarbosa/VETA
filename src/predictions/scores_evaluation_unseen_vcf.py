from preprocessing import *
from .apply_tools import apply_tool_predictions
from thresholds import threshold_list
from filters import filters
from plots.plots_unseen_vcf_analysis import *

def inspect_predictions(vcf, top_tools_file, n_top_tools, plot_tool, location, outdir):
    df = vep.get_df_ready(vcf, False, True, location)
    df = vep_cleaning(df)
    df = apply_tool_predictions(df, threshold_list).set_index('id')
    df_t = pd.concat(
        [df[[tool + "_prediction" for tool, *args in threshold_list]], df["HGVSc"], df["location"].to_frame()],
        axis=1).copy()
    df_t = df_t.rename(columns={col: col.split('_')[0] for col in df_t.columns})
    df_new = df_t.drop(["location", "HGVSc"], axis=1).apply(lambda x: x.value_counts(True, dropna=False),
                                                            axis=1).fillna(0).sort_values([True], ascending=False)
    df_new.rename(columns={False: 'is_benign', True: 'is_pathogenic', np.nan: "unpredictable"}, inplace=True)
    plot_area(df_new, outdir)
    df_new["unpredictable"] *= 100
    df_t.to_csv(os.path.join(outdir, "predictions.tsv"), sep="\t")
    plot_heatmap(df_new, outdir, False)
    plot_heatmap(df_new[df_new.is_pathogenic > 0.75], outdir, True)

    if top_tools_file and os.path.isfile(top_tools_file):
        tools = pd.read_csv(top_tools_file, sep="\t")['tool'].head(n_top_tools).tolist()
        tools.append("location")
        df_top = df_t[tools]
        booldic = {True: 1, False: -1}
        plot_heatmap_toptools(df_top.replace(booldic), filters, outdir)

    if plot_tool:
        for tool in plot_tool:
            plot_tool_score_distribution(df, tool, outdir)

