import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import gzip
import sys
import logging
import process_fromVCF
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
#plt.switch_backend('TkAgg')

def process_db(df,db):
    if db == "OMIM":
        db_ids = df["CLNDISDB"].str.extractall("(OMIM:[0-9]+)")
        return set(db_ids[0].str.replace(db + ":", "").tolist())
    else:
        toremove = {'CN169374', 'CN517202'}
        medgen = df["CLNDISDB"].str.extractall("(MedGen:[0-9A-Za-z]+)")
        s = set(medgen[0].str.replace("MedGen:", "").tolist())
        return s - toremove

def get_disease_df(db,df,id):
    df_=df[(df["CLNDISDB"].str.contains(":" + id + "[,|]")) | (df["CLNDISDB"].str.contains(":" + id + "$"))].copy()
    genes=df_['GENEINFO'].str.split("|").dropna().tolist() #dropna to remove some variants with no GENEINFO field: e.g  some non coding transcripts
    unique_genes = set([gene for sublist in genes for gene in sublist])
    npatho=np.sum(df_['class'])
    return df_, unique_genes,npatho

def get_disease_name(db,db_dir):
    fname = [filename for filename in os.listdir(db_dir) if db in filename]
    if fname:
        #dname_dic = pd.read_csv(os.path.join(db_dir, fname[0]),index_col=0,usecols=[0,1],header=0,compression="gzip",error_bad_lines=False).to_dict(orient='index')
        dict={}
        with gzip.open(os.path.join(db_dir, fname[0]),'rb') as fin:
            for line in fin:
                l=line.decode('utf-8')
                dict[l.split('\t')[0]] = l.split("\t")[1].rstrip()
        return dict
    else:
        print("{} file doesn't exist. IDs will be used.")
        return {}

def plot_db_density(outdir,dic,db,star):
    sns.set()
    f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(7,7))

    nvariants=[i[0] for i in dic.values()]
    ngenes=[i[1] for i in dic.values()]
    ratio=[i[2]/i[0]*100 for i in dic.values()]
    ratio_morethan_100 = [i[2]/i[0]*100 for i in dic.values() if i[0] > 100]

    sns.distplot(nvariants, color="red", ax=ax1)
    sns.distplot(ngenes, color="olive", ax=ax2)
    sns.distplot(ratio, bins=10, color="blue", ax=ax3)
    sns.distplot(ratio_morethan_100, bins=10, color="blue", ax=ax4)
    ax1.set_xlabel("Variants per {} id".format(db))
    ax1.set_ylabel("Frequency")
    #ax1.xaxis.set_ticks(np.arange(0, max(nvariants), 20))
    ax2.set_xlabel("Genes per {} id".format(db))
    ax3.set_xlabel("% of pathogenic variants")
    ax3.set_xlim(0,100)
    ax4.set_xlabel("% of pathogenic variants (>100 variants)")
    ax4.set_xlim(0,100)
    plt.savefig(os.path.join(outdir, "hist_{}star_{}.pdf".format(star,db)))
    plt.close()

    plt.figure()
    plt.xlabel("#Variants")
    plt.ylabel("#genes")
    plt.title("{} correlation".format(db))
    plt.scatter(nvariants,ngenes,color="dimgrey", alpha=0.8)
    plt.savefig(os.path.join(outdir,"scatter_{}star_{}.pdf".format(star,db)))
    plt.close()


def plot_50(outdir,ids_info,db,star):
    dname = [i[1][3] if i[1][3] != None else i[0] for i in ids_info]
    nvariants =  [i[1][0] for i in ids_info]
    y = [x + "_N=" + str(n) for x, n in zip(dname, nvariants)]
    patho_ratio = [i[1][2]/i[1][0]*100 for i in ids_info]
    npathogenic = [i[1][2] for i in ids_info]
    nbenign = [ntotal - npatho for ntotal, npatho in zip(nvariants, npathogenic)]

    sns.set(style="whitegrid")
    plt.figure(figsize=(11,8))
    sns.stripplot(x=patho_ratio,y=y,color="brown",size=8,orient="h", linewidth=1, edgecolor="w")
    ax = plt.axes()
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.set_xlabel("% of pathogenic variants")
    plt.tight_layout()
    #plt.title("{} balance of variant classes in the top 50 phenotypes with more variants associated".format(db))
    plt.savefig(os.path.join(outdir,"top50_strip_{}star_{}.pdf".format(star,db)))
    plt.close()

    fig,axes = plt.subplots(ncols=2,sharey=True,figsize=(15,10))
    axes[0].barh(dname, nbenign, align='center', color='navy')
    axes[1].barh(dname, npathogenic, align='center', color='firebrick')
    range=max(npathogenic + nbenign) + 2

    axes[0].set_xlim(0,range)
    axes[1].set_xlim(0,range)
    axes[0].invert_xaxis()
    axes[0].set_xlabel('#Benign variants')
    axes[1].set_xlabel('#Pathogenic variants')
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "top50_bar_{}star_{}.pdf".format(star, db)))

    #fig.suptitle("{} balance of benign/pathogenic variants in the top 50 phenotypes".format(db))
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.close()

def main(df,outdir):
    clinvar_stars = {1: process_fromVCF.filter_clinvar_1_star,
                     2: process_fromVCF.filter_clinvar_two_stars,
                     3: process_fromVCF.filter_clinvar_three_stars}
    dbs = ['OMIM', 'MedGen']
    dbname_dir = os.path.dirname(os.path.abspath(sys.argv[0]).replace("analysis", "dbs"))
    for db in dbs:
        logging.info("Started {} analysis by retrieving disease names.".format(db))
        dname_dict = get_disease_name(db, dbname_dir)
        for i in range(1, 4):
            out_dir = os.path.join(outdir, str(i) + "star")
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            logging.info("Clinvar {} star(s)..".format(i))
            df_ = clinvar_stars[i](df)
            logging.info("{} {} star(s) variants found".format(df_.shape[0], i))
            ids = process_db(df_, db)
            logging.info("{} different {} ids found".format(len(ids), db))
            logging.info("Phenotype specific analysis started.")
            density_dic, id_idx, old = {}, 1, -1
            for id in ids:
                percentage = int(id_idx / len(ids) * 100)
                if percentage in range(0, 100, 5) and percentage != old:
                    old = percentage
                    logging.info("Progress: {}%".format(percentage))
                df_id, genes, npatho = get_disease_df(db, df_, id)
                try:
                    dname = dname_dict[id]
                except KeyError:
                    dname = None

                density_dic[id] = (df_id.shape[0], len(genes), npatho, dname)
                id_idx += 1
                if np.sum(df_id['class']) >= 50 and np.sum(~df_id['class']) >= 50:
                    logging.info("{} {} id will be used for machine learning analysis".format(id, db))
                    out_dir_id = os.path.join(out_dir, db + "_" + id)
                    os.mkdir(out_dir_id)
                    os.mkdir(os.path.join(out_dir_id, "cv"))
                    os.mkdir(os.path.join(out_dir_id, "standard"))
                    # machine_learning_analysis(df_id, outdir=os.path.join(out_dir_id, "cv"))
                    # machine_learning_analysis(df_2s, df_test=df_id, outdir=os.path.join(out_dir_id, "standard"), robustness=True)

            plot_db_density(out_dir, density_dic, db, i)
            ids_with_more_variants = sorted(density_dic.items(), key=lambda x: x[1][0], reverse=True)[:50]
            plot_50(out_dir, ids_with_more_variants, db, i)