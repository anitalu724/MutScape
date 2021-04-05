############################################################################################
# FileName     [ total_mutated_burden.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Calculate TMB for each data. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']

class HRDScore:
    """HRD score

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        ref             {string}    -- The reference genome used, grch38 or grch37 or mouse (default: grch38)
        pic             {string}    -- The path especially for output figures(.pdf)

    Outputs:
        all_HRDresults.csv

    Pictures:
        HRD_Score.pdf
        high_HRD_pie.pdf

    """
    def __init__(self, file):
        print(colored(("\nStart analysing HRD Score...."), 'yellow'))
        self.list = ((pd.read_csv(file, sep="\t"))[['CNV_input']].values.T)[0]
    def data_analysis(self, folder, ref):
        scar_r = open(folder + "scar.r", "a")
        scar_r.write("library(\"scarHRD\")\n")
        meta_list = []
        for i in self.list:
            scar_r.write("scar_score(\"" + i + "\", reference = \""+ref+"\", seqz = FALSE, outputdir = \"" + folder[:-1] + "\")\n")
        scar_r.close()
        os.system("Rscript " + folder + "scar.r\n")
        os.system("rm "+ folder + "scar.r\n")
        
        for file in os.listdir(folder):
            if file.endswith("_HRDresults.txt"):
                meta_list.append(file)
        meta_list.sort(reverse = True)
        final_df = pd.DataFrame()
        for meta in meta_list:
            df = pd.read_csv(folder+meta, sep="\t", index_col=False)
            final_df = pd.concat([df, final_df]) if not final_df.empty else df
            os.system("rm " + folder + meta + "\n")
        final_df.columns = [['Sample_id','HRD_LOH','Telomeric_AI','LST','HRD-sum']]
        final_df.to_csv(folder + "all_HRDresults.csv",  index=False)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + folder + "all_HRDresults.csv"), 'green'))
    def plotting(self, folder, pic):
        LABEL_SIZE, TITLE_SIZE = 24,30
        #Bar Plot
        df = pd.read_csv(folder+"all_HRDresults.csv")
        size = df.shape[0]
        HRD_LOH = tuple(list(df['HRD_LOH']))
        TAI = tuple(list(df['Telomeric_AI']))
        LST = tuple(list(df['LST']))
        SUM = list(df["HRD-sum"])
        Sample = tuple(list(df['Sample_id']))
        ind = np.arange(size)
        # print(ind)
        width = 0.7
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, HRD_LOH, width, color=COLOR_MAP[7])
        ax.bar(ind, TAI, width, bottom=HRD_LOH, color=COLOR_MAP[2])
        ax.bar(ind, LST, width, bottom=np.array(TAI)+np.array(HRD_LOH), color=COLOR_MAP[6])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.set_ylabel('Scores', fontsize=LABEL_SIZE, fontweight='bold')
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_ylim(top = max(SUM)*1.25)
        # ax.set_title('HRD Scores',fontsize=TITLE_SIZE, fontweight='bold')
        # plt.xticks(ind, Sample,rotation=45,horizontalalignment='right',fontweight='light', fontsize=12)
        ax.set_xlim([-1,len(ind)])
        ax.xaxis.set_visible(False)
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax.set_yticks(np.arange(0, max(SUM)*1.25+3, 10))
        ax.legend(labels=['HRD_LOH','Telomeric_AI','LST'], fontsize=LABEL_SIZE-4, edgecolor='white')
        plt.savefig(pic+"HRD_Score.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "HRD_Score.pdf"), 'green'))
        #Pie Plot
        over = len([i for i in SUM if i >=42])
        data = [over, size-over]
        labels='≧ 42','< 42'
        fig1, ax1 = plt.subplots()
        ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[7],COLOR_MAP[2]] ,textprops={'fontsize': LABEL_SIZE}) #'size': LABEL_SIZE
        ax1.axis('equal')
        # plt.title("Propotion of Samples with HRD Phenotype", fontsize=TITLE_SIZE, fontweight='bold')
        plt.savefig(pic+"high_HRD_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "high_HRD_pie.pdf"), 'green'))
