############################################################################################
# FileName     [ wgd_cin.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Whole-genome doubling (WGD) and Chromosome instability (CIN) analysis. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']

class WGDnCIN:
    """MAF analysis: Whole-genome doubling (WGD) and Chromosome instability (CIN)

    Parameters
    ----------
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.
    pic : str
        The path for storing plots.


    Output files
    ------------
    output :
        WGD_result.csv
        CIN_result.csv

    pictures:
        CIN_Score.pdf
        WGD_pie.pdf

    """
    def __init__(self, maf_file):
        print(colored(("\nStart analysing WGD and CIN...."), 'yellow'))
        self.list = ((pd.read_csv(maf_file, sep="\t"))[['CNV_input']].values.T)[0]
    def data_analysis(self, output_folder):
        genome_length = [249250621, 243199373, 198022430, 191154276, 180915260,
                         171115067, 159138663, 146364022, 141213431, 135534747,
                         135006516, 133851895, 115169878, 107349540, 102531392,
                         90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
        whole_length = sum(genome_length)
        WGD_list, CIN_list, sample_list = [], [], []
        for sample in self.list:
            df = pd.read_csv(sample, sep="\t")
            sample_name = pd.unique(df['SampleID'])[0]
            sample_list.append(sample_name)
            WGD_selected = df.loc[(df['A_cn'] >= 2)]
            CIN_selected = df.loc[(df['A_cn'] != 1) | (df['B_cn'] != 1)]
            WGD_SUM, CIN_SUM = 0, 0
            for i in range(WGD_selected.shape[0]):
                WGD_SUM += WGD_selected.iloc[i]['End_position']-WGD_selected.iloc[i]['Start_position']
            for i in range(CIN_selected.shape[0]):
                CIN_SUM += CIN_selected.iloc[i]['End_position']-CIN_selected.iloc[i]['Start_position']
            WGD_list.append(WGD_SUM >= 0.5*whole_length)
            CIN_list.append(CIN_SUM/whole_length)
        WGD_df = (pd.DataFrame([sample_list, WGD_list])).T
        CIN_df = (pd.DataFrame([sample_list, CIN_list])).T
        WGD_df.columns,CIN_df.columns = [['SampleID','WGD']], [['SampleID','CIN']]
        WGD_df.to_csv(output_folder + "WGD_result.csv",  index=False)
        CIN_df.to_csv(output_folder + "CIN_result.csv",  index=False)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + output_folder + "WGD_result.csv"), 'green'))
        print(colored(("   " + output_folder + "CIN_result.csv"), 'green'))
    def plotting(self, output_folder, pic):
        # WGD Pie Plot
        LABEL_SIZE, TITLE_SIZE = 24,30
        wgd_df = pd.read_csv(output_folder+"WGD_result.csv")
        data = [0,0]
        for i in range(wgd_df[['WGD']].shape[0]):
            if wgd_df[['WGD']].iloc[i]['WGD'] == False:
                data[1]+=1
            elif wgd_df[['WGD']].iloc[i]['WGD'] == True:
                data[0]+=1
        labels = 'WGD','Non-WGD'
        fig1, ax1 = plt.subplots()
        _, _, autotexts = ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[1],COLOR_MAP[0]] ,textprops={'fontsize': LABEL_SIZE}) #, 'size': LABEL_SIZE
        ax1.axis('equal')
        autotexts[0].set_color('black')
        autotexts[1].set_color('white')
        plt.savefig(pic+"WGD_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "WGD_pie.pdf"), 'green'))
        
        # CIN Bar plot
        CIN_df = pd.read_csv(output_folder+"CIN_result.csv")
        size = CIN_df.shape[0]
        CIN = tuple(list(CIN_df['CIN']))
        Sample = tuple(list(CIN_df['SampleID']))
        ind = np.arange(size)
        width = 0.7
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, CIN, width, color=COLOR_MAP[0])
        ax.set_ylabel('Scores', fontsize=LABEL_SIZE, fontweight='bold')
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.set_xlim([-1,len(ind)])
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_yticks(np.arange(0, 1, 0.2))
        ax.xaxis.set_visible(False)
        plt.savefig(pic+"CIN_Score.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "CIN_Score.pdf"), 'green'))
