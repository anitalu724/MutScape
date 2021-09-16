############################################################################################
# FileName     [ hrd_score.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Analysis of homologous recombination deficiency (HRD). ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

# from mutscape.lib.analysis.mutational_sig import LABEL_SIZE
from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']
LABEL_SIZE = 12

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



###############################################

# python3 mafAnalysis.py -f examples/test_data/maf/hrd.maf -hrdc examples/tsv/hrd_compare.tsv grch37 -o examples/output -p examples/pic/


class HRDCompare:
    '''
    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        ref             {string}    -- The reference genome used, grch38 or grch37 or mouse (default: grch38)
        pic             {string}    -- The path especially for output figures(.pdf)

    Parameters:
        self.type       {list}      -- [type1, type2, ...]
        self.fileList   {list}      -- [[file1-1, file1-2, ...], [file2-1, file2-2, ...]]
        self.outputFile {list}      -- [file1.csv, file2.csv, ...]

    Outputs:
        all_HRDresults.csv

    Pictures:
        HRD_Score.pdf
        high_HRD_pie.pdf

    '''
    

    def __init__(self, file):
        print(colored(("\nStart analysing HRD Comparing...."), 'yellow'))
        df = (pd.read_csv(file, sep='\t', index_col=None)).dropna(axis='columns')
        self.type = list(df.columns)
        self.fileList = [list(df[i]) for i in self.type]
        self.hrdFile, self.wgdFile, self.cinFile = [], [], []
        
    def HRD(self, idx, fileList, output_folder, ref):
        scar_r = open(output_folder + "scar.r", "a")
        scar_r.write("library(\"scarHRD\")\n")
        meta_list, delete_list, sample_list = [], [], []
        
        for i in fileList:
            # check if chrx,y exists or total=2&A_cn=1&B_cn=1
            tmp_df = pd.read_csv(i, sep = '\t')
            sample_list.append(tmp_df['SampleID'][0])
            
            chrx = tmp_df.loc[tmp_df['Chromosome'] == 'chrX']
            chry = tmp_df.loc[tmp_df['Chromosome'] == 'chrY']
            total2 = tmp_df.loc[tmp_df['total_cn'] == 2].loc[tmp_df['A_cn'] == 1].loc[tmp_df['B_cn'] == 1]
            delete = pd.concat([chrx, chry, total2]).drop_duplicates().reset_index(drop=True)
            if delete.shape[0] != 0:
                tmp_df = tmp_df.drop(chrx.index).drop(chry.index).drop(total2.index)
            
            if tmp_df.shape[0] != 0:
                scar_r.write("scar_score(\"" + i + "\", reference = \""+ref+"\", seqz = FALSE, outputdir = \"" + output_folder[:-1] + "\")\n")
            else:
                delete_list.append(i)
        scar_r.close()
        
        os.system("Rscript " + output_folder + "scar.r\n")

        os.system("rm "+ output_folder + "scar.r\n")
        
        for file in os.listdir(output_folder):
            if file.endswith("_HRDresults.txt"):
                meta_list.append(file)
        meta_list.sort(reverse = True)
        final_df = pd.DataFrame()

        for sampleID in sample_list:
            if sampleID+'_HRDresults.txt' in meta_list:
                df = pd.read_csv(output_folder+sampleID+'_HRDresults.txt', sep="\t", index_col=False)
                final_df = pd.concat([final_df, df]) if not final_df.empty else df
            else:
                
                new_list = [sampleID, 0, 0, 0, 0]
                new_df = pd.DataFrame(new_list).T
                new_df.columns = final_df.columns
                final_df = pd.concat([final_df, new_df]) if not final_df.empty else new_df
                
        final_df.columns = [['Sample_id','HRD_LOH','Telomeric_AI','LST','HRD-sum']]
        output_file = output_folder + 'all_HRDresults_' + self.type[idx] + '.csv'
        final_df.to_csv(output_file,  index=False)
        self.hrdFile.append(output_file)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + output_file), 'green'))

    def WGD_CIN(self, idx, fileList, output_folder):
        genome_length = [249250621, 243199373, 198022430, 191154276, 180915260,
                         171115067, 159138663, 146364022, 141213431, 135534747,
                         135006516, 133851895, 115169878, 107349540, 102531392,
                         90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
        whole_length = sum(genome_length)
        WGD_list, CIN_list, sample_list = [], [], []
        for sample in fileList:
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
        WGD_df, CIN_df = (pd.DataFrame([sample_list, WGD_list])).T, (pd.DataFrame([sample_list, CIN_list])).T
        wgd_file, cin_file = output_folder + 'WGD_result_'+self.type[idx]+'.csv', output_folder + 'CIN_result_'+self.type[idx]+'.csv'
        WGD_df.columns, CIN_df.columns = [['SampleID','WGD']], [['SampleID','CIN']]
        WGD_df.to_csv(wgd_file, index=False)
        CIN_df.to_csv(cin_file, index=False)
        self.wgdFile.append(wgd_file)
        self.cinFile.append(cin_file)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + output_folder + 'WGD_result_'+self.type[idx]+'.csv'), 'green'))
        print(colored(("   " + output_folder + 'CIN_result_'+self.type[idx]+'.csv'), 'green'))

    def WGDheatmap(self, pic):
        wgdList = []
        for wgd_file in self.wgdFile:
            wgdList.append([int(elem) for elem in list(pd.read_csv(wgd_file)['WGD'])])
        wgdList.reverse()
        yLabel = self.type
        yLabel.reverse()

        M = np.array(wgdList)
        sns.set(font_scale=2)
        sns.set_style('white')
        
        f, ax = plt.subplots(1, 1, figsize=(20,6))
        ax = sns.heatmap(M, vmin=0, vmax = 1, square = True, yticklabels = yLabel, xticklabels = False, linewidth = 1, cmap=['#b7d5ea','#266199'], ax = ax, cbar_kws={'orientation': 'horizontal','shrink':1, 'aspect':70})
        
        colorbar = ax.collections[0].colorbar 
        r = M.max().max()
        colorbar.set_ticks([0.25*r, 0.75*r])
        colorbar.set_ticklabels(['Non-WGD' , 'WGD'])                       
        colorbar.ax.tick_params(labelsize = LABEL_SIZE + 12)    
        
        ax.tick_params(axis='both',length=0)
        ax.set_yticklabels(ax.get_yticklabels(), color='#222222', rotation = 'horizontal', fontsize=LABEL_SIZE + 14, fontweight = 'bold')
        plt.ylim(bottom=0, top=len(wgdList)+0.5)
        plt.savefig(pic+'WGD_heatmap.pdf',dpi = 300,bbox_inches='tight')
        plt.close()
        print(colored(('=> Generate WGD comparison Plot: '+pic+'WGD_heatmap.pdf'), 'green'))  
    
    def HRDheatmap(self, pic):
        hrdList = []
        for hrd_file in self.hrdFile:
            hrdList.append([int(elem >= 42) for elem in list(pd.read_csv(hrd_file)['HRD-sum'])])
        hrdList.reverse()
        yLabel = self.type
        # yLabel.reverse()

        M = np.array(hrdList)
        sns.set(font_scale=2)
        sns.set_style('white')
        
        f, ax = plt.subplots(1, 1, figsize=(20,6))
        
        ax = sns.heatmap(M, vmin=0, vmax = 1, square = True, yticklabels = yLabel, xticklabels = False, linewidth = 1, cmap=[ '#E0CADB',  '#695D73'], ax = ax, cbar_kws={'orientation': 'horizontal','shrink':1, 'aspect':70})
        colorbar = ax.collections[0].colorbar 
        
        r = M.max().max()
        colorbar.set_ticks([0.25*r, 0.75*r])
        colorbar.set_ticklabels(['HRD<42' , 'HRD≧42'])   
        colorbar.ax.tick_params(labelsize = LABEL_SIZE + 12)                        
        
        ax.tick_params(axis='both',length=0)
        ax.set_yticklabels(ax.get_yticklabels() , color='#222222', rotation = 'horizontal', fontsize=LABEL_SIZE + 14, fontweight = 'bold')
        plt.ylim(bottom=0, top=len(hrdList)+0.5)
        plt.savefig(pic+'HRD_heatmap.pdf',dpi=300,bbox_inches='tight')
        plt.close()
        print(colored(('=> Generate HRD comparison Plot: '+pic+'HRD_heatmap.pdf'), 'green'))  

    def CINbarplot(self, pic):
        CIN_COLOR_MAP = ['#266199','#b7d5ea']
        cinList = []
        for cin_file in self.cinFile:
            cinList.append(list(pd.read_csv(cin_file)['CIN']))
            
        col_num = len(cinList[0])
        
        for df in cinList:
            if len(df) != col_num:
                raise ValueError('[MutScape] The files\' format are not available.')

        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_axes([0,0,1,1])
        bar_width = 0.35
        index = np.arange(col_num)

        barList = []
        for idx, df in enumerate(cinList):
            tmp_bar = ax.bar(index+idx*0.35, tuple(df), width = bar_width, label = self.type[idx], color = CIN_COLOR_MAP[idx], linewidth = 0)
            barList.append(tmp_bar)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_xlim([-0.5,len(index)])
        ax.set_ylabel('CIN Score', fontsize=LABEL_SIZE, fontweight='bold')
        plt.yticks(fontsize=LABEL_SIZE - 2)
        ax.set_yticks(np.arange(0, 1, 0.2))
        ax.xaxis.set_visible(False)
        
        plt.legend()
        plt.savefig(pic+'CIN_barplot.pdf',dpi = 300, bbox_inches='tight')
        print(colored(('=> Generate CIN comparison Bar Plot: '+ pic +'CIN_barplot.pdf'), 'green'))  


    def HRDbarplot(self, pic):
        HRD_COLOR_MAP = ['#031A54','#266199','#b7d5ea', '#561729', '#C88984', '#EFCAC2', '#5D3E6B', '#695D73', '#E0CADB', '#71a0a5','#9Ab377','#ACd5AA']
        
        hrdList = []
        for hrd_file in self.hrdFile:
            tmp, df  = [], pd.read_csv(hrd_file)
            tmp.append(list(df['HRD_LOH']))
            tmp.append(list(df['Telomeric_AI']))
            tmp.append(list(df['LST']))
            tmp.append(list(df['HRD-sum']))
            hrdList.append(tmp)
        
        #Bar Plot
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_axes([0,0,1,1])
        for idx, type in enumerate(hrdList):
            size = len(type[0])
            HRD_LOH = tuple(hrdList[idx][0])
            TAI = tuple(hrdList[idx][1])
            LST = tuple(hrdList[idx][2])
            SUM = list(hrdList[idx][3])
        
            index = np.arange(size)
        
            width = 0.35
            ax.bar(index+idx*0.35, HRD_LOH, width, color=HRD_COLOR_MAP[0+idx*3])
            ax.bar(index+idx*0.35, TAI, width, bottom=HRD_LOH, color=HRD_COLOR_MAP[1+idx*3])
            ax.bar(index+idx*0.35, LST, width, bottom=np.array(TAI)+np.array(HRD_LOH), color=HRD_COLOR_MAP[2+idx*3])
        
        legend_list = []
        for i in range(3):
            tmp_list = []
            tmp_plot1,  = ax.plot([], [], c = HRD_COLOR_MAP[i+0] , marker='s', markersize=10, fillstyle='left', linestyle='none', mec = 'None')
            tmp_plot2,  = ax.plot([], [], c = HRD_COLOR_MAP[i+3] , marker='s', markersize=10, fillstyle='right', linestyle='none', mec = 'None')
            tmp_list.append(tmp_plot1)
            tmp_list.append(tmp_plot2)
            legend_list.append(tuple(tmp_list))
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.set_ylabel('HRD Score', fontsize=LABEL_SIZE, fontweight='bold')
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_ylim(top = max(SUM)*1.25)
        
        ax.set_xlim([-0.5,len(index)])
        ax.xaxis.set_visible(False)
        plt.yticks(fontsize=LABEL_SIZE - 2)
        ax.set_yticks(np.arange(0, max(SUM)*1.25+3, 10))
        ax.legend(tuple(legend_list), ('HRD_LOH','Telomeric_AI','LST'), labelspacing=0.5, loc='upper right', fontsize=LABEL_SIZE-2, edgecolor='white')

        plt.savefig(pic+"HRD_barplot.pdf", dpi = 300, bbox_inches='tight')
        print(colored(("=> Generate HRD Compare Bar Plot: " + pic + "HRD_barplot.pdf"), 'green'))


        