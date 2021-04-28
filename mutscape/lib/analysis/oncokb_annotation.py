############################################################################################
# FileName     [ oncokb_annotation.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Actionable mutation(drug) annotation ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']

class OncoKBAnnotator:
    '''MAF analysis: Actionable mutation(drug) annotation

    Parameters
    ----------
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.
    path : str
        The path for oncokb-annotator folder.
    token : str
        The personal token provided from OncoKB.
    clinical : str
        The path for clinical data.
    pic : str
        The path for storing plots.
    cna : str (Optional) 
        The path for cna data
    level : str (Optional) 
        The level the user chooses (default = 4)

    Output files
    ------------
    output :
        maf_oncokb_output.txt
        clinical_oncokb_output.txt

    pictures:
        oncokb_total_pie.pdf
        oncokb_freq_actionable_genes.pdf

    '''
    def __init__(self, maf_file):
        print(colored(("\nStart OncoKB annotator(drug)...."), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
    def data_analysis(self, output_folder, path, token, clinical, cna = ''):
        selected_df = (self.df[['NCBI_Build','Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'HGVSp',  'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']]).set_index("Hugo_Symbol")
        selected_df.to_csv(output_folder + "maf_oncokb_input.txt", sep="\t")
        if not os.path.isdir('oncokb-annotator'):
            os.system("git clone https://github.com/oncokb/oncokb-annotator.git\n")
        os.system('cp lib/auxiliary/autoChange.py oncokb-annotator\n')
        os.chdir("oncokb-annotator")
        os.system('python3 autoChange.py\n')
        os.system('pip3 install requests\n')
        
        p = os.popen("python3 MafAnnotator.py -i ../"+output_folder + "maf_oncokb_input.txt -o ../" + output_folder + "maf_oncokb_output.txt -c ../"+clinical+" -b " + token + "\n")
        print(p.read())
        p.close()
        p = os.popen("python3 ClinicalDataAnnotator.py -i ../"+clinical+" -o ../"+ output_folder +"clinical_oncokb_output.txt -a ../"+output_folder+"maf_oncokb_output.txt\n")
        print(p.read())
        p.close()
        if cna !='':
            p = os.popen("python3 CnaAnnotator.py -i ../"+cna+" -o ../"+output_folder+"cna_oncokb_output.txt -c ../"+clinical+" -b "+ token + "\n")
            print(p.read())
            p.close()
        os.chdir("..")
        # os.system("rm -rf oncokb-annotator\n")
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + output_folder + "maf_oncokb_output.txt"), 'green'))
        print(colored(("   " + output_folder + "clinical_oncokb_output.txt"), 'green'))
    def plotting(self, output_folder, pic, level='4'):
        LABEL_SIZE, TITLE_SIZE = 24,30
        self.file = output_folder + "clinical_oncokb_output.txt"
        df = pd.read_csv(self.file, sep="\t")
        df_level = df[['HIGHEST_LEVEL']]
        level_list = ['LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4']
        level_dict = dict.fromkeys(level_list,0)
        sample_size = df.shape[0]
        for i in range(sample_size):
            if df_level.iloc[i]['HIGHEST_LEVEL'] in level_list:
                level_dict[df_level.iloc[i]['HIGHEST_LEVEL']] += 1

        true_num = 0
        if level == '4':
            true_num = sum(level_dict.values())
        elif level == '3':
            true_num = sum(level_dict.values()) - level_dict['LEVEL_4']
            level_list = ['LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B']
        elif level == '2':
            true_num = level_dict['LEVEL_1'] + level_dict['LEVEL_2']
            level_list = ['LEVEL_1', 'LEVEL_2']
        elif level == '1':
            true_num = level_dict['LEVEL_1']
            level_list = ['LEVEL_1']
        # Pie Plot( Total pie plot )
        size = [true_num, sample_size - true_num]
        labels = "Actionable\nbiomarkers"," Current absence\n of new actionable\n biomarkers"
        fig1, ax1 = plt.subplots()
        _, _, autotexts = ax1.pie(size, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[3],COLOR_MAP[4]] ,textprops={'fontsize': LABEL_SIZE})
        autotexts[1].set_color('white')
        ax1.axis('equal')
        plt.savefig(pic+"oncokb_total_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "oncokb_total_pie.pdf"), 'green'))
        
        # Bar Plot( Frequency of Actionable Genes )
        df_drug_count = df[level_list]
        drug_total_dict = {}
        for i in range(sample_size):
            drug_list = [[], [], [], [], [], []]    # [total, level1, level2, level3A, level3B, level4]
            for idx, item in enumerate(level_list):
                data = df_drug_count.iloc[i][item]
                if not pd.isna(data):
                    new_drug_list = data.split("(")
                    new_drug_list.pop(0)
                    new_drug_list = [drug.split(" ")[0] for drug in new_drug_list]
                    for j in new_drug_list:
                        if j not in drug_list[0]:
                            drug_list[0].append(j)
                            drug_list[idx+1].append(j)
            for j in range(5):
                for item in drug_list[j+1]:
                    if item not in drug_total_dict:
                        drug_total_dict[item] = [0,0,0,0,0]
                        drug_total_dict[item][j] += 1
                    else:
                        drug_total_dict[item][j] += 1
        
        drug_df = pd.DataFrame(drug_total_dict)
        ALL_DRUG_LIST = drug_df.columns
        
        SUM = (drug_df.sum()).tolist()

        List1 = tuple(list(np.asarray(drug_df.iloc[0,:].tolist())/sample_size))
        List2 = tuple(list(np.asarray(drug_df.iloc[1,:].tolist())/sample_size))
        List3A = tuple(list(np.asarray(drug_df.iloc[2,:].tolist())/sample_size))
        List3B = tuple(list(np.asarray(drug_df.iloc[3,:].tolist())/sample_size))
        List4 = tuple(list(np.asarray(drug_df.iloc[4,:].tolist())/sample_size))
        LEVEL_COLOR = ['#359744', '#286ea0', '#8a4a8d', '#b490b5','#3d3d3b']
        width = 0.7
        fig2 = plt.figure(figsize=(10,5))
        ax2 = fig2.add_axes([0,0,1,1])
        
        ax2.bar(ALL_DRUG_LIST, List1 ,  width, color=LEVEL_COLOR[0])
        ax2.bar(ALL_DRUG_LIST, List2 ,  width, bottom=List1,color=LEVEL_COLOR[1])
        ax2.bar(ALL_DRUG_LIST, List3A , width, bottom = np.array(List1)+np.array(List2),color=LEVEL_COLOR[2])
        ax2.bar(ALL_DRUG_LIST, List3B , width, bottom = np.array(List1)+np.array(List2)+np.array(List3A),color=LEVEL_COLOR[3])
        ax2.bar(ALL_DRUG_LIST, List4 ,  width, bottom = np.array(List1)+np.array(List2)+np.array(List3A)+np.array(List3B), color=LEVEL_COLOR[4])
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_color('#cac9c9')
        ax2.spines['left'].set_color('#cac9c9')
        plt.ylabel("Percentage", fontsize=LABEL_SIZE, fontweight='bold')
        ax2.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax2.tick_params(axis='y',direction='in', color='#cac9c9')
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
        ax2.set_yticks(np.arange(0, max(SUM)/sample_size*1.25, 0.2))
        plt.xticks(color='#222222',rotation=90, fontsize=LABEL_SIZE-4,fontstyle='italic',horizontalalignment='center',verticalalignment='top')#verticalalignment='bottom',
        ax2.legend(labels=level_list, fontsize=LABEL_SIZE-4, edgecolor='white')
        plt.savefig(pic+"oncokb_freq_actionable_genes.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "oncokb_freq_actionable_genes.pdf\n"), 'green'))
