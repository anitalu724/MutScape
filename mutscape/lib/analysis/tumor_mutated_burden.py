############################################################################################
# FileName     [ tumor_mutated_burden.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Calculate TMB for each data. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 9 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd

#########################################################
#                                                       #
#   python3 mafAnalysis.py \                            #
#   -f examples/test_data/maf/TCGA_test.maf \           #
#   -tmb 60456963 \                                     #
#   -o examples/output \                                #
#   -p examples/pic/                                    #
#                                                       #
#########################################################

class TumorMutationBurden:
    '''Mutation burden statistics
    Arguments:
        maf_file            {string}        -- The input MAF file for all data.
        output_folder       {string}        -- The path for output files.
        length              {int}           -- The length of genome (WES = 60456963)

    Parameters:
        self.head           {string}        -- The column names of MAF file.
        self.df             {pd.DataFrame}  -- The data for the MAF file.
    
    Outputs:
        TMB_analysis.tsv
        TMB_statistic.tsv
    ''' 
    def __init__(self, maf_file):
        print(colored(('\nStart Tumor Mutation Burden....'), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
    def data_analysis(self, output_folder, length):
        select_df = self.df[['Variant_Classification', 'Tumor_Sample_Barcode']]
        non = ['Missense_Mutation','Splice_Site', 'Translation_Start_Site','Nonstop_Mutation','Frame_Shift_Ins','Frame_Shift_Del','In_Frame_Del','In_Frame_Ins', 'Nonsense_Mutation']
        sample = select_df['Tumor_Sample_Barcode'].unique()
        sample_dict = {s:[0,0] for s in sample}     #list[all, nonsynonmous]
        for i in range(select_df.shape[0]):
            barcode, variant = select_df.iloc[i,:]['Tumor_Sample_Barcode'], select_df.iloc[i,:]['Variant_Classification']
            sample_dict[barcode][0] += 1
            if variant in non:
                sample_dict[barcode][1] += 1
            
        for s in sample_dict:
            sample_dict[s].append(sample_dict[s][0]*1000000/length)
            sample_dict[s].append(sample_dict[s][1]*1000000/length)
        d = pd.DataFrame(sample_dict).T.reset_index()
        d.columns = ['sample', 'All', 'nonsynonmous', 'TMB_All', 'TMB_nonsynonmous']
        stat_dict = {'mean' : [d['nonsynonmous'].sum()/len(sample), d['TMB_All'].sum()/len(sample), d['TMB_nonsynonmous'].sum()/len(sample)], 
                     'median' : [d['nonsynonmous'].median(), d['TMB_All'].median(), d['TMB_nonsynonmous'].median()], 
                     'max' : [d['nonsynonmous'].max(), d['TMB_All'].max(), d['TMB_nonsynonmous'].max()], 
                     'min' : [d['nonsynonmous'].min(), d['TMB_All'].min(), d['TMB_nonsynonmous'].min()],
                     'Q1' : [d['nonsynonmous'].quantile(q=0.25), d['TMB_All'].quantile(q=0.25), d['TMB_nonsynonmous'].quantile(q=0.25)],
                     'Q3' : [d['nonsynonmous'].quantile(q=0.75), d['TMB_All'].quantile(q=0.75), d['TMB_nonsynonmous'].quantile(q=0.75)]}
        stat_dict_df = pd.DataFrame(stat_dict).T
        stat_dict_df.columns = ['nonsynonmous', 'TMB_All', 'TMB_nonsynonmous']
        d.to_csv(output_folder+'TMB_analysis.tsv', sep='\t')
        stat_dict_df.to_csv(output_folder+'TMB_statistic.tsv', sep='\t')
        print(colored('=> Generate analysis files: ', 'green'))
        print(colored(('   '+output_folder+'TMB_analysis.tsv'), 'green'))
        print(colored(('   '+output_folder+'TMB_statistic.tsv'+'\n'), 'green'))
