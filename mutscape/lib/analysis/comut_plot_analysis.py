############################################################################################
# FileName     [ comut_plot_analysis.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Calculate TMB for each data. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd

class CoMutAnalysis:
    '''MAF analysis: CoMut plot analysis

    Parameters
    ----------
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.

    Output files
    ------------
    output :
        mutation_data.tsv
        mutation_burden.tsv
    
    '''    
    def __init__(self, maf_file):
        print(colored(('\nStart CoMut_Plot_Analysis....'), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
    def data_analysis(self, output_folder):
        def mutation_type():
            maf = self.df
            chosen_col = maf[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification']]
            chosen_col = chosen_col.rename({'Tumor_Sample_Barcode':'sample', 'Hugo_Symbol':'category', 'Variant_Classification':'value'}, axis=1)
           
            value_old_list = ['Missense_Mutation', 'Nonsense_Mutation','In_Frame_Del', 'In_Frame_Ins','Splice_Site',
                              'Silent','Frame_Shift_Del','Frame_Shift_Ins','Nonstop_Mutation','Translation_Start_Site']
            remove_idx = []
            for i in range(len(chosen_col['value'])):
                if chosen_col['value'][i] not in value_old_list:
                    remove_idx.append(i)
                else:
                    if chosen_col['value'][i] == 'Missense_Mutation':
                        chosen_col['value'][i] = 'Missense'
                    elif chosen_col['value'][i] == 'Nonsense_Mutation':
                        chosen_col['value'][i] = 'Nonsense'
                    elif chosen_col['value'][i] == 'Nonstop_Mutation':
                        chosen_col['value'][i] = 'Nonstop'
                    elif chosen_col['value'][i] == 'In_Frame_Del' or chosen_col['value'][i] == 'In_Frame_Ins':
                        chosen_col['value'][i] = 'In frame indel'
                    elif chosen_col['value'][i] == 'Frame_Shift_Del' or chosen_col['value'][i] == 'Frame_Shift_Ins':
                        chosen_col['value'][i] = 'Frameshift indel'
                    elif chosen_col['value'][i] == 'Splice_Site':
                        chosen_col['value'][i] = 'Splice site'
                    elif chosen_col['value'][i] == 'Translation_Start_Site':
                        chosen_col['value'][i] = 'Translation start site'
            chosen_col = chosen_col.drop(chosen_col.index[remove_idx])
            # print(chosen_col)
            unique_chosen_col = chosen_col.drop_duplicates()
            # print(unique_chosen_col)
            # os._exit()
            unique_chosen_col.to_csv(output_folder+'mutation_data.tsv', sep = '\t', index = False)
            print(colored(('=> Generate CoMut_Analysis output files:'), 'green'))
            print(colored(('   '+output_folder+'mutation_data.tsv'), 'green'))
        def mutation_clonality():
            mutation_type_file = pd.read_csv(output_folder+'mutation_data.tsv', sep='\t', header=0)
            sample_dict = dict()
            for idx, data in mutation_type_file.iterrows():
                if data['sample'] not in sample_dict:
                    sample_dict[data['sample']] = [0, 0]
                if data['value'] == 'Silent':
                    sample_dict[data['sample']][1]+=1
                else:
                    sample_dict[data['sample']][0]+=1

            mutation_clone = pd.DataFrame.from_dict(sample_dict, orient='index')
            mutation_clone.reset_index(level=0, inplace=True)
            mutation_clone.to_csv(output_folder+'mutation_burden.tsv', sep='\t', header=['sample', 'Nonsynonymous', 'Synonymous'], index=False)
            print(colored(('   '+output_folder+'mutation_burden.tsv'+'\n'), 'green'))
        mutation_type()
        mutation_clonality() 
