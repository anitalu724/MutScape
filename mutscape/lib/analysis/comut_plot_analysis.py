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

from comut import comut
from comut import fileparsers

class CoMutPlot:
    def __init__(self, file, info):
        print(colored(("\nStart plotting CoMut Plot...."), 'yellow'))
        self.fd = (pd.read_csv(file, sep="\t")).to_dict('list')
        for item in self.fd:
            cleanedList = [x for x in self.fd[item] if str(x) != 'nan']
            self.fd[item] = cleanedList
        self.info = (pd.read_csv(info, sep="\t")).to_dict('list')
        for item in self.info:
            cleanedList = [x for x in self.info[item] if str(x) != 'nan']
            self.info[item] = cleanedList
    def plot(self, pic, theme, name):
        fixed_category = ['Same Patient','Copy Number Alteration','Mutation Type','Purity','Mutation Signature','Mutation Classification','Frequency', 'Whole Genome Doubling']
        feature_category = [item for item in self.fd if item not in fixed_category]
        fixed_path_list = [self.fd[item][0] if len(self.fd[item]) != 0 else "" for item in fixed_category ]
        
        feature_path_list = [self.fd[item][0] for item in feature_category ]
      
        print(colored(("*   Reading Data...."), 'yellow'))
        same_patient_df = pd.read_csv(fixed_path_list[0], sep = '\t') if fixed_path_list[0] != "" else pd.DataFrame(data={})
        copy_num_df     = pd.read_csv(fixed_path_list[1], sep = '\t') if fixed_path_list[1] != "" else pd.DataFrame(data={})
        mut_type_df     = pd.read_csv(fixed_path_list[2], sep = '\t') if fixed_path_list[2] != "" else pd.DataFrame(data={})
        purity_df       = pd.read_csv(fixed_path_list[3], sep = '\t') if fixed_path_list[3] != "" else pd.DataFrame(data={})
        mut_sig_df      = pd.read_csv(fixed_path_list[4], sep = '\t') if fixed_path_list[4] != "" else pd.DataFrame(data={})
        mut_clone_df    = pd.read_csv(fixed_path_list[5], sep = '\t') if fixed_path_list[5] != "" else pd.DataFrame(data={})
        freq_df         = pd.read_csv(fixed_path_list[6], sep = '\t') if fixed_path_list[6] != "" else pd.DataFrame(data={})
        wgd_df          = pd.read_csv(fixed_path_list[7], sep = '\t') if fixed_path_list[7] != "" else pd.DataFrame(data={})
        feature_df = [pd.read_csv(i, sep = '\t') for i in feature_path_list]
        
        # primary_df      = pd.read_csv(path_list[4], sep = '\t') if path_list[4] != "" else pd.DataFrame(data={})
        # best_res_df     = pd.read_csv(path_list[5], sep = '\t') if path_list[5] != "" else pd.DataFrame(data={})

        cold_set = {'mut':['#b7d5ea','#acc6aa','#266199','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22'],
                    'purity':'Blues',
                    'burden':['#266199','#b7d5ea'],
                    'sig':['#b7d5ea','#acc6aa','#266199','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22'],
                    'copy_num':['#E5E5E5','#b7d5ea','#266199', {'facecolor': 'black', 'alpha': 0.5}, '#7566AA','#000060', '#71a0a5'],
                    'feature':['#DDDDDD','#acc6aa','#71a0a5','#77628c', '#00033F', '#E0CADB','#b7d5ea','888888'],
                    'primary':{'Skin':'#E5E5E5', 'Acral': '#8ec6c5', 'Occult':'#6983aa', 'Mucosal': '#8566aa'},
                    'response':{'SD': '#DDDDDD', 'MR': '#acc6aa', 'PR': '#71a0a5', 'CR':'#77628c', 'PD':  '#00033F'}}

        warm_set = {'mut':['#fbc687','#D4E2C3','#ea907a','#AA1D22','#F4D764','#EBEDAA','#BBBBBB','#aacdbe','#AA7733'],
                    'purity':'Oranges',
                    'burden':['#ea907a','#fbc687'],
                    'sig':['#ea907a','#fbc687','#f4f7c5','#aacdbe'],
                    'copy_num':['#eeeeee','#EEBBBB','#AA1C19', {'facecolor': 'black', 'alpha': 0.2}, '#CCDD99','#575033'],
                    'primary':{'Skin':'#f7e7bd', 'Acral': '#d9c6a5', 'Occult':'#a35d6a', 'Mucosal': '#c26565'},
                    'response':{'SD': '#ea907a', 'MR': '#fbc687', 'PR': '#f4f7c5', 'CR':'#aacdbe', 'PD':  '#ACBF89'}}

        used_set = cold_set if theme == "0" else warm_set

    # info
        #wgd
        wgd_mapping = {'Yes': {'facecolor': 'black', 'alpha': 0.5}, 'No': 'white'}
        #same patient
        indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}
        #copy number
        cna = ['Baseline', 'Allelic amplification', 'Allelic deletion','aCN = 0', 'CN-LOH', 'Complex', 'Multiple']
        cna_mapping = dict(zip(cna, used_set['copy_num']))
        cna_order = self.info['Copy Number Alteration']
        #mut
        category_order = self.info['Mutation Type']
        item = ['Multiple','Silent','Missense','Frameshift indel','In frame indel','Translation start site','Nonstop','Splice site','Nonsense']
        mut_mapping = dict(zip(item, used_set['mut']))
        #signature
        sig_kwargs = {'edgecolor': 'white', 'width': 0.96, 'linewidth': 0.5}
        sig = self.info['Mutation Signature']
        if len(sig) >= 10:
            print(colored('ERROR: The number of mutation signature must be less than 10!\n', 'red'))
            print(colored('CoMut Plotting Halt\n', 'red'))
            return
        sig_mapping = {sig[i]:used_set['sig'][i] for i in range(len(sig))}
        #burden(clonality)
        burden = ['Nonsynonymous','Synonymous']
        bar_mapping = dict(zip(burden, used_set['burden']))
        bar_kwargs = {'width': 0.8}
    
    # features
        feature_mapping = []
        for feature in feature_category:
            fea = self.info[feature]
            if len(fea) > 8:
                print(colored(('ERROR: The number of type of \"'+feature+"\" must be less than 8!\n"), 'red'))
                print(colored('CoMut Plotting Halt\n', 'red'))
                return
            fea_mapping = {fea[i]:used_set['feature'][i] for i in range(len(fea))}
            feature_mapping.append(fea_mapping)
        from matplotlib import rcParams
        custom_rcParams = {'font.family': 'arial','font.size': 14}   
        rcParams.update(custom_rcParams)


    # Plot
        toy_comut = comut.CoMut()
        toy_comut.samples = list(wgd_df['sample'])
        width, fig_info = 0, 0
        FigCol = 2 if width > 10 else 1
        height = int(width/10)*8
        # print(width)
        FigSize = (10, 10)
        
        if not same_patient_df.empty:
        #     width += 1
        #     fig_info += 1
            toy_comut.add_sample_indicators(same_patient_df, name = fixed_category[0], plot_kwargs = indicator_kwargs) 
        if not wgd_df.empty:
            toy_comut.add_categorical_data(wgd_df, name=fixed_category[7], mapping = wgd_mapping)

        if not copy_num_df.empty:
        #     width += 2
        #     fig_info += 7
            toy_comut.add_categorical_data(copy_num_df, name=fixed_category[1], mapping = cna_mapping, category_order = cna_order, tick_style = 'italic') 
        if not mut_type_df.empty:
            # width += 5*int(len(self.info['Mutation Type'])/10)
            # fig_info += 9
            toy_comut.add_categorical_data(mut_type_df, name = fixed_category[2], category_order = category_order, mapping = mut_mapping,tick_style = 'italic')
        
        if not purity_df.empty:
        #     width += 1
            toy_comut.add_continuous_data(purity_df, name = fixed_category[3], mapping = used_set['purity'], value_range = (0, 1))
        
        # for i in range(len(feature_df)):
        #     width += 1
        #     toy_comut.add_categorical_data(feature_df[i], name=feature_category[i], mapping = feature_mapping[i])
        
        if not mut_sig_df.empty:
        #     width += 2
            toy_comut.add_bar_data(mut_sig_df, name = fixed_category[4], mapping = sig_mapping, stacked = True, ylabel = 'Mutational\nsignatures',bar_kwargs = sig_kwargs)
        if not mut_clone_df.empty:
        #     # width += 2
            toy_comut.add_bar_data(mut_clone_df, name = fixed_category[5], mapping = bar_mapping, stacked = True, bar_kwargs = bar_kwargs, ylabel = 'Muts/Mb')
        
        
        toy_comut.plot_comut(figsize = FigSize, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03)
        toy_comut.add_unified_legend(ncol = 1)
        toy_comut.axes[fixed_category[1]].set_xticklabels([])
        toy_comut.figure.savefig(pic+name, dpi = 300, bbox_inches = 'tight')
        # os._exit()
        
        # print(FigSize)
        
        
        # patient_leg = toy_comut.add_axis_legend('Same Patient', bbox_to_anchor = (1.025, 1.2), frameon = False, numpoints = 2)
        # toy_comut.plot_comut(figsize = FigSize, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03)
        # toy_comut.add_unified_legend(ncol = 2)
        # toy_comut.figure.savefig(folder+name, dpi = 300, bbox_inches = 'tight')
        print(colored(("=> Generate CoMut Plot: "+pic+name), 'green'))

 
