############################################################################################
# FileName     [ sig_mutated_gene_detect.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Detect significantly mutated gene. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 9 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import os

class SigMutatedGeneDetection:
    '''Significantly mutated gene detection 

    Arguments:
        maf_file            {string}        -- The input MAF file for all data.
        output_folder       {string}        -- The path for output files.

    Parameters:
        self.head           {string}        -- The column names of MAF file.
        self.df             {pd.DataFrame}  -- The data for the MAF file.
    
    Output:
        oncodriveCLUST.nonsyn.txt
        oncodriveCLUST.syn.txt
        oncodriveclust_results.tsv
    '''
    def __init__(self, maf_file):
        print(colored(("\nStart Significantly_Mutated_Gene_Detection...."), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
        
    def oncodriveCLUST(self, output_folder):
        def get_input():
            selected_col = self.df[['Hugo_Symbol', 'Variant_Classification','Tumor_Sample_Barcode', 'Transcript_ID', 'Gene', 'Protein_position']]
            selected_col.columns = ['symbol','Sample','Tumor_Sample_Barcode', 'transcript', 'gene', 'position']
            list1, list2 = [],[]
            for idx in range(selected_col.shape[0]):
                if type(selected_col.iloc[idx]['position']) != float:
                    pos = selected_col.iloc[idx]['position'].split("/")[0]
                    if "-" not in pos:
                        selected_col.iloc[idx]['position'] = pos
                        if selected_col.iloc[idx]['Sample'] == "Silent":
                            list2.append(idx)
                        else:
                            list1.append(idx)
            file_1, file_2 = selected_col.iloc[list1, :], selected_col.iloc[list2, :]
            for idx in range(file_1.shape[0]):
                if file_1.iloc[idx]['Sample'] == 'Nonsense_Mutation':
                    file_1.iloc[idx]['Sample'] = "stop"
                else:
                    file_1.iloc[idx]['Sample'] = "non-synonymous"
            for idx in range(file_2.shape[0]):
                file_2.iloc[idx]['Sample'] = "synonymous"
            file_1 = file_1.loc[:,['symbol', 'gene', 'transcript', 'Tumor_Sample_Barcode', 'Sample', 'position']]
            file_1.columns = ['symbol', 'gene', 'transcript', 'Sample', 'ct', 'position']
            file_2 = file_2.loc[:,['symbol', 'gene', 'transcript', 'Tumor_Sample_Barcode', 'Sample', 'position']]
            file_2.columns = ['symbol', 'gene', 'transcript', 'Sample', 'ct', 'position']
            file_1.to_csv(output_folder+"oncodriveCLUST.nonsyn.txt", sep = '\t', index = False)
            file_2.to_csv(output_folder+"oncodriveCLUST.syn.txt", sep = '\t', index = False)
            print(colored(("=> Generate OncodriveCLUST's input files:"), 'green'))
            print(colored(("   "+output_folder+"oncodriveCLUST.nonsyn.txt"), 'green'))
            print(colored(("   "+output_folder+"oncodriveCLUST.syn.txt\n"), 'green'))
        def implement():
            os.system("oncodriveclust -m 3 --cgc lib/auxiliary/CGC_phenotype.tsv "+\
                      output_folder+"oncodriveCLUST.nonsyn.txt "+output_folder+"oncodriveCLUST.syn.txt "+\
                      "lib/auxiliary/gene_transcripts.tsv -o "+output_folder+"oncodriveclust_results.tsv\n")
        get_input()
        implement()
        print("\n")
    