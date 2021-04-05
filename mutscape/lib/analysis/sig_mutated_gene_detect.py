############################################################################################
# FileName     [ sig_mutated_gene_detect.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Control all functions to implement MAF analysis and visualization. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from ..maf_filter import fast_read_maf

class SigMutatedGeneDetection:
    '''MAF analysis: Significantly mutated gene detection 

    Parameters
    ----------
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.
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
            os.system("oncodriveclust -m 3 --cgc src/auxiliary_file/CGC_phenotype.tsv "+\
                      output_folder+"oncodriveCLUST.nonsyn.txt "+output_folder+"oncodriveCLUST.syn.txt "+\
                      "src/auxiliary_file/gene_transcripts.tsv -o "+output_folder+"oncodriveclust_results.tsv\n")
        get_input()
        implement()
        print("\n")
    