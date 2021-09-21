############################################################################################
# FileName     [ known_cancer_gene_anno.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Annotate known cancer gene. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 9 ]
############################################################################################

from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd

class KnownCancerGeneAnnotation:
    ''' Known cancer gene annotation
    Arguments:
        maf_file            {string}        -- The input MAF file for all data.
        output_folder       {string}        -- The path for output files.

    Parameters:
 
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.  

    Output files
    ------------
    output :
        kcga.output.maf
            
    ''' 
    def __init__(self, maf_file):
        print(colored(('\nStart Known_Cancer_Gene_Annotation....'), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
    def annotation(self, output_folder):
        anno_file = 'lib/auxiliary/cancerGeneList.txt'
        origin_file = self.df
        output_file = output_folder+'kcga.output.maf'
        anno = pd.read_csv(anno_file, header = 0, sep = '\t')
        origin_file['#_of_occurrence_within_resources'] = 'No'
        origin_file['OncoKB_Annotated'] = 'No'
        origin_file['Is_Oncogene'] = 'No'
        origin_file['Is_Tumor_Suppressor_Gene'] = 'No'
        origin_file['MSK-IMPACT'] = 'No'
        origin_file['MSK-HEME'] = 'No'
        origin_file['FOUNDATION_ONE'] = 'No'
        origin_file['FOUNDATION_ONE_HEME'] = 'No'
        origin_file['Vogelstein'] = 'No'
        origin_file['SANGER_CGC(05/30/2017)'] = 'No'
        List = ['#_of_occurrence_within_resources', 'OncoKB_Annotated', 'Is_Oncogene','Is_Tumor_Suppressor_Gene','MSK-IMPACT','MSK-HEME','FOUNDATION_ONE','FOUNDATION_ONE_HEME','Vogelstein','SANGER_CGC(05/30/2017)']
        for item in range(len(anno)):
            origin_file.loc[origin_file['Hugo_Symbol'] == anno['Hugo Symbol'][item],List] = anno['# of occurrence within resources (Column D-J)'][item], anno['OncoKB Annotated'][item],anno['Is Oncogene'][item], anno['Is Tumor Suppressor Gene'][item], anno['MSK-IMPACT'][item], anno['MSK-HEME'][item], anno['FOUNDATION ONE'][item], anno['FOUNDATION ONE HEME'][item],anno['Vogelstein'][item],anno['SANGER CGC(05/30/2017)'][item]
        origin_file.to_csv(output_file, sep = '\t', index = False)
        print(colored('=> Generate output file: ', 'green'))
        print(colored(('   '+output_file+'\n'), 'green'))
