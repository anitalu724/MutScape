############################################################################################
# FileName     [ load_tsv.py ]
# PackageName  [ lib ]
# Synopsis     [ Loading data from TSV file ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

from termcolor import colored
from tabulate import tabulate
import pandas as pd

def read_tsv(tsv_file):
    """ Read TSV file to determine whether to implement VCF or MAF.

    Parameters
    ----------
    tsv_file : str

    Returns
    -------
    flag : str (vcf/maf)
    category : list
    category_caller : list

    Raises
    ------
    ValueError if the format of TSV file is wrong.
    """
    file = pd.read_csv(tsv_file, sep='\t')
    category, category_print, category_caller = [], [], []
    flag = ""
    if file.shape[1] == 1:      # MAF only
        for i in range(file.shape[0]):
            if file.loc[i, 'MAF'][-3:]!= "maf":
                raise ValueError('[MutScape] The format of TSV file is wrong!')
            category.append(file.loc[i, 'MAF'])
        length = [[str(len(category))]]
        print(colored("\nThe input tsv file\'s content:\n", "green"))
        print(tabulate(length, headers=['# of MAF files'], tablefmt='orgtbl'))
        print("\n")
        flag = "maf"
        
    elif file.shape[1] == 9:    # VCF and MAF
        for i in range(file.shape[0]):
            sample, vcf_file, caller_list = [], [], []
            sample.append(file.loc[i,'NORMAL'])
            sample.append(file.loc[i,'TUMOR'])
            variant_list = ['MuSe', 'Mutect2', 'SomaticSniper', 'Strelka2', 'VarScan2']
            for j in variant_list:
                if not isinstance(file.loc[i, j], float):
                    if file.loc[i, j][-3:] != "vcf":
                        raise ValueError('[MutScape] The format of TSV file is wrong!')
                    vcf_file.append(file.loc[i, j])
                    caller_list.append(j)
            sample.append(vcf_file)
            sample.append(file.loc[i,'At Least # CALLS'])
            sample.append(file.loc[i,'At Most # REJECT'])
            category.append(sample)
            category_caller.append(caller_list)
        for i in category:
            if i[3] == 0 or i[3] > len(i[2]):
                raise ValueError('[MutScape] The format of TSV file is wrong!')
            sample = []
            sample.extend((i[0], i[1], len(i[2]), i[3], i[4]))
            category_print.append(sample)
        print(colored("\nThe input tsv file\'s content:\n", "green"))
        print(tabulate(category_print, headers=['NORMAL', 'TUMOR', '# of VCF files', 'At least # of variants', 'At most # of REJECT'], tablefmt='orgtbl'))
        flag = "vcf"
    else:
        raise ValueError('[MutScape] The format of TSV file is wrong!')
    return flag, category, category_caller


def loading_tsv(tsv_file):
    print(colored("\nReading TSV file....", "yellow"))
    return read_tsv(tsv_file)