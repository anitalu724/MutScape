############################################################################################
# FileName     [ dataPreprocess.py ]
# PackageName  [ MutScape ]
# Synopsis     [ Control all functions to preprocess input data, including filtering, 
#                combination and transformation. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

from lib.load_tsv import *
# from lib.vcf_filter import *

import argparse

def main():
    parser = argparse.ArgumentParser(description='Data preprocessing', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f','--file', help='Input the tsv file.\n\n', required = True, metavar='TSV file')
    
    args = parser.parse_args()

    flag, category, category_caller = loading_tsv(args.file)
    

if __name__ == '__main__':
    main()