############################################################################################
# FileName     [ dataPreprocess.py ]
# PackageName  [ MutScape ]
# Synopsis     [ Control all functions to preprocess input data, including filtering, 
#                combination and transformation. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

from lib.load_tsv import *
from lib.vcf_filter import *

import argparse, textwrap

def main():
    ''' Implement data preprocessing in one single command.

    This function can support both VCFs and MAFs as input data.
    In this function, we provide filtering, combination, transf-
    ormation simultaneously.

    Examples
    --------
    python3 dataPreprocess.py \
    -f examples/tsv/testData.tsv \
    -vf GI "[1,3]" \
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta 

    python3 dataPreprocess.py \
    -f examples/tsv/testData.tsv \
    -vf GI "{1: [*,*], 2 : [1, 300000]}" CI "15,15,0,0,0,0.05,8,8" PA 0 AV 0.9\
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta
    '''
    parser = argparse.ArgumentParser(description='Data preprocessing', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f','--file', help='Input the tsv file.\n\n', required = True, metavar='tsv_file')
    parser.add_argument("-vf", "--vcf_filter",nargs='*', metavar='params',\
                                              help=textwrap.dedent("GI: Genome Interval\n"
                                                                   "CI: Caller Information\n"
                                                                   "PA: Keep or exclude non-PASS tag\n"
                                                                   "AV: Artifact variant filter: FFPE filter (Only for caller = Mutect2)\n\n"))
    parser.add_argument("-c","--combine", action='store_true',help="Input the path for combination VCF files.\nThe path must end with a folder.\n\n")
    parser.add_argument("-v2m","--vcf2maf", nargs=1, help="The first value is the path for transformed MAF files. This path must end with a folder.\n"
                                                             "The second value is \"max_filter_ac\" which is an important parameter(int) when transforming files to MAF.\n\n")
    
    parser.add_argument("-o", "--output", required = True, help="The path for storing output files.\nThis path must end with a folder.\n\n", metavar='out_folder')
    parser.add_argument("-m","--meta", required = True, help="The path for storing metafiles.\nThis path must end with a folder.\n\n", metavar='meta_folder')
    
    args = parser.parse_args()

    flag, category, category_caller = loading_tsv(args.file)
    folder = args.output if args.output[-1:] == '/' else (args.output + '/')
    meta = args.meta if args.meta[-1:] == '/' else (args.meta + '/')

    if flag == 'vcf':
        if not args.combine or not args.vcf2maf:
            raise ValeError('[MutScape] Command -c, -v2m must required if inputs are VCFs.')
        filter_list = []
        vcf_filter(args.vcf_filter, category, category_caller, meta)
        
if __name__ == '__main__':
    main()
