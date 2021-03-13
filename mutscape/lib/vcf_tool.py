############################################################################################
# FileName     [ vcf_tool.py ]
# PackageName  [ lib ]
# Synopsis     [ Assisting vcf_filter.py ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

import vcf

def read_vcf(vcf_file):
    ''' Read VCF file by module `vcf`

    Parameters
    ----------
    vcf_file : str

    Returns
    -------
    vcf_reader : vcf.parser.Reader
    '''
    vcf_reader = vcf.Reader(open(vcf_file),'r')
    return vcf_reader
    
# Check the caller of the vcf file
def check_caller(vcf_reader):
    caller = ""
    if "source" in vcf_reader.metadata:
        if vcf_reader.metadata["source"][0] == 'strelka':
            caller = "Strelka2"
        elif vcf_reader.metadata["source"][0] == 'VarScan2':
            caller = "VarScan2"
        elif vcf_reader.metadata["source"][1] == 'Mutect2':
            caller = "Mutect2"
    else:
        if "MuSE_version" in vcf_reader.metadata:
            caller = "MuSE"
        elif "BCOUNT" in vcf_reader.formats:
            caller = "SomaticSniper"
    return caller
