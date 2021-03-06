############################################################################################
# FileName     [ vcf2maf.py ]
# PackageName  [ lib ]
# Synopsis     [ Use vcf2maf module to convert VCF file to MAF ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021.3 ]
############################################################################################

import vcfpy
from vcfpy import header
import os
from termcolor import colored
import multiprocessing
import subprocess
from .vcf_tool import *

def vcf2vep2maf(vcf_file_list, maf_file_list, path, category, max_filter_ac):
    ''' Write a .sh file to automatically implement vcf2maf utility.

    Parameters
    ----------
    vcf_file_list : list
    maf_file_list : list
    path : str
    category : list
    max_filter_ac : str
    '''
    run_file = open(path+"run.sh", "w")
    run_file.write("YELLOW='\\033[0;33m'\n")
    run_file.write("GREEN='\\033[0;32m'\n")
    run_file.write("NC='\\033[0m'\n")

    perl_path = input("Please enter vcf2maf.pl's path (If the path is '../../vcf2maf.pl', just press ENTER.): ")
    vep_path = input("Please enter the path of vep (Folder containing the vep script): ")
    apath = (subprocess.check_output("which vep", shell=True)).decode("utf-8") 

    if perl_path == "":
        perl_path = "../../vcf2maf.pl"
    if vep_path == '':
        vep_path = apath[:-4]
        print(vep_path)
    fork = str(multiprocessing.cpu_count())
    if isinstance(vcf_file_list, list) and isinstance(maf_file_list, list) and len(vcf_file_list) == len(maf_file_list):
        for index, vcf_file in enumerate(vcf_file_list):
            vcf_read_samples = (read_vcf(vcf_file)).samples
            run_file.write("printf \"${YELLOW}\nStart transforming file:\nVCF: "+vcf_file+"\nMAF: "+maf_file_list[index]+"\n\n${NC}\"\n\n")
            run_file.write("printf \"VCF file's size: "+str(os.stat(vcf_file).st_size/10**6)+" MB\"\n\n")
            run_file.write("printf \"\n\"\n")
            run_file.write("perl "+perl_path+" \\\n"+ \
                            "--tumor-id "+vcf_read_samples[1]+" \\\n"+\
                            "--normal-id "+vcf_read_samples[0]+" \\\n"+\
                            "--vcf-tumor-id "+vcf_read_samples[1]+" \\\n"+\
                            "--vcf-normal-id "+vcf_read_samples[0]+" \\\n"+\
                            "--vep-path "+vep_path+" \\\n"+\
                            "--max-filter-ac "+max_filter_ac+" \\\n"+\
                            "--input-vcf "+vcf_file+" \\\n"+\
                            "--output-maf "+maf_file_list[index]+" \\\n")
            if fork != "":
                run_file.write("--vep-forks "+fork+" \\\n")
            run_file.write("\n")
            run_file.write("printf \"${GREEN}\n=> Finish transforming file: "+vcf_file+" to "+maf_file_list[index]+"${NC}\n\n\"\n\n")
    else:
        print(len(vcf_file_list),len(maf_file_list))
        print(colored("ERROR: Different number of VCF files and MAF files!\n", "red"))
    run_file.close()
    os.system("sh "+path+"/run.sh\n")
    os.system("rm "+path+"/run.sh")

def vcf2maf(if_vcf2maf, combine_filter_filelist, folder, category):
    ''' Transform VCF file to MAF file using vcf2maf utility.

    Parameters
    ----------
    combine_filter_filelist : list
        The list of VCF files.
    folder : str
    category : list
    num : str
        max_filter_ac

    Returns
    -------
    maf_output_list
    '''
    if if_vcf2maf:
        print(colored("Start transforming VCF to MAF....\n", "yellow"))
        print("WARNING: This transformation tool must be implemented in the direction of \"vcf2maf-1.6.20\"!\n")
        maf_output_list = []
        for idx, file in enumerate(combine_filter_filelist):
            fileName = file[:-4]+"_2maf.maf"
            maf_output_list.append(fileName)
        vcf2vep2maf(combine_filter_filelist, maf_output_list, folder, category, if_vcf2maf[0])
        return maf_output_list
    else:
        return []