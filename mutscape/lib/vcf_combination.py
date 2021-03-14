############################################################################################
# FileName     [ vcf_combination.py ]
# PackageName  [ lib ]
# Synopsis     [ Combine 2 or more vcf files for a sample with different variants ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

from .vcf_tool import *
from termcolor import colored
from datetime import datetime
import os

def generate_header(read_list, caller_list, NT):
    ''' Generate combined header in the VCF

    Parameters
    ----------
    read_list : list
        The type of all items in `read_list` is `vcf.parser.Reader`.
    caller_list : list
    NT : list
        The length of `NT` is 2 which includes the name for NORMAL and TUMOR.
    
    Returns
    -------
    read_header : vcf.parser.Reader
    '''
    def combine_metadata(read_list, caller_list, NT):
        ''' Combine metadata in VCF header

        Parameters
        ----------
        read_list : list
            The type of all items in `read_list` is `vcf.parser.Reader`.
        caller_list : list
        NT : list
            The length of `NT` is 2 which includes the name for NORMAL and TUMOR.

        Returns
        -------
        Combine_metadata : collections.OrderedDict
        '''
        if len(read_list) != len(caller_list):
            print(colored("Error when reading combination files!!", "red"))
            return
        Combine_metadata = read_list[0].metadata.copy()
        Combine_metadata.clear()
        # filedate
        Combine_metadata['filedate'] = []
        Combine_metadata['filedate'].append(datetime.now().strftime("%Y%m%d"))
        # Fileformat
        fileformat =[i.metadata['fileformat'] for i in read_list if 'fileformat' in i.metadata]
        if not any([i != fileformat[0] for i in fileformat]):
            Combine_metadata['fileformat'] = [fileformat[0]]
        # Source
        Combine_metadata['source'] = [caller_list]
        # Source version
        Combine_metadata['source_version'] = [[]]
        for i in range(len(caller_list)):
            if caller_list[i] == "MuSE":
                Combine_metadata['source_version'][0].append(read_list[i].metadata['MuSE_version'][0])
            elif caller_list[i] == "Mutect2":
                Combine_metadata['source_version'][0].append(read_list[i].metadata['MutectVersion'][0])
            elif caller_list[i] == "Strelka2":
                Combine_metadata['source_version'][0].append(read_list[i].metadata['source_version'][0])
            else:
                Combine_metadata['source_version'][0].append('.')
        # Reference
        Combine_metadata['reference'] = []
        for idx, file in enumerate(read_list):
            if 'reference' in file.metadata:
                Combine_metadata['reference'].append(file.metadata['reference'])
            else:
                Combine_metadata['reference'].append('.')
        # NORMAL and TUMOR
        Combine_metadata['NORMAL'] = [NT[0]]
        Combine_metadata['TUMOR'] = [NT[1]]
        print(type(Combine_metadata))
        return Combine_metadata
    def combine_infos_filters_formats(read_list):
        ''' Combine info, filter and format in VCF header

        Parameters
        ----------
        read_list : list
            The type of all items in `read_list` is `vcf.parser.Reader`.
        
        Returns
        -------
        combine_info : collections.OrderedDict
        combine_filter : collections.OrderedDict
        combine_format : collections.OrderedDict
        '''
        combine_info = read_list[[len(i.infos) for i in read_list].index(max([len(i.infos) for i in read_list]))].infos.copy()
        combine_filter = read_list[[len(i.filters) for i in read_list].index(max([len(i.filters) for i in read_list]))].filters.copy()
        combine_format = read_list[[len(i.formats) for i in read_list].index(max([len(i.formats) for i in read_list]))].formats.copy()
        for idx, file in enumerate(read_list):
            for info in file.infos:
                if info not in combine_info:
                    combine_info[info] = file.infos[info] 
            for flt in file.filters:
                if flt not in combine_filter:
                    combine_filter[flt] = file.filters[flt]  
            for fmt in file.formats:
                if fmt not in combine_format:
                    combine_format[fmt] = file.formats[fmt]
        combine_info['CALLS'] = vcf.parser._Info(id= "CALLS", num = 1, type= "String", desc="Sign caller names that contain this variant", source = None, version = None)
        combine_info['REJECT'] = vcf.parser._Info(id= "REJECT", num = 1, type= "String", desc="Sign caller names that reject this variant", source = None, version = None)
        print(type(combine_info), type(combine_filter),type(combine_format))
        return combine_info, combine_filter, combine_format
    read_header = read_list[0]
    read_header.metadata = combine_metadata(read_list, caller_list, NT)
    read_header.infos, read_header.filters, read_header.formats = combine_infos_filters_formats(read_list)
    return read_header

# Get Somatic list
# A few dictionaries in this list, one dictionary for a file
# If don't want to remove chrX, chrY, set Somatic to False
def get_somatic_list(read_list, caller_list, Somatic = True):
    Somatic_list = []
    chromo_num_list =[str(i) for i in range(1,23)]#+["M","X","Y"]
    for read in read_list:
        new_dict = dict(zip(chromo_num_list, [[] for i in range(len(chromo_num_list))]))
        for record in read:
            if str(record.CHROM) in chromo_num_list:
                new_dict[str(record.CHROM)].append(record)
        Somatic_list.append(new_dict)
    somatic_print = [([i]+[len(file[i]) for file in Somatic_list])for i in chromo_num_list]
    num_list = [sum(somatic_print[ic][1:]) for ic, i in enumerate(somatic_print)]
    somatic_print = [chrom for chrom in somatic_print if any(item != 0 for item in chrom[1:])]
    print(colored("The number of each CHROM in each input vcf file:\n", 'green'))
    print(tabulate(somatic_print, headers=['#']+caller_list, tablefmt='orgtbl'))
    print("\n")
    return Somatic_list, num_list


def vcf_combination(sample_content_list, caller_list, output_file):
    Read_list = [read_vcf(x) for x in sample_content_list[2]]
    Read_header = generate_header(Read_list, caller_list, [sample_content_list[0], sample_content_list[1]])
    vcf_writer = vcf.Writer(open(output_file,"w"), Read_header)
    INFO_list = [infos for infos in Read_header.infos]
    FORMAT_list = ['GT', 'AD', 'DP', 'AF']
    FORMAT_list.extend([x for x in Read_header.formats if x not in FORMAT_list])
    Sample_list = [sample_content_list[0], sample_content_list[1]]
    NEW_FORMAT_STR = ":".join(FORMAT_list)
    Somatic_list, num_list = get_somatic_list(Read_list, caller_list)
    Merged_list = merge_and_sort(Somatic_list, NEW_FORMAT_STR, INFO_list, caller_list, FORMAT_list, Sample_list, num_list)
    print(colored("\nWriting to ...."+output_file+"\n",'yellow'))
    for idx, chrom in enumerate(Merged_list):
        for record in chrom:
            vcf_writer.write_record(record)
    vcf_writer.close()
    print(colored(("=> Finish combining file: "+output_file+'\n'), 'green'))

def all_combine(category, category_caller, meta):
    print(colored("\nStart VCF combination....\n", "yellow"))
    combine_output_list = []
    for idx in range(len(category)):
        path_name = meta + category[idx][0] + '_' + category[idx][1] + '_combination.vcf'
        combine_output_list.append(path_name)

    for idx, sample in enumerate(category):
        vcf_combination(category[idx], category_caller[idx], combine_output_list[idx])
    ## Deal with "At Least...." condition
    combine_filter_filelist = [file[:-4]+"_f.vcf" for file in combine_output_list]
    for idx, file in enumerate(combine_output_list):
        del_cou = 0
        vcf_read = read_vcf(file)
        vcf_writer = vcf.Writer(open(combine_filter_filelist[idx],"w"), vcf_read)
        for record in vcf_read:
            calls = len(record.INFO['CALLS'].split('_')) if record.INFO['CALLS']!= None else 0
            reject = len(record.INFO['REJECT'].split('_')) if record.INFO['REJECT'] != None else 0
            if calls >= category[idx][3] and reject <= category[idx][4]:
                vcf_writer.write_record(record)
            else:
                del_cou+=1
        vcf_writer.close()
        print("NOTICE: "+str(del_cou)+" data have been removed from "+combine_output_list[idx])
    print("\n")
    return combine_filter_filelist
