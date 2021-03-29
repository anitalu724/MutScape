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
from tabulate import tabulate
from tqdm import tqdm
import os
import sys, operator
import numpy as np

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

def get_somatic_list(read_list, caller_list, Somatic = True):
    '''Use `read_list` to get somatic list

    Parameters
    ----------
    read_list : list
        The type of all items in `read_list` is `vcf.parser.Reader`.
    caller_list : list
    Somatic : bool(Optional)
        Set `False` if do not want to remove chrX, chrY, chrM.

    Returns
    -------
    Somatic_list : list
        The size of `Somatic_list` is equal to the size of `read_list`.
        In `Somatic_list`, each item is a dictionary which has classified 
        records by `record.CHROM`.
    num_list : list
        The size of `num_list` is 22 which calculate the total number of each
        # of CHROM.
    '''
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

def merge_and_sort(somatic_list, NEW_FORMAT_STR, INFO_list, caller_list, FORMAT_list, Sample_list, num_list):
    ''' Merge similar records in different VCFs, and sort them in the correct order.

    Parameters
    ----------
    somatic_list : list
        The size of `Somatic_list` is equal to the size of `read_list`.
        In `Somatic_list`, each item is a dictionary which has classified 
        records by `record.CHROM`.
    NEW_FORMAT_STR : str
        The combination of `FORMAT` in each VCFs. Every item is linked by `:`.
        Ex: GT:AD:AP:AF:F1R2:F2R1
    INFO_list : list
        A list of `INFO` which has combined those in each VCFs. We have added two 
        items `CALLS` and `REJECT` to store the origin of each record.
        Ex: ...;CALLS='MuSE_Mutect2';REJECT='Strelka2'
    caller_list : list
    FORMAT_list : list
        A list of `FORMAT` which has combined those in each VCFs.
    Sample_list : list
        Generally, the size of `Sample_list` is 2. 
        Ex: ['NORMAL', 'TUMOR']
    num_list : list
        The size of `num_list` is 22 which calculate the total number of each
        # of CHROM.

    Returns
    -------
    merged_list : list
        22 lists are in merged_list. Each of the list is related to 1~22 
        chromosome respectively.
    '''
    # Part2 in vcf_combination
# change info for each record
    def get_new_info(record_info, info_list, caller, PASS):
        '''Add `CALLS` and `REJECT` into `INFO` column in VCF.

        Parameters
        ----------
        record_info : dict
        info_list : list
        caller : str
            MuSE/ Mutect2/ SomaticSniper/ Strelka2/ VarScans
        PASS : bool

        Returns
        -------
        record_info : dict
            New `record_info` which has added
        '''
        for info in info_list:
            if info not in record_info:
                record_info[info] = None
        record_info['CALLS'] = caller
        record_info['REJECT'] = None if PASS else caller
        return record_info
    def get_new_samples(record, format_list, sample_list):
        origin_format = record.FORMAT.split(':')
        for sample in record.samples:
            new_list= []
            for item in format_list:
                if item not in origin_format:
                    new_list.append(None)
                else:
                    new_list.append(sample[item])
            sample.data = tuple(new_list)
        return record.samples
    def merge_info(candidate, INFO_list, caller_list, iMin, similar):
        # merge similar record to one
        merge_INFO = candidate[iMin].INFO.copy()
        for s in similar:
            for info in candidate[s].INFO:
                if info not in merge_INFO:
                    merge_INFO[info] =  candidate[s].INFO[info]   
                else:
                    if merge_INFO[info] != candidate[s].INFO[info]:
                        merge_INFO[info] = None
        for info in INFO_list:
            if info not in merge_INFO:
                merge_INFO[info] = None
        c = [iMin]+similar
        merge_INFO['CALLS'] = "_".join(list(np.array(caller_list)[c]))
        reject = "_".join(list(np.array(caller_list)[[i for i in c if len(candidate[i].FILTER) != 0]]))
        merge_INFO['REJECT'] = reject if reject != "" else None
        return merge_INFO
    def merge_samples(candidate, FORMAT_list, Sample_list, iMin, similar):
        origin_format = candidate[iMin].FORMAT.split(':')
        similar_format = [candidate[s].FORMAT.split(':') for s in similar]
        merge_Sample = candidate[iMin].samples.copy()
        for idx, sample in enumerate(merge_Sample):
            new_list= []
            for item in FORMAT_list:
                sim = []
                for ids, s in enumerate(similar):
                    if item in similar_format[ids]:
                        sim.append(candidate[s].samples[idx][item])
                if len(sim) == 0:
                    if item not in origin_format:
                        new_list.append(None)
                    else:
                        new_list.append(sample[item])
                else:
                    if item in origin_format:
                        sim.append(sample[item])
                        if any([s != sim[0] for s in sim]):
                            new_list.append(None)
                        else:
                            new_list.append(sample[item])
                    else:
                        if any([s != sim[0] for s in sim]):
                            new_list.append(None)
                        else:
                            new_list.append(sim[0])
            sample.data = tuple(new_list)
        return merge_Sample
    list_order = [str(i) for i in range(1,23)]
    merged_list, no = [[] for i in range(len(list_order))], len(somatic_list)
    print(colored("Merging....\n",'yellow'))
    for idx, CHROM in enumerate(list_order):
        if num_list[idx] != 0:
            pbar = tqdm(total = num_list[idx], desc = "chr"+CHROM+" ")
            candidate, candidate_POS, empty = [], [], False
            for file in somatic_list:
                if len(file[CHROM]) != 0:
                    candidate.append(file[CHROM][0])
                    candidate_POS.append(file[CHROM][0].POS)
                else:
                    candidate.append(sys.maxsize)
                    candidate_POS.append(sys.maxsize)
            while not empty:
                iMin, Min = min(enumerate([i for i in candidate_POS]), key=operator.itemgetter(1))
                if Min != sys.maxsize:
                    similar = list(set(np.where([(i.POS == Min and i.REF == candidate[iMin].REF and i.ALT == candidate[iMin].ALT) for i in candidate if i != sys.maxsize])[0])-set([iMin])-set([ic for ic, i in enumerate(candidate) if i == sys.maxsize]))
                    if len(similar) != 0:
                        #Deal with similar record
                        somatic_list[iMin][CHROM].pop(0)
                        for item in similar:
                            somatic_list[item][CHROM].pop(0)
                        candidate[iMin].INFO = merge_info(candidate, INFO_list, caller_list, iMin, similar)
                        candidate[iMin].samples = merge_samples(candidate, FORMAT_list, Sample_list, iMin, similar)
                        candidate[iMin].FORMAT = NEW_FORMAT_STR
                        merged_list[idx].append(candidate[iMin])
                        candidate[iMin] = somatic_list[iMin][CHROM][0] if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        candidate_POS[iMin] = somatic_list[iMin][CHROM][0].POS if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        for item in similar:
                            candidate[item] = somatic_list[item][CHROM][0] if len(somatic_list[item][CHROM]) > 0 else sys.maxsize
                            candidate_POS[item] = somatic_list[item][CHROM][0].POS if len(somatic_list[item][CHROM]) > 0 else sys.maxsize
                        pbar.update(1+len(similar))
                    else:   # move this chosen record into merged_list
                        somatic_list[iMin][CHROM].pop(0)
                        candidate[iMin].INFO = get_new_info(candidate[iMin].INFO, INFO_list, caller_list[iMin], len(candidate[iMin].FILTER) == 0 )
                        candidate[iMin].samples = get_new_samples(candidate[iMin], FORMAT_list, Sample_list)
                        candidate[iMin].FORMAT = NEW_FORMAT_STR
                        merged_list[idx].append(candidate[iMin])
                        candidate[iMin] = somatic_list[iMin][CHROM][0] if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        candidate_POS[iMin] = somatic_list[iMin][CHROM][0].POS if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize   
                        pbar.update(1)   
                else:
                    empty = True
    return merged_list


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
