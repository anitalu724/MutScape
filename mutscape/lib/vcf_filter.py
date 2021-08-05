############################################################################################
# FileName     [ vcf_filter.py ]
# PackageName  [ lib ]
# Synopsis     [ Sifting records from VCFs ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

import vcf
import numpy as np
from .vcf_tool import *
from termcolor import colored

def get_filter_data(vcf_flt, num=4):
    ''' Interpret vcf filtering parameters into a list.

    Parameters
    ----------
    vcf_flt : list
        The maximum size of `vcf_flt` is 8 and ths size of `vcf_flt`
        must be even.
    num : int, optional
        The number of filters provided in this function.
    
    Returns
    -------
    flt_list : list
        The size of `flt_list` is 4. If no parameter is entered, the 
        item will be recorder to `False`.

    Raises
    ------
    ValueError if the format for command -vf is wrong.
    '''
    flt_list = [False]*num
    genome_length = [249250621, 243199373, 198022430, 191154276, 180915260,
                     171115067, 159138663, 146364022, 141213431, 135534747,
                     135006516, 133851895, 115169878, 107349540, 102531392,
                     90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
    if len(vcf_flt)%2 != 0:
        raise ValueError('[MutScape] Wrong format for command -vf!')
    for idx, data in enumerate(vcf_flt):
        if idx%2 == 0:
            if data == "GI":
                info = vcf_flt[idx+1]
                info = info.replace(" ","")
                if "{" not in info and "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    new_info = []
                    for i in info:
                        if ":" not in i:
                            new_info.append(int(i))
                        else:
                            s, e = i.split(":")
                            for no in range(int(s),int(e)+1):
                                new_info.append(no)
                    flt_list[0] = new_info
                elif "{" in info:
                    info = info.strip('{').strip('}').split("],")
                    new_dict = {}
                    for i in info:
                        i = i.split(":")
                        s, e = i[1].split(",")
                        s = 1 if (s.strip('[') == '*') else int(s.strip('['))
                        e = genome_length[int(i[0])-1] if (e.strip(']') == '*') else int(e.strip(']'))
                        new_dict[str(i[0])] = [s, e]
                    flt_list[0] = new_dict
            elif data == "CI":
                flt_list[1] = [x for x in (vcf_flt[idx+1].replace(" ", "").replace("[", "").replace("]", "")).split(",")]
                for i in range(len(flt_list[1])):
                    flt_list[1][i] = float(flt_list[1][i]) if flt_list[1][i] != '*' else flt_list[1][i]
            elif data == "PA":
                flt_list[2] = bool(int(vcf_flt[idx+1]))
            elif data == "AV":
                flt_list[3] = float(vcf_flt[idx+1])
    return flt_list

def nt_test(vcf_read, col_number, caller):
    col_0, col_1 = [], []
    if caller == "SomaticSniper" or caller == "MuSE":
        for record in vcf_read:
            col_0.append(record.samples[0]['SS'])
            col_1.append(record.samples[1]['SS'])
        if all((i == 2) for i in col_0) and not all ((i == 2) for i in col_1):
            return True
        else:
            return False
    elif caller == "Strelka2":
        return False
    else:
        for record in vcf_read:
            col_0.append(record.samples[0]['GT'])
            col_1.append(record.samples[1]['GT'])
        if all((i == "0/0" or i == "0|0") for i in col_1) and not all ((i == "0/0" or i == "0|0") for i in col_0):
            return True
        else:
            return False

def genome_interval(record, select):
    ''' Genome Interval(GI) filter

    Parameters
    ----------
    record : vcf.model._Record
        A example of `record` is `Record(CHROM=1, POS=1560973, REF=C, ALT=[T])`.
    select : list/dict
        Two type of `select` are supported. The type `list` means 
        `chromosome selection`, while type `dict` means `interval-
        selection`.
    
    Returns
    -------
    PASS : bool
    '''
    
    PASS = True
    if type(select) == list:
        selected_list = []
        for i in select:
            selected_list.append(str(i))
        PASS = True if record.CHROM in selected_list else False
    elif type(select) == dict:
        if record.CHROM in select.keys():
            if select[record.CHROM][0] <= record.POS and record.POS <= select[record.CHROM][1]:
                PASS = True
            else:
                PASS = False
        else:
            PASS = False
    return PASS
    

def caller_info(call, record, para_list, nt_name):
    ''' Caller Information(CI) filter

    Parameters
    ----------
    call : str (MuSE / Mutect2 / SomaticSniper / Strelka2 / VarScan2)
    record : vcf.model._Record
        A example of `record` is `Record(CHROM=1, POS=1560973, REF=C, ALT=[T])`.
    para_list : list
        DP_N : float/*
            Coverage of the normal.
        DP_T : float/*
            Coverage of the tumor.
        AD_N : float/*
            Number of reads supporting the alternative allele in the normal.
        AD_T : float/*
            Number of reads supporting the alternative allele in the tumor.
        AF_N : float/*
            Fraction of variant supporting reads in the normal.
        AF_T : float/*
            Fraction of variant supporting reads in the tumor.
        NLOD : float/*
            Normal LOD score.
        TLOD : float/*
            Tumor LOD score.
    nt_name : list
    

    Returns
    -------
    PASS : bool

    Raises
    ------
    ValueError if `call` is not in the list below.
        [MuSE / Mutect2 / SomaticSniper / Strelka2 / VarScan2]
    '''
    if call == 'Mutect2' or call == 'Dragen':
        for idx, para in enumerate(para_list):
            if idx == 0 and para != '*':    # DP_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['DP'] <= para:
                        return False
            elif idx == 1 and para != '*':    # DP_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['DP'] <= para:
                        return False
            elif idx == 2 and para != '*':    # AD_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['AD'][1] > para:
                        return False
            elif idx == 3 and para != '*':    # AD_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['AD'][1] < para:
                        return False
            elif idx == 4 and para != '*':    # AF_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and type(sample['AF']) is float:
                        if sample['AF'] > para:
                            return False
            elif idx == 5 and para != '*':    # AF_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and type(sample['AF']) is float:
                        if sample['AF'] < para:
                            return False
            elif idx == 6 and para != '*':    # NLOD
                if record.INFO['NLOD'][0] < para:
                    return False
            elif idx == 7 and para != '*':    # TLOD
                if record.INFO['TLOD'][0] < para:
                    return False
    elif call == 'MuSE':
        for idx, para in enumerate(para_list):
            if idx == 0 and para != '*':    # DP_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['DP'] <= para:
                        return False
            elif idx == 1 and para != '*':    # DP_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['DP'] <= para:
                        return False
            elif idx == 2 and para != '*':    # AD_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['AD'][1] > para:
                        return False
            elif idx == 3 and para != '*':    # AD_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['AD'][1] < para:
                        return False
            elif idx == 4 and para != '*':    # AF_N
                for sample in record.samples:
                    if sample.sample == nt_name[0]:
                        AF_N = round(sample['AD'][1]/sample['DP'],3)
                        if AF_N > para:
                            return False
            elif idx == 5 and para != '*':    # AF_T
                for sample in record.samples:
                    if sample.sample == nt_name[1]:
                        AF_T = round(sample['AD'][1]/sample['DP'],3)
                        if AF_T < para:
                            return False
    elif call == 'SomaticSniper':
        for idx, para in enumerate(para_list):
            if idx == 0 and para != '*':    # DP_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['DP'] <= para:
                        return False
            elif idx == 1 and para != '*':    # DP_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['DP'] <= para:
                        return False
            elif idx == 2 and para != '*':    # AD_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['AD'][1] > para:
                        return False
            elif idx == 3 and para != '*':    # AD_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['AD'][1] < para:
                        return False
            elif idx == 4 and para != '*':    # AF_N
                for sample in record.samples:
                    if sample.sample == nt_name[0]:
                        AF_N = round(sample['AD'][1]/sample['DP'],3)
                        if AF_N > para:
                            return False
            elif idx == 5 and para != '*':    # AF_T
                for sample in record.samples:
                    if sample.sample == nt_name[1]:
                        AF_T = round(sample['AD'][1]/sample['DP'],3)
                        if AF_T < para:
                            return False
    elif call == 'Strelka2':
        AD_N, AD_T = 0,0
        for sample in record.samples:
            _AD = []
            if sample["AU"][0] > 0:
                _AD.append(sample["AU"][0]) 
            if sample["CU"][0] > 0:
                _AD.append(sample["CU"][0]) 
            if sample["GU"][0] > 0:
                _AD.append(sample["GU"][0]) 
            if sample["TU"][0] > 0:
                _AD.append(sample["TU"][0])
            if len(_AD) == 2:
                _AD = [_AD[1],_AD[0]] if _AD[1] > _AD[0] else _AD
            if sample.sample == nt_name[0]:
                AD_N = _AD
            else:
                AD_T = _AD
        for idx, para in enumerate(para_list):
            if idx == 0 and para != '*':    # DP_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['DP'] <= para:
                        return False
            elif idx == 1 and para != '*':    # DP_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['DP'] <= para:
                        return False
            elif idx == 2 and para != '*':    # AD_N
                for sample in record.samples:
                    if sample.sample == nt_name[0]:
                        if AD_N[1] < para:
                            return False
            elif idx == 3 and para != '*':    # AD_T
                for sample in record.samples:
                    if sample.sample == nt_name[1]:
                        if AD_T[1] < para:
                            return False
            elif idx == 4 and para != '*':    # AF_N
                for sample in record.samples:
                    if sample.sample == nt_name[0]:
                        AF_N = round(AD_N[1]/sample['DP'],3)
                        if AF_N > para:
                            return False
            elif idx == 5 and para != '*':    # AF_T
                for sample in record.samples:
                    if sample.sample == nt_name[1]:
                        AF_T = round(AD_T[1]/sample['DP'],3)
                        if AF_T < para:
                            return False
    elif call == 'VarScan2':
        for idx, para in enumerate(para_list):
            if idx == 0 and para != '*':    # DP_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['DP'] <= para:
                        return False
            elif idx == 1 and para != '*':    # DP_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['DP'] <= para:
                        return False
            elif idx == 2 and para != '*':    # AD_N
                for sample in record.samples:
                    if sample.sample == nt_name[0] and sample['AD'] > para:
                        return False
            elif idx == 3 and para != '*':    # AD_T
                for sample in record.samples:
                    if sample.sample == nt_name[1] and sample['AD'] < para:
                        return False
            elif idx == 4 and para != '*':    # AF_N
                for sample in record.samples:
                    if sample.sample == nt_name[0]:
                        AF_N = round(sample['AD']/sample['DP'],3)
                        if AF_N > para:
                            return False
            elif idx == 5 and para != '*':    # AF_T
                for sample in record.samples:
                    if sample.sample == nt_name[1]:
                        AF_T = round(sample['AD']/sample['DP'],3)
                        if AF_T < para:
                            return False
    else:
        raise ValueError('[MutScape] The name for caller variant is not defined.')
    return True
    

def pass_filter(record):
    ''' PASS(PA) filter 
    Return  `False` if `record.FILTER` is not `PASS`

    Parameters
    ----------
    record : vcf.model._Record
        A example of `record` is `Record(CHROM=1, POS=1560973, REF=C, ALT=[T])`.
    
    Returns
    -------
    PASS : bool
    '''
    PASS = True
    if type(record.FILTER) == list:
        PASS = True if len(record.FILTER) == 0 else False
    elif type(record.FILTER) == type(None):
        PASS = True
    return PASS

def artifact_variant(call,record,delta):
    ''' Artifact variant(AV) filter

    Parameters
    ----------
    call : str
    record : vcf.model._Record
        A example of `record` is `Record(CHROM=1, POS=1560973, REF=C, ALT=[T])`.
    delta : float

    Returns
    -------
    True/False : bool

    Raises
    ------
    ValueError if the name of variant caller is not defined.
    '''
    if call == 'Mutect2' or call == 'Dragen':
        F1R2 = record.samples[1]['F1R2'][1]
        F2R1 = record.samples[1]['F2R1'][1]
        d = abs((F1R2-F2R1)/(F1R2+F2R1))
        return True if d <= delta else False
    else:
        raise ValueError('[MutScape] The name of variant caller is not defined.')

def vcf_filter(if_filter, category, category_caller, meta):
    '''Implement VCF formalizing and filtering simultaneously.

    Parameters
    ----------
    if_filter : bool
    category : list
    category_caller : list
    meta : str

    Returns
    -------
    category : list
        New category must be updated.
    
    Outputs
    -------
    The formalized VCFs or filtered VCFs will be outputed in the path of `meta`.
    '''
    filter_list = []
    if if_filter:
        print(colored("\nFormalizing and filtering VCF files....\n", "yellow"))
        filter_list = get_filter_data(if_filter)
    else:
        print(colored("\nFormalizing VCF files....\n", "yellow"))
    
    for s_idx, sample in enumerate(category):
        formalized_data_list = []
        vcf_read_set = [read_vcf(x) for x in sample[2]]
        for vcf_idx, vcf_read in enumerate(vcf_read_set):
            new_filter_list = list(filter_list)
            if len(filter_list) != 0:
                FFPE = (category_caller[s_idx][vcf_idx] != 'Mutect2' or category_caller[s_idx][vcf_idx] != 'Dragen')
                if FFPE and new_filter_list[3]:
                    new_filter_list[3] = None
                    print(colored(("Warning: FFPE filter does not apply to variant = "+ category_caller[s_idx][vcf_idx]), 'yellow'))
                    print("Skip the FFPE filter for " + sample[2][vcf_idx])
            del_count, change = 0, False
            output_path = meta+(sample[2][vcf_idx].split("/")[-1])[:-4]+"_formalized.vcf" if not if_filter else meta+(sample[2][vcf_idx].split("/")[-1])[:-4]+"_formalized_filter.vcf"
            formalized_data_list.append(output_path)
            if len(vcf_read.samples) == 2:
                change = nt_test(vcf_read, len(vcf_read.samples), category_caller[s_idx][vcf_idx])
            vcf_read = read_vcf(sample[2][vcf_idx])
            vcf_read.samples = [sample[0], sample[1]] if len(vcf_read.samples) == 2 else [sample[1]]
            vcf_writer = vcf.Writer(open(output_path,"w"), vcf_read)
            for record in vcf_read:
                PASS = True
                if change and PASS:
                    x = record.samples[0]
                    record.samples[0] = record.samples[1]
                    record.samples[1] = x
                if "chr" in record.CHROM and PASS:
                    record.CHROM = record.CHROM[3:]
                if len(record.ALT) != 1 and PASS:
                    del_count += 1
                    PASS = False
                
                if len(filter_list) != 0 and PASS:
                    filter_score = np.zeros(4, dtype = bool)
                    for i in range(len(new_filter_list)):
                        nt_name = [category[0][0], category[0][1]]
                        if PASS:
                            if new_filter_list[i]:
                                if i == 0:
                                    PASS = genome_interval(record, new_filter_list[i])
                                elif i == 1:
                                    call = category_caller[s_idx][vcf_idx]
                                    [DP_N, DP_T, AD_N, AD_T, AF_N, AF_T, NLOD, TLOD] = new_filter_list[i]
                                    PASS = caller_info(call, record, new_filter_list[i], nt_name)
                                elif i == 2:
                                    PASS = pass_filter(record)
                                elif i == 3:
                                    call = category_caller[s_idx][vcf_idx]
                                    PASS = True if call != 'Mutect2' or call != 'Dragen' else artifact_variant(call, record, new_filter_list[i])
                        else:
                            break
                    del_count = del_count + 1 if not PASS else del_count
                if PASS:
                    vcf_writer.write_record(record)
            if del_count != 0:
                print("NOTICE: " + str(del_count)+ " data have been removed from "+ sample[2][vcf_idx])
            vcf_writer.close()
        listToStr = ', '.join([str(elem) for elem in formalized_data_list]) 
        sample[2] = formalized_data_list
        print(colored("\n=> Finish formalizing files: \n[ "+ listToStr+" ]\n", "green"))
    return category