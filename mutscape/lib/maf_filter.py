############################################################################################
# FileName     [ maf_filter.py ]
# PackageName  [ lib ]
# Synopsis     [ Sifting records from MAFs ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

import csv
from termcolor import colored
import pandas as pd
from tqdm import tqdm, trange
import os

def get_maf_filter_data(maf_flt, num = 5):
    '''Interpret the input filtering parameters.

    Parameters
    ----------
    maf_flt : list
        The maximum size of `maf_flt` is 10 and ths size of `maf_flt`
        must be even.
    num : int, optional
        The number of filters provided in this function.
    
    Returns
    -------
    flt_list : list
        The size of `flt_list` is 5. If no parameter is entered, the 
        item will be recorder to `False`.

    Raises
    ------
    ValueError if the format for command -mf is wrong.
    '''
    if len(maf_flt)%2 != 0:
        raise ValueError('[MutScape] Wrong format for command -mf!')
    flt_list = [False]*num
    for idx, data in enumerate(maf_flt):
        if idx%2 == 0:
            if data == "GI":
                info = maf_flt[idx+1]
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
                        s, e = int(s.strip('[')), int(e.strip(']'))
                        new_dict[str(i[0])] = [s, e]
                    flt_list[0] = new_dict
            elif data == "CI":
                info = maf_flt[idx+1]
                if "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    new_info = []
                    for i in info:
                        new_info.append(int(i))
                    flt_list[1] = new_info
            elif data == "TE":
                info = maf_flt[idx+1]
                if "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    flt_list[2] = info
            elif data == "PF":
                flt_list[3] = bool(int(maf_flt[idx+1]))
            elif data == "H":
                flt_list[4] = int(maf_flt[idx+1])
    return flt_list

def read_TSV(tsv_file):
    '''Read lib/auxiliary/rna_tissue_consensus.tsv

    Parameters
    ----------
    tsv_file : str
        The path of `rna_tissue_consensus.tsv`.

    Returns
    -------
    ALL_DICT : dict

    '''
    print(colored("Reading TSV auxiliary file....\nIt may take a while....\n", "yellow"))
    ALL_DICT = {}
    df = pd.read_csv(tsv_file, sep="\t", header = 0)
    pbar = tqdm(total = df.shape[0])
    for i in range(df.shape[0]):
        pbar.update(1)
        if df.iloc[i]['Gene name'] not in ALL_DICT:
            ALL_DICT[df.iloc[i]['Gene name']] = {}
        ALL_DICT[df.iloc[i]['Gene name']][df.iloc[i]['Tissue']] = df.iloc[i]['NX']
    pbar.close()
    return ALL_DICT

def maf_filter(maf_file, flt_list, ALL_DICT, output_file):
    ''' Write a new MAF file to output_file

    Parameters
    ----------
    maf_file : str
    flt_list : list
        All the conditions and parameters of the 5 filters.
    ALL_DICT : dict
    output_file : str

    '''
    def fast_read_maf(maf_file):
        '''Read MAF(Mutation Annotation Format) file

        Parameters
        ----------
        maf_file : str

        Returns
        -------
        head : str
        df : pandas.dataFrame
        '''
        file = open(maf_file, 'r')
        head = file.readline()
        if head.startswith('#'):
            df = pd.read_csv(maf_file, sep="\t", skiprows = 1, header = 0,encoding = 'gb2312', low_memory=False)
        else:
            df = pd.read_csv(maf_file, sep="\t", header = 0, encoding = 'gb2312', low_memory=False)
            head = ""
        return head, df

    def genome_interval(data, interval):
        '''Genome Interval (GI) filter

        Parameters
        ----------
        data :
        interval : bool/list/dict

        Returns
        -------
        True/False : bool
        '''
        if interval == False:
            return True
        if isinstance(interval, list):  # Only select Genome Number
            if data['Chromosome'] in interval:
                return True
            else:
                return False
        elif isinstance(interval, dict):
            chromo = str(data['Chromosome'])
            start, end = data['Start_Position'], data['End_Position']
            if chromo in interval:
                i_start, i_end = interval[chromo][0], interval[chromo][1]
                if (start >= i_start and start <= i_end) or (end >= i_start and end <= i_end):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def caller_info(data, info):
        '''Caller Information (CI) filter

        Parameters
        ----------
        data :
        info : bool/list

        Returns
        -------
        True/False : bool
        '''
        if info == False:
            return True
        DP_N, DP_T, AD_T = info[0], info[1], info[2]
        if data['t_depth'] > DP_T and data['t_alt_count'] >= AD_T and data['n_depth'] > DP_N and data['n_alt_count'] == '.':
            return True
        else:
            return False
    
    def tissue_expression(data, tissue, ALL_DICT):
        '''Tissue Expression (TE) filter

        Parameters
        ----------
        data :
        tissue : list
        ALL_DICT : dict

        Returns
        -------
        True/False : bool
        '''
        if tissue == False:
            return True
        Hugo = data['Hugo_Symbol']
        PASS = True
        if Hugo in ALL_DICT:
            for item in tissue:
                if item in ALL_DICT[Hugo]:
                    if float(ALL_DICT[Hugo][item]) <= 0:
                        PASS = False
                else:
                    PASS = False
        else:
            PASS = False
        return PASS
    
    def population_frequency(data, info):
        '''Population Frequency (PF) filter

        Parameters
        ----------
        data :
        info : bool

        Returns
        -------
        True/False : bool
        '''
        if info == False:
            return True
        return data['FILTER'] == 'PASS'
    
    def Hypermutator(df, COUNT):
        '''Hypermutator (HY) filter

        Parameters
        ----------
        df : pandas.dataFrame
        COUNT : int

        Returns
        -------
        df : pandas.dataFrame
        '''
        sample_dict = {}
        rm_list = []
        for i in range(df.shape[0]):
            data = df.iloc[i]
            if data['Tumor_Sample_Barcode'] not in sample_dict:
                sample_dict[data['Tumor_Sample_Barcode']] = [i]
            else:
                sample_dict[data['Tumor_Sample_Barcode']].append(i)
        for item in sample_dict:
            if len(sample_dict[item]) > COUNT:
                rm_list += sample_dict[item]
        df = df.drop(df.index([rm_list]))
        return df

    head, df = fast_read_maf(maf_file)
    rm_list = []
    pbar = tqdm(total = df.shape[0])
    for i in range(df.shape[0]):
        pbar.update(1)
        GI = genome_interval(df.iloc[i], flt_list[0])
        CI = caller_info(df.iloc[i], flt_list[1])
        TE = tissue_expression(df.iloc[i], flt_list[2], ALL_DICT)
        PF = population_frequency(df.iloc[i], flt_list[3])
        if not (GI and CI and TE and PF):
            rm_list.append(i)
    pbar.close()
    df = df.drop(df.index[rm_list])
    if flt_list[4] != False:
        df  = Hypermutator(df, flt_list[4])
    df.to_csv(output_file, sep="\t", index= False, header = list(df.columns.values))
    if head != "":
        with open(output_file, "r+") as filtered_file:
            lines = filtered_file.readlines()
            filtered_file.seek(0)
            filtered_file.write(head)
            filtered_file.writelines(lines)

def all_maf_filter(input_params_list, maf_output_list):
    ''' Implement MAF filtering and combination

    Parameters
    ----------
    input_params_list : list
    maf_output_list : list

    Returns
    -------
    '''
    if input_params_list:
        print(colored("Start MAF filtering....\n", "yellow"))
        maf_flt_list = get_maf_filter_data(input_params_list)
        maf_filtered_list = [x[:-4]+"_filtered.maf" for x in maf_output_list]
        ALL_DICT = {}
        if maf_flt_list[2] != False:
            if not os.path.isfile('lib/auxiliary/rna_tissue_consensus.json'):
                ALL_DICT = read_TSV("lib/auxiliary/rna_tissue_consensus.tsv")
                with open("lib/auxiliary/rna_tissue_consensus.json", "w") as jsonfile:  
                    json.dump(ALL_DICT, jsonfile) 
            else:
                with open("lib/auxiliary/rna_tissue_consensus.json", "r") as jsonfile:  
                    ALL_DICT = json.load(jsonfile)
        for idx, maf in enumerate(maf_output_list):
            maf_filter(maf, maf_flt_list, ALL_DICT, maf_filtered_list[idx])
            print(colored(("\n=> Finish filtering MAF file: "+maf_filtered_list[idx]+"\n"), 'green'))
        # MAF combination
        if len(maf_filtered_list) > 1:
            print(colored("Start MAF combination....\n", "yellow"))
            maf_df, head = pd.DataFrame(), ""
            for maf_file in maf_filtered_list:
                head, maf = fast_read_maf(maf_file)
                maf_df = pd.concat([maf_df, maf])
            maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
            if head != "":
                with open(folder+"maf_combination.maf", "r+") as filtered_file:
                    lines = filtered_file.readlines()
                    filtered_file.seek(0)
                    filtered_file.write(head)
                    filtered_file.writelines(lines)
            print(colored(("=> Finish combining MAF files to " + folder + "maf_combination.maf" + "\n"), 'green'))
        else:
            if len(category) > 1:
                print(colored("Start MAF combination....\n", "yellow"))
                maf_df, head = pd.DataFrame(), ""
                for maf_file in category:
                    head, maf = fast_read_maf(maf_file)
                    maf_df = pd.concat([maf_df, maf])
                maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
                if head != "":
                    with open(folder+"maf_combination.maf", "r+") as filtered_file:
                        lines = filtered_file.readlines()
                        filtered_file.seek(0)
                        filtered_file.write(head)
                        filtered_file.writelines(lines)
                print(colored(("=> Finish combining MAF files to " + folder + "maf_combination.maf" + "\n"), 'green'))