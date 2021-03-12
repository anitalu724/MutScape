############################################################################################
# FileName     [ vcf_filter.py ]
# PackageName  [ lib ]
# Synopsis     [ Sifting records from VCFs ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 3 ]
############################################################################################

import vcf
from .vcf_tool import *
from termcolor import colored

#Genome Interval
class GenemoInterval:
