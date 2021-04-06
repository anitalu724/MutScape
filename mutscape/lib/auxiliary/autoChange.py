import os
mypath = "./"

if 'AnnotatorCore.py' in os.listdir(mypath):
    rfile = open('AnnotatorCore.py', 'r')
    lines = rfile.readlines()
    if 'def getsampleid(rawsampleid):\n' in lines:
        ifunc = lines.index('def getsampleid(rawsampleid):\n')
        if '    if rawsampleid.startswith("TCGA"):\n' == lines[ifunc+1]:
            lines[ifunc+1] = '#    if rawsampleid.startswith("TCGA"):\n'
        if '        return rawsampleid[:15]\n' == lines[ifunc+2]:
            lines[ifunc+2] = '#        return rawsampleid[:15]\n'

with open('AnnotatorCore.py', 'w') as wfile:
    wfile.writelines( lines )
       
    