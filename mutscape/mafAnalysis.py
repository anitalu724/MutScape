############################################################################################
# FileName     [ mafAnalysis.py ]
# PackageName  [ MutScape ]
# Synopsis     [ Control all functions to implement MAF analysis and visualization. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

import argparse, textwrap, ast, os
from lib.analysis.sig_mutated_gene_detect import SigMutatedGeneDetection
from lib.analysis.known_cancer_gene_anno import KnownCancerGeneAnnotation
from lib.analysis.tumor_mutated_burden import TumorMutationBurden
from lib.analysis.comut_plot_analysis import CoMutAnalysis, CoMutPlot
from lib.analysis.mutational_sig import MutationalSignature
from lib.analysis.hrd_score import HRDScore, HCWCompare
from lib.analysis.wgd_cin import WGDnCIN
from lib.analysis.oncokb_annotation import OncoKBAnnotator

def main():
    ''' Implement MAF analysis and visualization in one single command.

    This function contains 8 different analyses and some of them generate
    plots after analysis:

    1. Significantly mutated gene detection
    2. Known cancer gene annotation
    3. Mutation burden statistics
    4. CoMut plot analysis
    5. Mutational signature
    6. HRD Score
    7. Whole-genome doubling (WGD) and Chromosome instability (CIN)
    8. HRD_CIN_WGD Comparison
    9. Actionable mutation(drug) annotation

    '''
    parser = argparse.ArgumentParser(description="MAF analysis", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file", nargs=1, metavar="MAF file")
    parser.add_argument("-smg", "--significantly_mutated_gene", action="store_true")
    parser.add_argument("-kcga", "--known_cancer_gene_annotaiton", action="store_true")
    parser.add_argument("-tmb","--tumor_mutation_burden",nargs=1,help="One item must be entered:\n \
                                                                       1. Sequencing Length\n",)
    parser.add_argument("-cm", "--comut_analysis", nargs=1,help="One item must be entered:\n \
                                                                       1. Sequencing Length\n",)
    parser.add_argument('-cmp', '--comut_plot',nargs=4, help="Four items need to be entered:\n\
                                                              1. TSV file for data paths.\n\
                                                              2. TSV file which contain all information for image.\n\
                                                              3. plot color theme(0: cold, 1: warm)\n\
                                                              4. CoMut_plot picture file name")
    parser.add_argument("-ms", "--mutational_signature", nargs=2, metavar=('step_#', 'list_of_params'),
                        help=textwrap.dedent('''\
                        Two steps are included: Estimation(0) and Plotting(1).
                        For each step, two values are required:
                        * Estimation
                          1. 1
                          2. [<Rank start number>, <Rank end number>, <Epoch number>]
                        * Plotting
                          1. 2
                          2. [Signature Number]'''))
    parser.add_argument("-hrd","--hrd_score",nargs=2,help="Two items must be entered:\n\
                                                           1. The CSV_input file.\n\
                                                           2. The reference for HRD Score.\n")

    parser.add_argument("-hcwc","--hrd_cin_wgd_comparison",nargs=2,help="Two items must be entered:\n\
                                                           1. The CSV_input file1.\n\
                                                           2. The CSV_input file2.\n\
                                                           3. The reference for HRD Score.\n")
    parser.add_argument("-wgdcin", "--wgd_cin", nargs=1)
    parser.add_argument("-oncokb","--oncokb_annotator",nargs='*',help='Three items must be entered:\n \
                                                                     1. The relative path of the folder "oncokb-annotator".\n \
                                                                     2. The token of your OncoKB account.\n \
                                                                     3. The level of the drug. \n\
                                                                     4. The file path of clinical input.\n\
                                                                     5. (Optional) The file path of cna input.\n')

    parser.add_argument("-o","--output",required=True,metavar="OUTPUT folder",help="The path for storing every generated file.\n\
                                                                                    This path must end with a folder.\n")
    parser.add_argument("-p", "--picture", required=True,metavar="picture folder" ,help="The path for storing every picture.\n")

    args = parser.parse_args()

    folder = args.output if args.output[-1:] == '/' else (args.output + '/')
    pic = args.picture if args.picture[-1:] == '/' else (args.picture + '/')

    if args.significantly_mutated_gene:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -significantly_mutated_gene needs -file as a parameter.')
        df = SigMutatedGeneDetection(args.file[0])
        df.oncodriveCLUST(folder)
    
    if args.known_cancer_gene_annotaiton:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -known_cancer_gene_annotaiton needs -file as a parameter.')
        df = KnownCancerGeneAnnotation(args.file[0])
        df.annotation(folder)
    
    if args.tumor_mutation_burden:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -tumor_mutation_burden needs -file as a parameter.')
        df = TumorMutationBurden(args.file[0])
        df.data_analysis(folder, int(args.tumor_mutation_burden[0]))
    
    if args.comut_analysis:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -comut_analysis needs -file as a parameter.')
        df = CoMutAnalysis(args.file[0])
        df.data_analysis(folder, int(args.comut_analysis[0]))
    
    if args.comut_plot:
        plot1 = CoMutPlot(args.comut_plot[0], args.comut_plot[1])
        plot1.plot(pic, args.comut_plot[2], args.comut_plot[3])
    
    if args.mutational_signature:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -mutational_signature needs -file as a parameter.')
        df = MutationalSignature(args.file[0])
        # if args.mutational_signature[0] != '0':
        # params = ast.literal_eval(args.mutational_signature[1])
        if args.mutational_signature[0] != '2':
            df.get_input_file(folder)
        if args.mutational_signature[0] != '0':
            params = ast.literal_eval(args.mutational_signature[1])

        if args.mutational_signature[0] == '0':
            df.sig_refitting()
            df.getParams(args.mutational_signature[1])
            df.SBSplot(df.cosmic, pic)
            df.SigDistribution(df.contribution, folder, pic)
            df.DonutPlot(df.contribution, pic)
        elif args.mutational_signature[0] == '1':
            df.estimation(folder, pic, params[0], params[1], params[2])
        elif args.mutational_signature[0] == '2':
            df.plotting(folder, pic, params[0])
        else:
            raise ValueError('[MutScape] Command -ms only support 0/1/2 parameters.')
        
    if args.hrd_score:
        df = HRDScore(args.hrd_score[0])
        df.data_analysis(folder, args.hrd_score[1])
        df.plotting(folder, pic)

    if args.hrd_cin_wgd_comparison:
        df = HCWCompare(args.hrd_cin_wgd_comparison[0])
        for idx, fileList in enumerate(df.fileList):
            df.HRD(idx, fileList, folder, args.hrd_cin_wgd_comparison[1])
            df.WGD_CIN(idx, fileList, folder)
        df.CINbarplot(pic)
        df.HRDbarplot(pic)
        df.WGDheatmap(pic)
        df.HRDheatmap(pic)

    if args.wgd_cin:
        df = WGDnCIN(args.wgd_cin[0])
        df.data_analysis(folder)
        df.plotting(folder, pic)
    
    if args.oncokb_annotator:
        if args.file == None or len(args.file) == 0:
            raise ValueError('[MutScape] Command -oncokb_annotator needs -file as a parameter.')
        df = OncoKBAnnotator(args.file[0])
        if len(args.oncokb_annotator) == 5:
            df.data_analysis(folder,args.oncokb_annotator[0],args.oncokb_annotator[1],args.oncokb_annotator[3],args.oncokb_annotator[4])
        else:
            df.data_analysis(folder,args.oncokb_annotator[0],args.oncokb_annotator[1],args.oncokb_annotator[3])
        df.plotting(folder,pic, args.oncokb_annotator[2])

if __name__ == '__main__':
    main()