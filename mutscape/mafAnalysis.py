############################################################################################
# FileName     [ mafAnalysis.py ]
# PackageName  [ MutScape ]
# Synopsis     [ Control all functions to implement MAF analysis and visualization. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

import argparse, textwrap, ast
from lib.analysis.sig_mutated_gene_detect import SigMutatedGeneDetection
from lib.analysis.known_cancer_gene_anno import KnownCancerGeneAnnotation
from lib.analysis.total_mutated_burden import TotalMutationBurden
from lib.analysis.comut_plot_analysis import CoMutAnalysis, CoMutPlot
from lib.analysis.mutational_sig import MutationalSignature

def main():
    ''' Implement MAF analysis and visualization in one single command.

    This function contains 8 different analyses and some of them generate
    plots after analysis.

    Examples
    --------
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -smg \
    -o examples/output \
    -p examples/pic/




    '''
    parser = argparse.ArgumentParser(description="MAF analysis", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file", nargs=1, metavar="MAF file", required=True)
    parser.add_argument("-smg", "--significantly_mutated_gene", action="store_true")
    parser.add_argument("-kcga", "--known_cancer_gene_annotaiton", action="store_true")
    parser.add_argument("-tmb","--total_mutation_burden",nargs=1,help="One item must be entered:\n \
                                                                       1. Sequencing Length\n",)
    parser.add_argument("-cm", "--comut_analysis", action="store_true")
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

    parser.add_argument("-o","--output",required=True,metavar="OUTPUT folder",help="The path for storing every generated file.\n\
                                                                                    This path must end with a folder.\n")
    parser.add_argument("-p", "--picture", required=True,metavar="picture folder" ,help="The path for storing every picture.\n")

    args = parser.parse_args()

    folder = args.output if args.output[-1:] == '/' else (args.output + '/')
    pic = args.picture if args.picture[-1:] == '/' else (args.picture + '/')

    if args.significantly_mutated_gene:
        df = SigMutatedGeneDetection(args.file[0])
        df.oncodriveCLUST(folder)
    if args.known_cancer_gene_annotaiton:
        df = KnownCancerGeneAnnotation(args.file[0])
        df.annotation(folder)
    if args.total_mutation_burden:
        df = TotalMutationBurden(args.file[0])
        df.data_analysis(folder, int(args.total_mutation_burden[0]))
    if args.comut_analysis:
        df = CoMutAnalysis(args.file[0])
        df.data_analysis(folder)
    if args.comut_plot:
        plot1 = CoMutPlot(args.comut_plot[0], args.comut_plot[1])
        plot1.plot(pic, args.comut_plot[2], args.comut_plot[3])
    if args.mutational_signature:
        df = MutationalSignature(args.file[0])
        params = ast.literal_eval(args.mutational_signature[1])
        if args.mutational_signature[0] == '1':
            df.data_analysis(folder, pic, params[0], params[1], params[2])
        elif args.mutational_signature[0] == '2':
            df.plotting(folder, pic, params[0])

if __name__ == '__main__':
    main()