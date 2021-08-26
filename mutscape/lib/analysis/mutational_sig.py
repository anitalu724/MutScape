############################################################################################
# FileName     [ mutational_sig.py ]
# PackageName  [ lib/analysis ]
# Synopsis     [ Implement mutational signature analysis. ]
# Author       [ Cheng-Hua Lu ]
# Copyright    [ 2021 4 ]
############################################################################################

from numpy.core.numeric import outer
from ..maf_filter import fast_read_maf
from termcolor import colored
import pandas as pd
import numpy as np
import math
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.ticker as mtick
import matplotlib.style
import matplotlib

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']

class MutationalSignature:
    '''MAF analysis: Mutational signature

    Parameters
    ----------
    maf_file : str
        A MAF file path.
    output_folder : str
        The path for every output file.
    pic : str
        The path for storing comut plot.
    rank1, rank2 : int
        The range for estimate # signature.
    epoch : int
        # estimation running.
    sig : int
        The final factorization rank(# signature)
    
    Output files
    ------------
    output :
        ms_input.tsv
        96_sig.csv
        sig_sample.csv
        SBS.tsv
        
    pictures:
        Estimation.pdf
        SBS_96_plots.pdf
        S2S.pdf
        SigContribution.pdf
        SigSamHeatmap.pdf
        Donut_plot.pdf

    '''
    def __init__(self, maf_file):
        print(colored(('\nStart Mutational_Signature....'), 'yellow'))
        self.head, self.df = fast_read_maf(maf_file)
    def data_analysis(self, output_folder, pic, rank1, rank2, epoch):
        def get_input_file():
            output_file = output_folder+'ms_input.tsv'
            selected_col = self.df[['Tumor_Sample_Barcode','flanking_bps', 'Reference_Allele', 'Tumor_Seq_Allele2']]
            selected_col.columns = ['SampleID', 'Three_Allele', 'Ref', 'Mut']
            sample_list = selected_col.SampleID.unique()
            grouped = selected_col.groupby(selected_col['SampleID'])
            df_list = [grouped.get_group(sample).reset_index(drop=True) for sample in sample_list]
            final_dict = {}
            for d, df in enumerate(df_list):
                # order: 'C>A','C>G','C>T','T>A','T>C','T>G'
                cata_list = [[],[],[],[],[],[]]
                for i in range(len(df)):
                    item = df.loc[i]
                    if (item['Ref'] == 'C' and item['Mut'] == 'A') or (item['Ref'] == 'G' and item['Mut'] == 'T'):
                        cata_list[0].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'G') or (item['Ref'] == 'G' and item['Mut'] == 'C'):
                        cata_list[1].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'T') or (item['Ref'] == 'G' and item['Mut'] == 'A'):
                        cata_list[2].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'A') or (item['Ref'] == 'A' and item['Mut'] == 'T'):
                        cata_list[3].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'C') or (item['Ref'] == 'A' and item['Mut'] == 'G'):                       
                        cata_list[4].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'G') or (item['Ref'] == 'A' and item['Mut'] == 'C'):
                        cata_list[5].append(item)
                list_96 = []
                for cata in range(len(cata_list)):
                    cata_sum_list = [int(0)]*16
                    if cata in [0,1,2]:
                        three_allele_dict={'ACA':0,     'TGT':0,    'ACC':1,    'GGT':1,    'ACG':2,    'CGT':2,    'ACT':3,    'AGT':3, \
                                           'CCA':4,     'TGG':4,    'CCC':5,    'GGG':5,    'CCG':6,    'CGG':6,    'CCT':7,    'AGG':7, \
                                           'GCA':8,     'TGC':8,    'GCC':9,    'GGC':9,    'GCG':10,   'CGC':10,   'GCT':11,   'AGC':11,\
                                           'TCA':12,    'TGA':12,   'TCC':13,   'GGA':13,   'TCG':14,   'CGA':14,   'TCT':15,   'AGA':15 }   
                    elif cata in [3,4,5]:
                        three_allele_dict={'ATA':0,     'TAT':0,    'ATC':1,    'GAT':1,    'ATG':2,    'CAT':2,    'ATT':3,    'AAT':3, \
                                           'CTA':4,     'TAG':4,    'CTC':5,    'GAG':5,    'CTG':6,    'CAG':6,    'CTT':7,    'AAG':7, \
                                           'GTA':8,     'TAC':8,    'GTC':9,    'GAC':9,    'GTG':10,   'CAC':10,   'GTT':11,   'AAC':11,\
                                           'TTA':12,    'TAA':12,   'TTC':13,   'GAA':13,   'TTG':14,   'CAA':14,   'TTT':15,   'AAA':15 }  

                    for j in range(len(cata_list[cata])):
                        if (cata_list[cata][j])['Three_Allele'] in three_allele_dict:
                            cata_sum_list[three_allele_dict[(cata_list[cata][j])['Three_Allele']]] += 1;
                    list_96 += cata_sum_list
                final_dict[sample_list[d]] = list_96

            new_df = pd.DataFrame.from_dict(final_dict)
            list_a = ['A.A', 'A.C', 'A.G', 'A.T', 'C.A', 'C.C', 'C.G', 'C.T',\
                      'G.A', 'G.C', 'G.G', 'G.T', 'T.A', 'T.C', 'T.G', 'T.T']
            list_b = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
            new_row_name = []
            for item in list_b:
                for allele in list_a:
                    new_str = allele[0]+'['+item+']'+allele[2]
                    new_row_name.append(new_str)
            new_df.index = new_row_name
            new_df.to_csv(output_file, sep = '\t', index = True)
            print(colored('=> Generate input file: ', 'green'))
            print(colored(('   '+output_file), 'green'))
        def estimation():
            os.system('git clone https://github.com/mims-harvard/nimfa.git\n')
            os.chdir('nimfa')
            os.system('python3 setup.py install --user')
            code = open('nimfa.py', 'w')
            code.write("import nimfa\nfrom collections import defaultdict, Counter\nimport urllib\nimport numpy as np\nfrom matplotlib import pyplot as plt\nimport matplotlib.gridspec as gridspec\nfrom sklearn import preprocessing\nimport scipy.cluster.hierarchy as sch\nimport pandas as pd\n")
            code.write("df = (pd.read_csv(\"../" + output_folder + "ms_input.tsv\", sep=\"\t\")).T\n")
            code.write("data = (df.to_numpy())[1:]\n")
            code.write("rank_cands = range("+str(rank1)+","+ str(rank2)+", 1)\n")
            code.write("snmf = nimfa.Snmf(data, seed='random_vcol', max_iter=100)\n")
            code.write("summary = snmf.estimate_rank(rank_range=rank_cands, n_run="+str(epoch)+", what='all')\n")
            code.write("rss = [summary[rank]['rss'] for rank in rank_cands]\n")
            code.write("coph = [summary[rank]['cophenetic'] for rank in rank_cands]\n")
            code.write("disp = [summary[rank]['dispersion'] for rank in rank_cands]\n")
            code.write("spar = [summary[rank]['sparseness'] for rank in rank_cands]\n")
            code.write("spar_w, spar_h = zip(*spar)\n")
            code.write("evar = [summary[rank]['evar'] for rank in rank_cands]\n")
            code.write("fig, axs = plt.subplots(2, 3, figsize=(12,8))\n")
            code.write("axs[0,0].plot(rank_cands, rss, 'o-', color='#266199', label='RSS', linewidth=3)\n")
            code.write("axs[0,0].set_title('RSS', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,0].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,0].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[0,1].plot(rank_cands, coph, 'o-', color='#695D73', label='Cophenetic correlation', linewidth=3)\n")
            code.write("axs[0,1].set_title('Cophenetic', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,1].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,1].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[0,2].plot(rank_cands, disp,'o-', color='#71a0a5', label='Dispersion', linewidth=3)\n")
            code.write("axs[0,2].set_title('Dispersion', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,2].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,2].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,0].plot(rank_cands, spar_w, 'o-', color='#B88655', label='Sparsity (Basis)', linewidth=3)\n")
            code.write("axs[1,0].set_title('Sparsity (Basis)', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,0].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,0].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,1].plot(rank_cands, spar_h, 'o-', color='#E08B69', label='Sparsity (Mixture)', linewidth=3)\n")
            code.write("axs[1,1].set_title('Sparsity (Mixture)', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,1].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,1].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,2].plot(rank_cands, evar,  'o-', color='#841D22', label='Explained variance', linewidth=3)\n")
            code.write("axs[1,2].set_title('Explained variance', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,2].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,2].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("fig.tight_layout(pad=1.0)\n")
            code.write("plt.savefig(\"../"+pic+"Estimation.pdf\",dpi=300,bbox_inches = 'tight')\n")
            code.close()
            print(colored(('\nStart Estimation (may need a few minutes)....'), 'yellow'))
            p = os.popen('python3 nimfa.py\n')
            x = p.read()
            print(x)
            p.close()
            print(colored('=> Generate estimation figure: ', 'green'))
            print(colored(('   '+pic+'Estimation.pdf\n'), 'green'))
            os.chdir('..')
            os.system('rm -rf nimfa\n')
        get_input_file()
        estimation()  
    def plotting(self, output_folder, pic, sig):
        LABEL_SIZE, TITLE_SIZE = 24,30
        print(colored(('\nStart Mutational_Signature Plotting(signature number must be in the range of 2 to 9)....'), 'yellow'))
        def nmf():
            print(colored(('\nStart NMF....'), 'yellow'))
            from sklearn.decomposition import NMF
            if not os.path.isfile(output_folder+'ms_input.tsv'):
                raise ValueError('[MutScape] Mutational Signature: Step 1 must be done before step 2.')
            df = (pd.read_csv(output_folder+'ms_input.tsv', sep='\t')).T
            sample_list = df.index[1:]
            index_96 = df.to_numpy()[0]
            data = (df.to_numpy())[1:]
            model = NMF(n_components=int(sig),init='random', random_state=0)
            W = model.fit_transform(data)
            H = model.components_
            Hdf, Wdf = pd.DataFrame(H.T), pd.DataFrame(W.T)
            Hdf.columns = ['Signature '+str(i+1) for i in range(int(sig))]
            Wdf.columns = sample_list
            Hdf.index = index_96
            Wdf.index = ['Signature '+str(i+1) for i in range(int(sig))]
            Hdf.to_csv(output_folder+'96_sig.csv')
            Wdf.to_csv(output_folder+'sig_sample.csv')
            print(colored('=> Generate file: ', 'green'))
            print(colored(('   '+output_folder+'96_sig.csv'), 'green'))
            print(colored(('   '+output_folder+'sig_sample.csv'), 'green'))
        def SBSPlot():
            df = (pd.read_csv(output_folder+'96_sig.csv'))
            df = df.set_index(list(df.columns[[0]]))
            fig_x = tuple([ ' '+i[0]+' '+i[6] for i in list(df.index)])
            y_pos = np.arange(len(fig_x))
            fig_name = list(df.columns)
            fig, axes = plt.subplots(df.shape[1], 1, figsize=(12,2*df.shape[1]))#
            if df.shape[1] == 1:
                return
            for r in range(df.shape[1]):
                color_set = ['#02bdee', '#010101','#e32925','#cac9c9', '#a1cf63', '#ecc7c4']
                color_96 = [ c for c in color_set for i in range(16)]
                all_data = df.iloc[:, r]
                all_data /= (all_data.sum())
                maximum = max(all_data)*1.25
                data_list = all_data.tolist()
                axes[r].text(0.01, 0.86, fig_name[r], horizontalalignment='left',verticalalignment='center', transform=axes[r].transAxes, fontweight='bold')
                axes[r].bar(y_pos, data_list, color=color_96, width=0.4)
                axes[r].spines['bottom'].set_color('#cac9c9')
                axes[r].spines['top'].set_color('#cac9c9') 
                axes[r].spines['right'].set_color('#cac9c9')
                axes[r].spines['left'].set_color('#cac9c9')
                if r != df.shape[1]-1:
                    axes[r].xaxis.set_visible(False)
                    axes[r].set_xticklabels([])
                axes[r].tick_params(axis='x',length=0)
                axes[r].set_xlim([-0.8,len(data_list)-.8])

                axes[r].tick_params(axis='y',direction='in', color='#cac9c9', labelsize=10)
                axes[r].set_ylabel('Percentage', fontweight='bold')
                axes[r].tick_params(axis='y', labelsize=10)
                axes[r].set_ylim(top = max(all_data)*1.25)
                axes[r].yaxis.set_major_locator(ticker.LinearLocator(5))
                axes[r].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=1))
                for i in range(6):
                    axes[r].add_patch(matplotlib.patches.Rectangle((0+16*i ,maximum*0.95), 15.6 , 0.01, color=color_set[i],transform=axes[r].transData))
            mut_list = ['C>A','C>G','C>T','T>A','T>C','T>G']
            for i in range(6):
                plt.text(0.19+0.13*i,0.916-df.shape[1]*0.0029, mut_list[i], horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, fontweight='bold', fontsize=14)
            plt.xticks(y_pos, fig_x, color='#999999',rotation=90, fontsize=9,horizontalalignment='center',verticalalignment='top',fontname='monospace')#verticalalignment='bottom',
            space = 0.008075
            y_scale = [0.072, 0.084, 0.09, 0.094, 0.097, 0.0987, 0.1, 0.1013, 0.1023]
            for i in range(6):
                for j in range(16):
                    if i < 3:
                        plt.text((0.131+space*16*i)+space*j, y_scale[df.shape[1]-2], 'C',horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, color=color_set[i], fontsize=9, rotation=90,fontname='monospace', fontweight='bold')
                    else:
                        plt.text((0.131+space*16*i)+space*j, y_scale[df.shape[1]-2], 'T',horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, color=color_set[i], fontsize=9, rotation=90,fontname='monospace', fontweight='bold')
            plt.savefig(pic+'SBS_96_plots.pdf',dpi=300, bbox_inches='tight')
            print(colored(('=> Generate SBS Plot: '+pic+'SBS_96_plots.pdf'), 'green'))
        def CosineSimilarity():
            from sklearn.metrics.pairwise import cosine_similarity
            my_file, aux_file = output_folder+'96_sig.csv', 'lib/auxiliary/COSMIC_72.tsv'
            my_df, aux_df = pd.read_csv(my_file, index_col=0), pd.read_csv(aux_file, sep='\t',index_col=0)
            my_list, aux_list = my_df.columns, aux_df.columns
            X = np.array(my_df.T.to_numpy())
            Y = np.array(aux_df.T.to_numpy())
            M = cosine_similarity(X, Y, dense_output=True)
            Mdf= pd.DataFrame(M)
            Mdf.index, Mdf.columns = my_list, aux_list
            Mdf.to_csv(output_folder+'SBS.tsv', sep='\t')
            print(colored('=> Generate file: ', 'green'))
            print(colored(('   '+output_folder+'SBS.tsv'), 'green'))
            
            height, length = len(my_list), len(aux_list)
            sns.set(font_scale=2)
            sns.set_style('white')
            grid_kws = {'height_ratios': (.9, .2),'hspace': 0.3}  
            f, (ax, cbar_ax) = plt.subplots(2,figsize=(20,6), gridspec_kw=grid_kws)
            ax = sns.heatmap(M, vmin=0, vmax=1, xticklabels =aux_list, yticklabels = my_list, square=False, linewidth=1, cbar_ax=cbar_ax,ax=ax,
                                cmap='Blues',cbar_kws={'orientation': 'horizontal','shrink':1, 'aspect':70})
            # ax.set_title('Cosine Similarity',fontsize=TITLE_SIZE,weight='bold',pad=0,verticalalignment='bottom')
            ax.set_xticklabels(ax.get_xticklabels(),rotation=90, horizontalalignment='center', fontsize=LABEL_SIZE-6, color='#222222')
            ax.tick_params(axis='both',length=0)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=LABEL_SIZE-6,color='#222222',verticalalignment='center')
            plt.ylim(bottom=0, top=height+0.5)
            plt.savefig(pic+'S2S.pdf',dpi=300,bbox_inches='tight')
            plt.clf()
            print(colored(('=> Generate Cosine Similarity Plot: '+pic+'S2S.pdf'), 'green'))  
        def SigDistribution():
            df = pd.read_csv(output_folder+'sig_sample.csv', index_col=0)
            sample_list, sig_list = list(df.columns),list(df.index)
            SUM = (df.sum(axis = 0, skipna = True)).tolist()
            df = df/SUM
            dft = df.T
            # dft.columns = ['sample']+dft.columns
            dft.to_csv(output_folder+'SigContribution.tsv',index_label='sample', sep='\t')
            print(colored(('   '+output_folder+'SigContribution.tsv'), 'green'))
            ind = np.arange(df.shape[1])
            data = []
            for i in range(df.shape[0]):
                d = tuple(df.iloc[i].tolist())
                data.append(d)
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_axes([0,0,1,1])
            
            

            for i in range(len(data)):
                if i == 0:
                    ax.bar(ind, data[i], 0.8, color = COLOR_MAP[i])
                else:
                    b = np.array(data[0])
                    for k in range(1,i):
                        b = b+np.array(data[k])
                    ax.bar(ind, data[i], 0.8, bottom=b,color = COLOR_MAP[i])
            # ax.set_title('Relative Contribution',fontsize=TITLE_SIZE, fontweight='bold')
            ax.spines['bottom'].set_color('#cac9c9')
            ax.spines['top'].set_color('#FFFFFF') 
            ax.spines['right'].set_color('#FFFFFF')
            ax.spines['left'].set_color('#cac9c9')
            ax.set_xlim([-1,len(ind)])
            ax.tick_params(axis='y',direction='in', color='#cac9c9', labelsize=LABEL_SIZE-4)
            ax.tick_params(axis='x',direction='in', length=0)
            ax.xaxis.set_visible(False)
            ax.set_yticks(np.arange(0, 1+0.1, 0.25))
            ax.legend(title='',labels=sig_list,loc='lower center',ncol=3, fontsize=LABEL_SIZE-4, edgecolor='white',
                      labelspacing=0.5, bbox_to_anchor=(0.5, (-0.1-(math.ceil(len(sig_list)/3)*0.065))))
            plt.savefig(pic+'SigContribution.pdf', dpi=300,bbox_inches='tight')
            print(colored(('=> Generate Bar Plot: ' + pic+'SigContribution.pdf'), 'green')) 
            
            height, length = len(sig_list), len(sample_list)  
            h_data = np.array(df.to_numpy())
            sns.set(font_scale=2)
            f,ax = plt.subplots(figsize=(9+length/20,2+height*0.3))
            ax = sns.heatmap(data, vmin=0, vmax=1, yticklabels = sig_list, linewidths=1,
                             square=False, cmap='Blues',cbar_kws={'orientation': 'horizontal','shrink':1, 'aspect':50})
            # ax.set_title('Signature Sample Heatmap', fontsize=TITLE_SIZE,weight='bold',va='bottom')
            ax.xaxis.set_visible(False)
            ax.set_xticklabels([])
            ax.tick_params(axis='both',length=0)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=LABEL_SIZE-4,color='#222222')
            plt.savefig(pic+'SigSamHeatmap.pdf',dpi=300,bbox_inches='tight')
            print(colored(('=> Generate Heatmap: '+pic+'SigSamHeatmap.pdf\n'), 'green'))
        def DonutPlot():
            df = pd.read_csv(output_folder+'sig_sample.csv', index_col=0)
            raw_data = df.sum(axis=1)/df.shape[1]
            SUM = raw_data.sum(axis=0)
            raw_data = raw_data/SUM
            names, sizes = list(raw_data.index), list(raw_data.iloc[:])
            names = [names[i]+': '+'{:.1%}'.format(sizes[i]) for i in range(len(sizes))]
            fig, ax = plt.subplots(figsize=(6, 8), subplot_kw=dict(aspect='equal'))
            wedges, texts = ax.pie(sizes, colors=COLOR_MAP[:len(names)],wedgeprops=dict(width=0.6,edgecolor='w',linewidth=2), startangle=-40) #,normalize=False

            bbox_props = dict(boxstyle='square,pad=0.3', fc='w', ec='k', lw=0)
            kw = dict(arrowprops=dict(arrowstyle='-'),bbox=bbox_props, zorder=0, va='center')

            for i, p in enumerate(wedges):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = np.sin(np.deg2rad(ang))
                x = np.cos(np.deg2rad(ang))
                horizontalalignment = {-1: 'right', 1: 'left'}[int(np.sign(x))]
                connectionstyle = 'angle,angleA=0,angleB={}'.format(ang)
                kw['arrowprops'].update({'connectionstyle': connectionstyle})
                ax.annotate(names[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw, fontsize=LABEL_SIZE)
            plt.savefig(pic+'Donut_plot.pdf', dpi=300, bbox_inches='tight')
            print(colored(('=> Generate Donut Plot: '+pic+'Donut_plot.pdf'), 'green'))
        nmf()
        SBSPlot()
        DonutPlot()
        CosineSimilarity()
        SigDistribution()



    def lsqnonneg(self, input):
        print(input)
        mut_matrix = pd.read_csv(input, sep = '\t', index_col = 0)
        signatures = pd.read_csv('lib/auxiliary/COSMIC_72.tsv', sep = '\t', index_col = 0)

        for i in range(mut_matrix.shape[1]):
            d, C = mut_matrix.iloc[:,i], signatures
            colName = list(C.columns)
            m, n = C.shape[0], C.shape[1]
            import sys
            from scipy import linalg
            
            tol = 10 * sys.float_info.epsilon * linalg.norm(C, ord=2) * (max(n, m)+1)
            
            x = [0]*n
            P, Z = [False]*n, [True]*n

            resid = d - np.matmul(C, x)
            w = np.matmul(C.T, resid)
            wz = [0]*n

            # iteration params
            outeriter, it = 0, 0
            itmax, exitFlag = 3*n, 1

            while pd.Series(Z).any() and (w[Z]>tol).any():
                outeriter += 1
                z = [0]*n
                wz = [-np.inf]*n
                wz = w
                im = colName.idx(wz.idxmax())
                print(im)
                
                os._exit(0)
            
            

            
            
            

            
            
        


        