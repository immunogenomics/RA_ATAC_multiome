import scrublet as scr
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import isdir, isfile
import argparse
from argparse import ArgumentParser

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

def getOneShellOutput(cmd):
    print(cmd)
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]

def mkdir(x):
    if isdir(x):
        return True
    
    cmd = "mkdir %s"%(x)
    ret = getOneShellOutput(cmd)
    
    return isdir(x)


arg_obj = ArgumentParser(description='find doublets in gene space')
arg_obj.add_argument('counts_inFile', metavar = 'counts_inFile', help='imput genes x cells matrix in MM format')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('prefix', metavar = 'prefix', help="prefix for output files")
arg_obj.add_argument('-doublet_rate', metavar='doublet_rate', type=float, default=0.15, help="expected doublet rate; default: 0.15") 
arg_obj.add_argument('-min_counts', metavar='min_counts', type=int, default=2, help="minimum counts; default: 2")
arg_obj.add_argument('-min_cells', metavar='min_cells', type=int, default=3, help="minumum cell counts; default: 3")
arg_obj.add_argument('-min_geneVar', metavar='min_geneVar',type=int, default=85, help="minimum gene variability pctl; default: 85")
arg_obj.add_argument('-nPCs', metavar='nPCs', type=int, default=30, help="number of PCs; default: 30")
arg_obj.add_argument('-nNeighbors', metavar='n_neighbors',type=int, default=10, help="number of neighbors for UMAP; default 10")
arg_obj.add_argument('-min_dist', metavar='min_dist', type=float, default=0.3, help="minimum distance for UMAP; default=0.3")
args = arg_obj.parse_args()

if not isfile(args.counts_inFile):
    raise Exception("Input file must exist.")

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))


counts_matrix = scipy.sparse.coo_matrix.transpose(scipy.io.mmread(args.counts_inFile))
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.doublet_rate)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts, 
                                                          min_cells=args.min_cells, 
                                                          min_gene_variability_pctl=args.min_geneVar, 
                                                          n_prin_comps=args.nPCs)

scrub.plot_histogram()
plt.savefig('%s%s_scrublet_histogram.png'%(args.outDir,args.prefix))

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, args.nNeighbors, min_dist=args.min_dist))

scrub.plot_embedding('UMAP', order_points=True)
plt.savefig('%s%s_scrublet_UMAP.png'%(args.outDir,args.prefix))

np.save('%s%s_scrubletRes.npy'%(args.outDir,args.prefix), doublet_scores)
