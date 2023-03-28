import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
from os import listdir
import re
import decimal
from submitERISone import SubmitJobCompact, SubmitJobsAndWait, WhichJobsAreRunning, WaitForJobs
import time

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

def print_time(msg):
    print("\n%s - %s\n"%(time.strftime("%c"),msg))

def exists(f):
    if not isfile(f):
        raise Exception("File doesn't exist: %s"%(f))

def iter_exists(l):
    for x in l:
        exists(x)


"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"

"""Scripts"""
idvf_lm_script = homeDir+'identity_vs_functional_linear_model.R'
exists(idvf_lm_script)


arg_obj = ArgumentParser()
arg_obj.add_argument('data_file', metavar='data_file', help='Matrix file: binary peaks x cells')
arg_obj.add_argument('meta_file', metavar='meta_file', help='Cell metadata file')
arg_obj.add_argument('-sample_col', metavar='sample_col', help='sample/donor column in meta file; default in Rscript: sample')
arg_obj.add_argument('-frag_col', metavar='frag_col', help='number of fragments column in meta file; default in Rscript: nFrags')
arg_obj.add_argument('-cellType_col', metavar = 'cellType_col', help='cell type column; default in Rscript: cluster_name')
arg_obj.add_argument('-simplifyCT', action='store_true', help='if given, change X-#: description to X#') 
arg_obj.add_argument('lineage_file', metavar='lineage_file', help='cells x lineage file: 1 column: 1 if CD4+ CD8A-; 0 if double pos/neg; -1 if CD4- CD8A+')
arg_obj.add_argument('outDir', metavar='outDir', help='Output directory')
arg_obj.add_argument('prefix', metavar='prefix', help='Prefix for output files')
arg_obj.add_argument('peak2gene_file', metavar='peak2gene_file', help='peak to gene file')
arg_obj.add_argument('subsetPeak_file', metavar='subsetPeak_file', help='subsetPeak_file')
arg_obj.add_argument('numToRun', metavar='numToRun', type=int, help='number of peaks to run in one go')
arg_obj.add_argument('-lsf', metavar = 'lsf', help='comma-delimited list: queue,time (min),mem (MB)')
args = arg_obj.parse_args()


if not isfile(args.data_file) or not isfile(args.meta_file) or not isfile(args.lineage_file) or not isfile(args.peak2gene_file) or not isfile(args.subsetPeak_file):
    raise Exception("Input files must exist.")    


if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)
clusterLogs = args.outDir+'clusterLogs/'
mkdir(clusterLogs)


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')


useDefaults = args.lsf is None
if not useDefaults:
    avail_queues = {'vshort':[15,4000],'short':[60,4000],'medium':[1440,8000],'normal':[4320,8000],'long':[10080,8000],'vlong':[40320,8000],'big':[40320,498000],'big-multi':[40320,498000]}
    queue, Wtime, mem = args.lsf.split(',')
    if queue not in avail_queues.keys():
        print("Queue doesn't exist on ERISone; using defaults")
        useDefaults = True
    if int(Wtime)>avail_queues.get(queue)[0]:
        print("Time limit must be less than or equal to %s for queue %s; using defaults"%(Wtime,avail_queues.get(queue)[0]))
        useDefaults = True
    if int(mem)>avail_queues.get(queue)[1]:
        print("Memory limit must be less than or equal to %s for queue %s; using defaults"%(mem,avail_queues.get(queue)[1]))
        useDefaults = True
    swap=None

if useDefaults:
    queue='medium'
    Wtime='900'
    mem='8000'


intermInDir = args.outDir+'peaksToRunFiles/'
mkdir(intermInDir)
intermOutDir = args.outDir+'intermediateOutputs/'
mkdir(intermOutDir)


input_files = []
input_prefix = intermInDir + args.prefix + "_peaksToRun_"
currEnd = args.numToRun
currOut_file = input_prefix+"0_%s.txt"%(currEnd-1)
input_files.append(currOut_file)
OUT = open(currOut_file,'w')
IN = open(args.subsetPeak_file,'r')
for i,line in enumerate(IN):
    if i<currEnd:
        OUT.write(line.strip()+'\n')
    elif i==currEnd:
        OUT.close()
        currEnd = currEnd+args.numToRun
        currOut_file = input_prefix+"%s_%s.txt"%(i,currEnd-1)
        input_files.append(currOut_file)
        OUT = open(currOut_file,'w')
        OUT.write(line.strip()+'\n')
    
OUT.close()
IN.close()

submitted_jobs = []
for input_file in input_files:
    currPrefix = basename(input_file).rsplit('.',1)[0]
    
    extraArgs = ''
    if args.sample_col is not None:
        extraArgs+=' --sample_col %s'%(args.sample_col)
    if args.frag_col is not None:
        extraArgs+=' --frag_col %s'%(args.frag_col)
    if args.cellType_col is not None:
        extraArgs+=' --cellType_col %s'%(args.cellType_col)
    if args.simplifyCT:
        extraArgs+=' --simplifyCT'
    extraArgs+=' '
    
    cmd = "/usr/bin/time Rscript %s%s%s %s %s %s %s --peak2gene_file %s --subsetPeak_file %s"%(idvf_lm_script, extraArgs, args.data_file, args.meta_file, args.lineage_file, intermOutDir, currPrefix, args.peak2gene_file, input_file)
    
    job_name='LLM_%s'%(currPrefix)
    queue_call = 'bsub -R "select[hname!=cn007]" -R "select[hname!=cn001]" -R "select[hname!=cn002]" -R "select[hname!=cn003]" -R "select[hname!=cn004]" -R "select[hname!=cn005]" -q %s -W %s -M %s -J %s -o %s-%%J.out -e %s-%%J.out -cwd %s "%s"'%(queue,Wtime,mem,job_name,job_name,job_name,clusterLogs,cmd)
    #print(queue_call)
    job_id = SubmitJobCompact(queue_call) #prints job_string within
    submitted_jobs.append(job_id)
    print('\n')
    

    
