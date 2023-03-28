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


arg_obj = ArgumentParser()
arg_obj.add_argument('frag_listing_file', metavar='frag_listing_file', help='fragment listing file')
arg_obj.add_argument('cell_listing_file', metavar='cell_listing_file', help='cell listing file')
arg_obj.add_argument('-header',action='store_true', help='If given, skip the header in the good cells file')
arg_obj.add_argument('outDir', metavar='outDir', help='Output directory')
arg_obj.add_argument('prefix', metavar='prefix', help='Prefix for output files; will add donor"_fragments.tsv" to end!')
arg_obj.add_argument('-lsf', metavar = 'lsf', help='comma-delimited list: queue,time (min),mem (MB)')
arg_obj.add_argument('-tabix', action='store_true', help='if given, also make tabix file.')
args = arg_obj.parse_args()


if not isfile(args.frag_listing_file) or not isfile(args.cell_listing_file):
    raise Exception("Input file(s) must exist.")


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

if useDefaults:
    queue='medium' 
    Wtime='300'
    mem='2000'


#read in all fragment files and verify they exist
with open(args.frag_listing_file) as f:
    frag_files = f.readlines()
frag_files = [f.strip() for f in frag_files]
if not all([isfile(f) for f in frag_files]):
    raise Exception("Not all fragment files exist.")


#read in all cell files and verify they exist
with open(args.cell_listing_file) as f:
    cell_files = f.readlines()
cell_files = [f.strip() for f in cell_files]
if not all([isfile(f) for f in cell_files]):
    raise Exception("Not all cell files exist.")


#match up IDs
submitted_jobs = []
for frag_file in frag_files:
    
    #get ID from fragment file
    thisID = set(re.findall('(BRI-[0-9]+)',frag_file))
    if len(thisID)!=1:
        raise Exception("more than 1 ID in fragment file: %s"%(frag_file))
    thisID = thisID.pop()
    
    #find cell files with this ID
    theseFiles = [f for f in cell_files if thisID in f]
    
    for cell_file in theseFiles:
        
        this_outFile = args.outDir+args.prefix+'_'+thisID+'_fragments.tsv'
        
        internal_cmd = 'awk \'BEGIN{OFS=FS="%s"} FNR==NR{%sa[$1]=$1;next} {if($4 in a) print $0}\' %s %s | sort -k1,1 -k2,2n > %s'%(r"\t","if(NR!=1) " if args.header else "", cell_file,frag_file,this_outFile)
        
        if args.tabix:
            bgzip_outFile = this_outFile+'.gz'
            internal_cmd = internal_cmd+"\n\n"+"bgzip %s"%(this_outFile)
            
            tabix_outFile = this_outFile+'.gz.tbi'
            internal_cmd = internal_cmd+"\n\n"+"tabix -p bed %s"%(bgzip_outFile)
            
        this_script = clusterLogs+args.prefix+'_'+thisID+'_fragments.sh'
        
        with open(this_script,'w') as f:
            f.write(internal_cmd+'\n')
        
        cmd = '/usr/bin/time sh %s'%(this_script)

        job_name='subFrag_%s'%(args.prefix+'_'+thisID)
        queue_call = 'bsub -R "select[hname!=cn007]" -R "select[hname!=cn001]" -R "select[hname!=cn002]" -R "select[hname!=cn003]" -R "select[hname!=cn004]" -R "select[hname!=cn005]" -q %s -W %s -M %s -J %s -o %s-%%J.out -e %s-%%J.out -cwd %s "%s"'%(queue,Wtime,mem,job_name,job_name,job_name,clusterLogs,cmd)
        #print(queue_call)
        job_id = SubmitJobCompact(queue_call) #prints job_string within
        submitted_jobs.append(job_id)
        #print('\n')
    

    
