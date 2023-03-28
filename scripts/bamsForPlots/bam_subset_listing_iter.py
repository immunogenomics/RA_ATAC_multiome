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
arg_obj.add_argument('bam_listing_file', metavar='bam_listing_file', help='bam listing file')
arg_obj.add_argument('cell_listing_file', metavar='cell_listing_file', help='cell listing file')
arg_obj.add_argument('-header',action='store_true', help='If given, skip the header in the good cells file')
arg_obj.add_argument('outDir', metavar='outDir', help='Output directory')
arg_obj.add_argument('-lsf', metavar = 'lsf', help='comma-delimited list: queue,time (min),mem (MB)')
args = arg_obj.parse_args()


if not isfile(args.bam_listing_file) or not isfile(args.cell_listing_file):
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
    Wtime='1200'
    mem='1000'


#read in all bam files and verify they exist
with open(args.bam_listing_file) as f:
    bam_files = f.readlines()
bam_files = [f.strip() for f in bam_files]
if not all([isfile(f) for f in bam_files]):
    raise Exception("Not all BAM files exist.")


#read in all cell files and verify they exist
with open(args.cell_listing_file) as f:
    cell_files = f.readlines()
cell_files = [f.strip() for f in cell_files]
if not all([isfile(f) for f in cell_files]):
    raise Exception("Not all cell files exist.")


#match up IDs
submitted_jobs = []
for bam_file in bam_files:
    
    #get ID from BAM file
    thisID = set(re.findall('(BRI-[0-9]+)',bam_file))
    if len(thisID)!=1:
        raise Exception("more than 1 ID in BAM file: %s"%(bam_file))
    thisID = thisID.pop()
    
    #find cell files with this ID
    theseFiles = [f for f in cell_files if thisID in f]
    
    for cell_file in theseFiles:
        
        currPrefix = basename(cell_file).rsplit('.',1)[0]
        
        this_outFile = args.outDir+currPrefix+'.bam'
        
        internal_cmd = 'samtools view -h %s | awk \'BEGIN{OFS=FS="%s"} FNR==NR{%sa[$1]=$1;next} {if(substr($0,1,1)=="@") {print $0} else {for (i=1;i<=NF;i++) {if(substr($i,1,4)=="CB:Z") currCB=substr($i,6,length($i)-5)};if(currCB in a) print $0}}\' %s - | samtools view -b | samtools sort > %s'%(bam_file,r"\t","if(NR!=1) " if args.header else "", cell_file,this_outFile)
        
        internal_cmd2 = 'samtools index %s'%(this_outFile)
        
        this_script = clusterLogs+currPrefix+'.sh'
        
        with open(this_script,'w') as f:
            f.write(internal_cmd+'\n')
            f.write(internal_cmd2+'\n')
        
        cmd = '/usr/bin/time sh %s'%(this_script)

        job_name='subBam_%s'%(currPrefix)
        queue_call = 'bsub -R "select[hname!=cn044]" -R "select[hname!=cn007]" -R "select[hname!=cn001]" -R "select[hname!=cn002]" -R "select[hname!=cn003]" -R "select[hname!=cn004]" -R "select[hname!=cn005]" -q %s -W %s -M %s -J %s -o %s-%%J.out -e %s-%%J.out -cwd %s "%s"'%(queue,Wtime,mem,job_name,job_name,job_name,clusterLogs,cmd)
        #print(queue_call)
        job_id = SubmitJobCompact(queue_call) #prints job_string within
        submitted_jobs.append(job_id)
        print('\n')
    

    
