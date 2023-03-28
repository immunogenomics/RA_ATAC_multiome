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

script_file = './bam_subset_merge_iter.py'
exists(script_file)


arg_obj = ArgumentParser()
arg_obj.add_argument('gene', metavar='gene', help='gene name')
arg_obj.add_argument('region', metavar='region', help='region to grep')
arg_obj.add_argument('state_listing_file', metavar='state_listing_file', help='cell state listing file: 2 columns: cell state in filename format and read count scaling coefficient')
arg_obj.add_argument('listing_dir', metavar='listing_dir', help='listing directory to find cell state bam file listing files')
arg_obj.add_argument('-byCS',action='store_true', help='If given, do individual cell state files')
arg_obj.add_argument('tmpDir', metavar='tmpDir', help='temporary directory')
arg_obj.add_argument('genome_file', metavar='genome_file', help='genome file')
arg_obj.add_argument('script_outDir', metavar='script_outDir', help='script output directory')
arg_obj.add_argument('script_prefix', metavar='script_prefix', help='script prefix')
arg_obj.add_argument('outDir', metavar='outDir', help='Output directory')
arg_obj.add_argument('prefix', metavar='prefix', help='Prefix for output files')
args = arg_obj.parse_args()


if not isfile(args.state_listing_file) or not isdir(args.listing_dir) or not isdir(args.tmpDir) or not isfile(args.genome_file):
    raise Exception("Input file(s) must exist.")

xx = re.match('chr[0-9XY]+:([0-9]+)-([0-9]+)',args.region)
if xx is None or len(xx.groups())!=2 or int(list(xx.groups())[1])-int(list(xx.groups())[0])<=0:
    raise Exception("Input region must be in form chr[0-9XY]+:[0-9]+-[0-9]+ and have the stop value greater than the start value")
    
if args.tmpDir[-1]!='/':
    args.tmpDir+='/'
    
if args.script_outDir[-1]!='/':
    args.script_outDir+='/'
mkdir(args.script_outDir)

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')


state_dic = {}
with open(args.state_listing_file) as f:
    lines = f.readlines()
for ll in lines:
    tabs = ll.strip().split('\t')
    if not tabs[1].replace('.','',1).isdigit():
        raise Exception('second number of state listing file must be a float')
    state_dic[tabs[0]] = float(tabs[1])


listing_files = listdir(args.listing_dir)

if not args.byCS:
    outFile = args.outDir+args.prefix+'.sh'
    OUT = open(outFile,'w')

for cs,coef in state_dic.items():
    
    if args.byCS:
        outFile = args.outDir+args.prefix+'_'+cs+'.sh'
        OUT = open(outFile,'w')
    
    thisListing = [xx for xx in listing_files if cs in xx]
    if len(thisListing)!=1:
        raise Exception("more than 1 listing file associated to cell state: %s"%(cs))
    thisListing = args.listing_dir+thisListing[0]
    
    thisBamPrefix = "%s%s_%s_%s"%(args.script_outDir,args.script_prefix,cs,args.gene)
    
    command1 = "/usr/bin/time python %s %s %s %s%s_%s/ %s.bam -genome_file %s &> %s%s_%s_%s.out"%(script_file,thisListing,args.region,args.tmpDir,args.gene,cs.split('_')[0].replace('-',''), thisBamPrefix,args.genome_file,args.outDir,args.script_prefix,cs,args.gene)
    OUT.write(command1+'\n\n')
    
    command2 = "bedtools genomecov -ibam %s.bam -bg > %s.bedgraph"%(thisBamPrefix,thisBamPrefix)
    OUT.write(command2+'\n\n')
    
    command3 = "awk 'BEGIN{OFS=FS=\"%s\"} {print $1,$2,$3,$4/%s}' %s.bedgraph > %s_scaleTo1e7.bedgraph"%(r"\t",coef,thisBamPrefix,thisBamPrefix)
    OUT.write(command3+'\n\n\n')
    
    if args.byCS:
        OUT.close()

if not args.byCS:
    OUT.close()
    
    
