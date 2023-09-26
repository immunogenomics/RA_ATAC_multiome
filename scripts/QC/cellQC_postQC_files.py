import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import time

def getOneShellOutput(cmd):
    print(cmd)
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]

def GenerateTempFilename(temp_dir = '', base_string = None):
    import time
    import hashlib
    from random import random
    from os.path import isdir, isfile
    from subprocess import call
    
    if temp_dir!="" and not isdir(temp_dir):
        temp_dir=""
    
    if temp_dir!='' and temp_dir[-1]!="/":
        temp_dir=temp_dir+'/'
    
    fileTmp= temp_dir + hashlib.md5((time.ctime() + str(random())).encode('utf-8')).hexdigest()
    while isfile(fileTmp):
        fileTmp= temp_dir + hashlib.md5((time.ctime() + str(random())).encode('utf-8')).hexdigest()
    
    call('touch %s'%fileTmp,shell=True)
    
    return fileTmp

def mkdir(x):
    if not isdir(x):
        cmd = "mkdir %s"%(x)
        ret = getOneShellOutput(cmd)
        print(ret)
        
        if not isdir(x):
            return False
        else:
            return True
    else:
        return True


def getCBct(aFile,bFile,interFile,outFile):
    print(time.strftime("%c"))
    cmd = "intersectBed -nonamecheck -abam %s -b %s > %s"%(aFile,bFile,interFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(interFile):
        raise Exception("intersectBed failed.\nERROR: %s"%(ret[1]))
    print(time.strftime("%c"))
    
    cmd = "python %s %s %s"%(CB_counts_script,interFile,outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(outFile):
        raise Exception("Count CB failed.\nERROR: %s"%(ret[1]))
    print(time.strftime("%c"))


def parseCells(critFile,colNum,cellFile,outFile):
    print(time.strftime("%c"))
    if not isfile(critFile):
        print("SKIPPING: file doesn't exist: %s"%(critFile))
        return
    
    cmd = 'awk \'BEGIN{OFS=FS="%s"} FNR==NR{if(NR!=1) a[$1]=$1;next} {if($%s in a) print $0}\' %s %s > %s'%(r"\t",colNum,cellFile,critFile,outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(outFile):
        raise Exception("intersectBed failed.\nERROR: %s"%(ret[1]))
    print(time.strftime("%c"))
    

def exists(f):
    if not isfile(f):
        raise Exception("File doesn't exist: %s"%(f))

def iter_exists(l):
    for x in l:
        exists(x)



"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"

files_required = []
"""Scripts"""
CB_counts_script = homeDir+'CB_counts.py'
cellQC_stats_script = homeDir+'cellQC_stats.R'
files_required.extend([CB_counts_script,cellQC_stats_script])

iter_exists(files_required)

"""NOTE: this script assumes you have already done cellQC_files.py with outputs in QCdir using the same prefix and cellQC_filter.R to get the QCcells input file."""

arg_obj = ArgumentParser(description='subset files by post cell QC cells')
arg_obj.add_argument('bam', metavar = 'bam', help='deduplicated bam file')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix for output')
arg_obj.add_argument('QCdir', metavar = 'QCdir', help='cellQC input directory')
arg_obj.add_argument('postDir', metavar = 'postDir', help='postQC output directory')
arg_obj.add_argument('peaks', metavar = 'peaks', help='peak file; ideally consensus peaks')
arg_obj.add_argument('fragments', metavar = 'fragments', help='deduplicated fragments file')
arg_obj.add_argument('QCcells', metavar = 'QCcells', help='QCcells file; output in cellQC_filter.R')
arg_obj.add_argument('-statsTable', metavar = 'statsTable', help='filename with a criteria vs (cutoff, maxGenome, maxCrit) table')
arg_obj.add_argument('-samePeaks', action = 'store_true', help=argparse.SUPPRESS) #if same peaks used for cellQC_files.py, this script with subset that output. if not, it will regenerate the counts with the new peaks.
args = arg_obj.parse_args()

if not isfile(args.bam) or not isdir(args.QCdir) or not isfile(args.peaks) or not isfile(args.fragments) or not isfile(args.QCcells):
    raise Exception('Check inputs.')

if args.QCdir[-1]!='/':
    args.QCdir+='/'
    
if args.postDir[-1]!='/':
    args.postDir+='/'
mkdir(args.postDir)

"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))


outPrefix = args.postDir+args.prefix+"_cellQC"

"""deduped fragments"""
print("PARSING normal chromosomes")
PARSE_normal_fragments = outPrefix+'_fragments.tsv'
parseCells(args.fragments,4,args.QCcells,PARSE_normal_fragments)


"""Other Chromosomes; chrM in fragments file."""
print("PARSING other chromosomes")
otherChr_fragments = args.QCdir+args.prefix+'_otherChr_fragments.tsv'
PARSE_otherChr_fragments = outPrefix+'_otherChr_fragments.tsv'
parseCells(otherChr_fragments,4,args.QCcells,PARSE_otherChr_fragments)


"""blacklist"""
print("PARSING blacklist")
blacklist_CB_ct = args.QCdir+args.prefix+'_blacklist_byCell.txt'
PARSE_blacklist_CB_ct = outPrefix+'_blacklist_byCell.txt'
parseCells(blacklist_CB_ct,1,args.QCcells,PARSE_blacklist_CB_ct)


"""Peaks"""
if args.samePeaks:
    print("PARSING peaks")
    peaks_CB_ct = args.QCdir+args.prefix+'_peaks_byCell.txt'
    PARSE_peaks_CB_ct = outPrefix+'_peaks_byCell.txt'
    parseCells(peaks_CB_ct,1,args.QCcells,PARSE_peaks_CB_ct)
else:
    print("Running NEW Peaks")
    PARSE_peaks_CB_ct = outPrefix+'_peaks_byCell.txt'
    getCBct(args.bam,args.peaks,outPrefix+'_peaks.bam',outPrefix+'_peaks_byCell.txt')


"""Promoters"""
print("PARSING promoters")
promoters_CB_ct = args.QCdir+args.prefix+'_promoters_byCell.txt'
PARSE_promoters_CB_ct = outPrefix+'_promoters_byCell.txt'
parseCells(promoters_CB_ct,1,args.QCcells,PARSE_promoters_CB_ct)


"""CB_counts"""
print("PARSING NOT deduped whole genome cell counts file")
WG_CB_ct = args.QCdir+args.prefix+'_CB_counts.txt'
PARSE_WG_CB_ct = outPrefix+'_CB_counts.txt'
parseCells(WG_CB_ct,1,args.QCcells,PARSE_WG_CB_ct)


"""Making Stat Plots"""
print("Making Stat Plots")
plotDir = args.postDir+"plots_minC-%s/"%(0)
mkdir(plotDir)
print(time.strftime("%c"))
cmd = "/usr/bin/time Rscript %s %s %s %s %s %s %s %s %s %s %s %s"%(cellQC_stats_script,PARSE_normal_fragments,
                                                                   PARSE_otherChr_fragments,PARSE_WG_CB_ct,
                                                                   PARSE_blacklist_CB_ct,PARSE_peaks_CB_ct,PARSE_promoters_CB_ct,
                                                                   outPrefix+'_WG_noDedup_readCtByCell.RDS',
                                                                   outPrefix+'_WG_dedup_readCtByCell.RDS',
                                                                   0,plotDir,args.prefix+"_cellQC")
if args.statsTable is not None and isfile(args.statsTable):
    cmd+='--stats_table_file %s'%(args.statsTable)
ret = getOneShellOutput(cmd)
if "Execution halted" in ret[1]:
    print("cellQC likely failed.\n%s"%(ret[1]))
print(time.strftime("%c"))


