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
bam_dedup_script = homeDir+'bam_dedup_4ends.py'
bam_dedup_2ends_script = homeDir+'bam_dedup.py'
make_fragments_script = homeDir+'make_fragments.py'
CB_counts_script = homeDir+'CB_counts.py'
cellQC_stats_script = homeDir+'cellQC_stats.R'
files_required.extend([bam_dedup_script,bam_dedup_2ends_script,make_fragments_script,CB_counts_script,cellQC_stats_script])

"""Data files"""
blacklist_file = homeDir+'downloadedInputs/hg38-blacklist.v2.bed'
otherChr_file = homeDir+'downloadedInputs/hg38_other_chromosome.txt'
promoters_file = homeDir+'downloadedInputs/hg38_GENCODEv28_promoters.bed'
files_required.extend([blacklist_file,otherChr_file])

iter_exists(files_required)


"""GOAL:
    Get files ready for cellQC_filter.R
    other chromosome deduplicated fragments file
    # non-deduplicated reads in blacklist per cell - deduplication removes blacklisted reads
    # deduplicated reads in promoters per cell
    # deduplicated reads in promoters per cell
    # non-deduplicated reads per cell - if not already done
"""

arg_obj = ArgumentParser(description='Counting reads from BAM per QC criteria for use in cellQC_filter.R')
arg_obj.add_argument('oriBam', metavar = 'oriBam', help='original bam file')
arg_obj.add_argument('bam', metavar = 'bam', help='deduplicated bam file')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix for output')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('peaks', metavar = 'peaks', help='peak file; ideally consensus peaks')
arg_obj.add_argument('-blacklist', metavar = 'blacklist', default=blacklist_file, help="Blacklist BED file to exclude from resulting BAM file. Default = Boyle Lab hg38 v2")
arg_obj.add_argument('-otherChr', metavar = 'otherChr', default=otherChr_file, help='file of extra regions to put through filtering: chr:start-stop per line; order matters here! Default is hg38 NOT 1-22XY')
arg_obj.add_argument('-promoters', metavar = 'promoters', default=promoters_file, help="promoters BED file. Default 2kb upstream of GENCODE v28 HAVANA protein-coding transcripts")
arg_obj.add_argument('-normFrag', metavar = 'normFrag', help='chr1-22XY chromosome fragments; if included, will do baseline stat plots')
arg_obj.add_argument('-statsTable', metavar = 'statsTable', help='filename with a criteria vs (cutoff, maxGenome, maxCrit) table')
arg_obj.add_argument('-bam_dedup_2ends', action='store_true', help='if given, use bam_dedup.py script instead of bam_dedup_4ends.py')
args = arg_obj.parse_args()

if not isfile(args.oriBam) or not isfile(args.bam) or not isfile(args.peaks) or not isfile(args.blacklist) or not isfile(args.otherChr) or not isfile(args.promoters):
    raise Exception('Check inputs.')

if args.normFrag is not None:
    if not isfile(args.normFrag):
        raise Exception("chr1-22XY chromosome fragment file doesn't exist.")
    
if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)

"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))


"""Other Chromosomes; chrM in fragments file."""
print("Starting other chromosomes")
print(time.strftime("%c"))
otherChr_bam = args.outDir+args.prefix+'_otherChr.bam'
if args.bam_dedup_2ends:
    cmd = "python %s %s %s %s"%(bam_dedup_2ends_script,args.oriBam,args.otherChr,otherChr_bam.rsplit('.',1)[0]+'.sam')
else:
    cmd = "python %s %s %s %s"%(bam_dedup_script,args.oriBam,args.otherChr,otherChr_bam)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(otherChr_bam):
    raise Exception("bam_dedup script failed.\nERROR: %s"%(ret[1]))
print(time.strftime("%c"))


otherChr_fragments = args.outDir+args.prefix+'_otherChr_fragments.tsv'
cmd = "python %s %s %s"%(make_fragments_script,otherChr_bam,otherChr_fragments)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(otherChr_fragments):
    raise Exception("make_fragments script failed.\nERROR: %s"%(ret[1]))
print(time.strftime("%c"))


"""blacklist"""
print("Running blacklist")
getCBct(args.oriBam,args.blacklist,args.outDir+args.prefix+'_blacklist.bam',args.outDir+args.prefix+'_blacklist_byCell.txt')


"""Peaks"""
print("Running Peaks")
getCBct(args.bam,args.peaks,args.outDir+args.prefix+'_peaks.bam',args.outDir+args.prefix+'_peaks_byCell.txt')


"""Promoters"""
print("Running promoters")
getCBct(args.bam,args.promoters,args.outDir+args.prefix+'_promoters.bam',args.outDir+args.prefix+'_promoters_byCell.txt')


"""CB_counts"""
WG_CB_counts_file = args.outDir+args.prefix+"_CB_counts.txt"
if not isfile(WG_CB_counts_file):
    print("Generating NOT deduped whole genome cell counts file")
    print(time.strftime("%c"))
    cmd = "python %s %s %s"%(CB_counts_script,args.oriBam,WG_CB_counts_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(WG_CB_counts_file):
        raise Exception("Count CB failed.\nERROR: %s"%(ret[1]))
    print(time.strftime("%c"))


if args.normFrag is not None:
    """Making Stat Plots"""
    print("Making Stat Plots")
    plotDir = args.outDir+"plots_minC-0/"
    mkdir(plotDir)
    print(time.strftime("%c"))
    cmd = "/usr/bin/time Rscript %s %s %s %s %s %s %s %s %s %s %s %s"%(cellQC_stats_script,args.normFrag,
                                                                       otherChr_fragments,WG_CB_counts_file,
                                                                       args.outDir+args.prefix+'_blacklist_byCell.txt',
                                                                       args.outDir+args.prefix+'_peaks_byCell.txt',
                                                                       args.outDir+args.prefix+'_promoters_byCell.txt',
                                                                       args.outDir+args.prefix+'_WG_noDedup_readCtByCell.RDS',
                                                                       args.outDir+args.prefix+'_WG_dedup_readCtByCell.RDS',
                                                                       0,plotDir,args.prefix)
    if args.statsTable is not None and isfile(args.statsTable):
        cmd+='--stats_table_file %s'%(args.statsTable)
    ret = getOneShellOutput(cmd)
    if "Execution halted" in ret[1]:
        print("cellQC likely failed.\n%s"%(ret[1]))
    print(time.strftime("%c"))


