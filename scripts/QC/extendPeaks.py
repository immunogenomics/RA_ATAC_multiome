import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import re
import decimal
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
plusMinus_bed_script = homeDir+'plusMinus_bed.py'
files_required.extend([plusMinus_bed_script])

"""Data files"""
hg38_chrom_sizes_file = homeDir+'downloadedInputs/hg38.chrom.sizes'
files_required.extend([hg38_chrom_sizes_file])

iter_exists(files_required)


arg_obj = ArgumentParser(description='extend peaks into peak neighborhoods')
arg_obj.add_argument('input_peak_file', metavar = 'input_peak_file', help='peak input file (narrowPeak)')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('peak_prefix', metavar = 'peak_prefix', help='prefix for peak files - do not need the _peaks.extension')
arg_obj.add_argument('-extend', metavar = 'extend', default = '5,Y,both,A,default', help='extend peaks options: flank,percentageYN,type,bodyYNA,chromSizesFile;\ndefault: 5,Y,both,A,default')
args = arg_obj.parse_args()

if not isfile(args.input_peak_file):
    raise Exception("input peak file must exist.")

if args.outDir[-1]!='/':
    args.outDir+='/'

#validating peak extension argument
if args.extend.count(',')!=4:
    raise Exception("Not the right amount of commas in extend argument.")

flank, percentageYN, typ, bodyYNA, chromSizesFile = args.extend.split(',')
if percentageYN not in ['Y','N']:
    raise Exception("Extend percentage must be Y (yes) or N (no).")
if bodyYNA not in ['Y','N','A']:
    raise Exception("Extend including gene body must be Y (yes) or N (no) or A (not applicable).")
if typ not in ['start','end','both','center','summit']:
    raise Exception("Extend type must be either start, end, both, center (or summit if given in column 10)")
if bodyYNA=="Y" and typ not in ['start','end']:
    raise Exception("Extend including gene body must have a type of start or end")
if chromSizesFile!="default" and not isfile(chromSizesFile):
    raise Exception("Extend: chromosome sizes file must be default or a valid filename")

if chromSizesFile=="default":
    chromSizesFile = hg38_chrom_sizes_file

if typ=='both':
    phrase = 'both sides'
elif typ in ['start','end']:
    phrase = typ+' side'
if bodyYNA=="Y":
    phrase += ' (including gene bodies)'
elif typ in ['center','summit']:
    phrase = typ

print("Extending peaks by %s %s from %s using chromosize sizes in %s"%(flank,"times" if percentageYN=="Y" else "bp",phrase, chromSizesFile))


"""removing subpeaks"""
noSub_file = args.outDir+args.peak_prefix+"_noSub_peaks.narrowPeak"
cmd = "grep -P 'peak_[0-9]+a?%s' %s > %s"%(r"\t",args.input_peak_file,noSub_file)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(noSub_file):
    raise Exception("ERROR: removing subpeaks.\n%s"%(ret[1]))


"""Extending"""
extend_file = args.outDir+args.peak_prefix+"_extended_peaks.narrowPeak"
extraArgs = ''
if percentageYN=="Y":
    extraArgs+=' -perc'
if bodyYNA=="Y":
    extraArgs+=' -keepBody'
extraArgs+=' '
cmd = "python %s%s%s %s %s %s %s"%(plusMinus_bed_script,extraArgs,noSub_file,extend_file,flank,typ,chromSizesFile)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(extend_file):
    raise Exception("ERROR: extending narrowPeaks.\n%s"%(ret[1]))


"""Merging extended peaks"""
merge_extend_file = args.outDir+args.peak_prefix+"_extended_peaks.bed"
cmd = "cut -f1-3 %s | sort -k1,1 -k2,2n | mergeBed -i - > %s"%(extend_file,merge_extend_file)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(merge_extend_file):
    raise Exception("ERROR: merging extended narrowPeaks.\n%s"%(ret[1]))

    
