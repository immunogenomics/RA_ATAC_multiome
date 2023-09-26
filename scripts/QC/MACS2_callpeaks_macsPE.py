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
maxSubPeak_script = homeDir+'maxSubPeak.py'
peak_iterative_overlap_script = homeDir+'peak_iterative_overlap.py'
files_required.extend([maxSubPeak_script,peak_iterative_overlap_script])

"""Data files"""
control_PE_file = homeDir+'downloadedInputs/GSE92674_control_hg38_macsPE_sorted.bed'
files_required.extend([control_PE_file])

iter_exists(files_required)

"""This version was used for consensus peak calling after merging all donors' BEDPE files"""

arg_obj = ArgumentParser(description='call MACS2 peaks starting from BED PE file format')
arg_obj.add_argument('input', metavar = 'input', help='BED PE file for MACS2')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('-addToPrefix', metavar = 'addToPrefix', help='suffix to add to peak prefix')
arg_obj.add_argument('-peakPrefix', metavar = 'peakPrefix', help=argparse.SUPPRESS)
arg_obj.add_argument('-overlap_style', metavar = 'overlap_style', default="iter", choices=["iter","BTMT"], help='Overlapping peak method; options: BTMT - Best subpeak-trim-merge-trim; iter - choose 1 peak per overlapping region; Default = iter')
arg_obj.add_argument('-trim', metavar = 'trim', type=int, default = 100, help="Half the size of peaks; summit+/-trim; Default = 100 bp")
args = arg_obj.parse_args()

print("Using input style: %s"%('macsPE'))

overlap_style = args.overlap_style
print("Using peak overlap style: %s"%(overlap_style))

if not isfile(args.input):
    raise Exception("input file must exist.")

if args.outDir[-1]!='/':
    args.outDir+='/'

input_file = args.input
print("Input file: %s"%(input_file))

half_trim = args.trim
whole_trim = half_trim*2
print("Trimming peaks to %s bp"%(whole_trim))

control_file = control_PE_file
print("Control file: %s"%(control_file))
if not isfile(control_file):
    raise Exception("Control file for MACS2 must exist: %s"%(control_file))    



"""MACS2 callpeaks"""
peak_prefix = basename(input_file).rsplit('.',1)[0]
if args.addToPrefix is not None:
        if peak_prefix.rsplit('_',1)[-1]=='macsPE':
            peak_prefix = peak_prefix.rsplit('_',1)[0]+'_'+args.addToPrefix+'_macsPE'
        else:
            peak_prefix += '_'+args.addToPrefix
if args.peakPrefix is not None:
    peak_prefix = args.peakPrefix
narrowPeak_file = args.outDir+peak_prefix+"_peaks.narrowPeak"
summits_file = args.outDir+peak_prefix+"_summits.bed"

print(time.strftime("%c"))
cmd = "macs2 callpeak -t %s -c %s --outdir %s -n %s --call-summits -f BEDPE --keep-dup all"%(input_file,control_file,args.outDir,peak_prefix)
ret = getOneShellOutput(cmd)
print("\macs2 callpeak error stream (i.e. output):\n"+ret[1]+"\n")
if not isfile(narrowPeak_file) or not isfile(summits_file):
    raise Exception("ERROR: getting MACS2 peaks.")
print(time.strftime("%c"))


"""MACS2 choose/merge subpeaks + trimming"""
print(time.strftime("%c"))
if overlap_style=="BTMT":
    
    """Choose subpeaks by qvalue"""
    best_peak_file = narrowPeak_file.rsplit('.',1)[0]+'_bestSubPeak.'+narrowPeak_file.rsplit('.',1)[1]
    cmd = "python %s %s -peak"%(maxSubPeak_script,narrowPeak_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(best_peak_file):
        raise Exception("ERROR: getting best subpeak file.\n%s"%(ret[1]))
    
    """trim, merge, trim"""
    toUse_peaks_file = args.outDir+peak_prefix+"_peaks_toUse.bed"
    cmd = 'awk \'BEGIN{OFS=FS="%s"} {peak=$2+$10;print $1,peak-%s,peak+%s}\' %s | sort -k1,1 -k2,2n | mergeBed -i - | awk \'BEGIN{OFS=FS="%s"} {mid=int(($3+$2)/2);print $1,mid-%s,mid+%s}\' > %s'%(r"\t",half_trim,half_trim,best_peak_file,r"\t",half_trim,half_trim,toUse_peaks_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(toUse_peaks_file):
        raise Exception("ERROR: trim, merge, trim peaks.\n%s"%(ret[1]))
elif overlap_style=="iter":
    
    """trim & sort"""
    trim_sorted_peaks_file = narrowPeak_file.rsplit('.',1)[0]+'_trim%s_sort.'%(whole_trim)+narrowPeak_file.rsplit('.',1)[1]
    cmd = 'awk \'BEGIN{OFS=FS="%s";half=%s} {peak=$2+$10;$2=peak-half;$3=peak+half;$10=half;print $0}\' %s | sort -k1,1 -k2,2n > %s'%(r"\t",half_trim,narrowPeak_file,trim_sorted_peaks_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(trim_sorted_peaks_file):
        raise Exception("ERROR: getting trimmed and sorted peaks.\n%s"%(ret[1]))
    
    """iterative over peaks, choosing best 1 out of any overlap"""
    toUse_peaks_file = args.outDir+peak_prefix+"_peaks_toUse.bed"
    cmd = "python %s %s -outFile %s"%(peak_iterative_overlap_script,trim_sorted_peaks_file,toUse_peaks_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(toUse_peaks_file):
        raise Exception("ERROR: getting peak iterative overlap file.\n%s"%(ret[1]))
else:
    raise Exception("Peak Overlap style needs to be set!")
print(time.strftime("%c"))



