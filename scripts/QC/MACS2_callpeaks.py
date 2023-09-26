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
control_BG_file = homeDir+'downloadedInputs/GSE92674_control_hg38_sorted.bedGraph'
files_required.extend([control_PE_file,control_BG_file])

iter_exists(files_required)

"""This version was used for donor-specific peak calling, which we only used to calculate overlap with the consensus peaks"""

arg_obj = ArgumentParser(description='call MACS2 peaks starting from BAM/SAM file format')
arg_obj.add_argument('input', metavar = 'input', help='BAM/SAM file for MACS2')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
arg_obj.add_argument('-input_style', metavar = 'input_style', default="PE", choices=["PE","BG"], help='Input file format for MACS2; options: PE - MACS2 BED Paired End; BG - BedGraph; Default = PE')
arg_obj.add_argument('-overlap_style', metavar = 'overlap_style', default="iter", choices=["iter","BTMT"], help='Overlapping peak method; options: BTMT - Best subpeak-trim-merge-trim; iter - choose 1 peak per overlapping region; Default = iter')
arg_obj.add_argument('-trim', metavar = 'trim', type=int, default = 100, help="Half the size of peaks; summit+/-trim; Default = 100 bp")
args = arg_obj.parse_args()

input_style = args.input_style
print("Using input style: %s"%(input_style))

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


"""convert to BAM file if not already"""
if input_file.rsplit('.',1)[1]!='bam':
    input_BAM = input_file.rsplit('.',1)[0]+'.bam'
    cmd = "samtools view -b %s > %s"%(input_file,input_BAM)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(input_BAM):
        raise Exception("ERROR: converting SAM to BAM.\n%s"%(ret[1]))

    input_file = input_BAM

if input_style=="BG":
    control_file = control_BG_file
    print("Control file: %s"%(control_file))
    if not isfile(control_file):
        raise Exception("Control file for MACS2 must exist: %s"%(control_file))
    
    print(time.strftime("%c"))
    """Convert BAM->bedgraph."""
    bedgraph_input_file = input_file.rsplit(".",1)[0]+'.bedgraph'
    cmd = "bedtools genomecov -ibam %s -bg > %s"%(input_file,bedgraph_input_file)
    ret = getOneShellOutput(cmd)
    if not isfile(bedgraph_input_file):
        raise Exception("ERROR: getting bedgraph file.\n%s"%(ret[1]))
    print(time.strftime("%c"))
    
    
    """bedgraph stats"""
    cmd = "cut -f1 %s | uniq"%(bedgraph_input_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or len(ret[0])==0:
        print("ERROR: getting bedgraph stats.\n%s"%(ret[1]))
    else:
        print("\nBEDGRAPH STATS:\n"+ret[0]+"\n")
    print(time.strftime("%c"))
elif input_style=="PE":
    control_file = control_PE_file
    print("Control file: %s"%(control_file))
    if not isfile(control_file):
        raise Exception("Control file for MACS2 must exist: %s"%(control_file))    
    
    
    print(time.strftime("%c"))
    """convert to MACS2 BEDPE file format"""
    input_BEDPE = input_file.rsplit('.',1)[0]+'_macsPE.bed'
    cmd="/usr/bin/time macs2 randsample -i %s -f BAMPE -p 100 -o %s"%(input_file,input_BEDPE)
    ret = getOneShellOutput(cmd)
    print("\macs2 BEDPE error stream (i.e. output):\n"+ret[1]+"\n")
    if not isfile(input_BEDPE):
        raise Exception("ERROR: getting MACS2 input BAMPE file.")
else:
    raise Exception("Input style needs to be set!")



"""MACS2 callpeaks"""
peak_prefix = basename(input_file).rsplit('.',1)[0]+"_"+input_style
narrowPeak_file = args.outDir+peak_prefix+"_peaks.narrowPeak"
summits_file = args.outDir+peak_prefix+"_summits.bed"

print(time.strftime("%c"))
if input_style=="BG":
    cmd = "macs2 callpeak -t %s -c %s --outdir %s -n %s --call-summits --nomodel --extsize %s -f BED --keep-dup all"%(bedgraph_input_file,control_file,args.outDir,peak_prefix,whole_trim)
elif input_style=="PE":
    cmd = "macs2 callpeak -t %s -c %s --outdir %s -n %s --call-summits -f BEDPE --keep-dup all"%(input_BEDPE,control_file,args.outDir,peak_prefix)
else:
    raise Exception("Input style needs to be set!")
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



