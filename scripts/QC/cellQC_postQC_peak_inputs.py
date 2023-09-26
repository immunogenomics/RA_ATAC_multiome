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


"""Input Arguments"""
arg_obj = ArgumentParser(description='subset deduplicated bam file by post cellQC cell barcodes & make a corresponding BEDPE file for MACS2 peak calling')
arg_obj.add_argument('bam_dedup',metavar = 'bam_dedup', help='deduplicated bam file')
arg_obj.add_argument('prefix',metavar = 'prefix', help='prefix')
arg_obj.add_argument('goodCells',metavar = 'goodCells', help='post cellQC file; cell barcodes in first column; with header')
arg_obj.add_argument('outDir', metavar = 'outDir', help='output directory')
args = arg_obj.parse_args()



"""Verifying Input Arguments"""
if not isfile(args.bam_dedup) or not isfile(args.goodCells):
    raise Exception("Inputs must exist.")
    
if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')


print(time.strftime("%c"))
bam_cellQC_file = '%s%s_cellQC.bam'%(args.outDir,args.prefix)
cmd = 'samtools view -h %s | awk \'BEGIN{OFS=FS="%s"} FNR==NR{if(NR!=1) a[$1]=$1;next} {if(substr($0,1,1)=="@") {print $0} else {for (i=1;i<=NF;i++) {if(substr($i,1,4)=="CB:Z") currCB=substr($i,6,length($i)-5)};if(currCB in a) print $0}}\' %s - | samtools view -b > %s'%(args.bam_dedup,r"\t",args.goodCells,bam_cellQC_file)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(bam_cellQC_file):
    raise Exception("ERROR: getting bam cell QC file\n%s"%(ret[1]))

cmd = 'samtools index %s'%(bam_cellQC_file)
ret = getOneShellOutput(cmd)
if not isfile(bam_cellQC_file+'.bai'):
    raise Exception("index command failed.\nERROR: %s"%(ret[1]))


print(time.strftime("%c"))
"""convert to MACS2 BEDPE file format"""
input_BEDPE = bam_cellQC_file.rsplit('.',1)[0]+'_macsPE.bed'
cmd="/usr/bin/time macs2 randsample -i %s -f BAMPE -p 100 -o %s"%(bam_cellQC_file,input_BEDPE)
ret = getOneShellOutput(cmd)
print("\macs2 BEDPE error stream (i.e. output):\n"+ret[1]+"\n")
if not isfile(input_BEDPE):
    raise Exception("ERROR: getting MACS2 input BAMPE file.")




