import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename
import re
import time

def getOneShellOutput(cmd):
    print(cmd)
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]


arg_obj = ArgumentParser(description='concatenate multiple BAM files into one BAM file and make a corresponding index file')
arg_obj.add_argument('regions', metavar = 'regions', help='file of regions to grep for chromosome; chr:start-stop per line; order matters here!')
arg_obj.add_argument('inDir', metavar = 'inDir', help='input directory of per chromosome deduplicated bam files of format: prefix+chr+".bam"')
arg_obj.add_argument('outFile', metavar = 'outFile', help='file to be outputted; concatenated bam file; .bam only')
arg_obj.add_argument('prefix', metavar = 'prefix', help='prefix used in the input directory')
args = arg_obj.parse_args()

if not isfile(args.regions) or not isdir(args.inDir):
    raise Exception('Input directory and regions file must exist. Please try again.')


to_concat=[]
IN = open(args.regions,'r')
for i,line in enumerate(IN):
    line = line.strip()
    
    x = re.match(r"(chr.*?):([0-9]+)\-([0-9]+)",line)
    if x is None:
        raise Exception("Region must be in the format of chr:start-stop")
    chrom,start,stop = x.groups()
    
    chr_out = args.inDir+args.prefix+"_"+chrom+'.bam'
    if not isfile(chr_out):
        raise Exception("%s: output file not generated. Exitting."%(chr_out))
    to_concat.append(chr_out)
    
IN.close()

cmd = "samtools merge -cp -O BAM %s %s"%(args.outFile," ".join(to_concat))
print(time.strftime("%c"))
ret = getOneShellOutput(cmd)
print(time.strftime("%c"))
print(ret[1])
if not isfile(args.outFile):
    raise Exception("cat command failed.\nERROR: %s"%(ret[1]))

cmd = 'samtools index %s'%(args.outFile)
ret = getOneShellOutput(cmd)
if not isfile(args.outFile+'.bai'):
    raise Exception("index command failed.\nERROR: %s"%(ret[1]))


