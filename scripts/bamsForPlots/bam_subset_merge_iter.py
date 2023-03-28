import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname
from os import listdir
import re
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


arg_obj = ArgumentParser()
arg_obj.add_argument('listing', metavar = 'listing', help='listing file of SAM/BAM files to subset and then merge')
arg_obj.add_argument('region', metavar = 'region', help='chr:start-stop')
arg_obj.add_argument('tmpDir', metavar = 'tmpDir', help='temporary directory')
arg_obj.add_argument('outFile', metavar = 'outFile', help='file to be outputted; .bam only')
arg_obj.add_argument('-genome_file', metavar = 'genome_file', help='genome chrom sizes file to verify region is legit')
args = arg_obj.parse_args()


if not isfile(args.listing) or not isdir(dirname(args.outFile)) or not isdir(dirname(args.tmpDir[:-1])):
    raise Exception('Input file(s) must exist. Please try again.')


x = re.match(r"(chr.*?):([0-9]+)\-([0-9]+)",args.region)
if x is None:
    raise Exception("Region must be in the format of chr:start-stop")
chrom,start,stop = x.groups()
if int(start)<0:
    raise Exception("Region start can't be less than zero.")

if args.genome_file is not None:
    IN = open(args.genome_file,'r')
    for i,line in enumerate(IN):
        tabs=line.strip().split('\t')
        if tabs[0]==chrom and int(stop)>int(tabs[1]):
            IN.close()
            print(tabs[0],chrom,tabs[1],stop)
            raise Exception("Region stop can't be greater than the chromosome stop.")
    IN.close()
#else just assuming it's file.


if args.tmpDir[-1]!='/':
    args.tmpDir += '/'

tmpDir = args.tmpDir
mkdir(args.tmpDir)


to_concat=[]
IN=open(args.listing,'r')
for i,line in enumerate(IN):
    line = line.strip()
    if not isfile(line) or (line.rsplit('.',1)[-1]!='bam' and line.rsplit('.',1)[-1]!='sam'):
        raise Exception("ERROR: File doesn't exist or isn't a BAM/SAM: %s"%(line))
    
    if not isfile(line+'.bai'):
        cmd = "samtools index %s"%(line)
        ret = getOneShellOutput(cmd)
        if len(ret[1])!=0 or not isfile(line+'.bai'):
            raise Exception("ERROR: getting bam index file: %s"%(line+'.bai'))
    
    if i==0:
        #make the header file! (only need to do this once or else I'd use the samtools view -h option)
        header_file = tmpDir+'header.sam'
        cmd = 'samtools view -H %s > %s'%(line,header_file)
        ret = getOneShellOutput(cmd)
        if len(ret[1]) != 0 or not isfile(header_file):
            raise Exception("ERROR: getting header file: %s"%(ret[1]))
    
    #subset each file by the region
    listing_outFile = tmpDir+basename(line).rsplit('.',1)[0]+'_%s_%s_%s.bam'%(chrom,start,stop)
    cmd = 'samtools view -b %s %s > %s'%(line,args.region,listing_outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1]) != 0 or not isfile(listing_outFile):
        raise Exception("ERROR: getting listing output file: %s\n%s\n"%(listing_outFile,ret[1]))
    to_concat.append(listing_outFile)
    
IN.close()

#merge all the individually subsetted files
cmd = "samtools merge -cp -O BAM %s %s"%(args.outFile," ".join(to_concat))
print(time.strftime("%c"))
ret = getOneShellOutput(cmd)
print(time.strftime("%c"))
print(ret[1])
if not isfile(args.outFile):
    raise Exception("ERROR: cat command failed: %s"%(ret[1]))

#do the bam index file as well!
cmd = "samtools index %s"%(args.outFile)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(args.outFile+'.bai'):
    print("ERROR: getting bam index file: %s"%(args.outFile+'.bai'))

#clean up temp directory
cmd = 'rm %s %s'%(" ".join(to_concat),header_file)
ret = getOneShellOutput(cmd)
if len(listdir(tmpDir))==0:
    cmd = 'rmdir %s'%(tmpDir)
    ret = getOneShellOutput(cmd)

