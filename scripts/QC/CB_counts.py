import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import re
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


arg_obj = ArgumentParser(description='counting cell barcodes from a BAM file')
arg_obj.add_argument('bam', metavar = 'bam', help='bam')
arg_obj.add_argument('outFile', metavar = 'outFile', help='output file; cell barcode, count')
args = arg_obj.parse_args()

if not isfile(args.bam):
    raise Exception("BAM file must exist.")

"""Convert BAM->SAM for given region"""
print(time.strftime("%c"))
temp = GenerateTempFilename()
cmd = "samtools view %s > %s"%(args.bam,temp)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0:
    raise Exception("Samtools view failed.\nERROR: %s"%(ret[1]))
print(time.strftime("%c"))

dic = {}

IN = open(temp,'r')
for i,line in enumerate(IN):
    tabs = line.strip().split('\t')
    if 'CB:Z' not in line:
        continue
    CB = [p for p in tabs if len(p)>3 and p[0:4]=='CB:Z'][0][5:]
    
    if CB in dic.keys():
        dic[CB]+=1
    else:
        dic[CB]=1
IN.close()
print(time.strftime("%c"))
getOneShellOutput("rm %s"%(temp))

print(time.strftime("%c"))
OUT = open(args.outFile,'w')
for k,v in dic.items():
    OUT.write("%s\t%s\n"%(k,v))
OUT.close()
print(time.strftime("%c"))
    
    


