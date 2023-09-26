import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import os
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
    
    if base_string is not None:
        temp_prefix = temp_dir + base_string + '_'
    else:
        temp_prefix = temp_dir
    
    fileTmp= temp_prefix + hashlib.md5((time.ctime() + str(random())).encode('utf-8')).hexdigest()
    while isfile(fileTmp):
        fileTmp= temp_prefix + hashlib.md5((time.ctime() + str(random())).encode('utf-8')).hexdigest()
    
    call('touch %s'%fileTmp,shell=True)
    
    return fileTmp

def makeFragment(working):
    if len(working)!=2:
        print("IGNORING QNAME: %s not exactly 2 instances!"%(working[0][3]))
        return None
    
    #Tn5 9bp trimming (soft clipping aware)
    x = re.search(r"^([0-9]+)S",working[0][7])
    if x is not None:
        if int(x.group(1))<4:
            working[0][1] = int(working[0][1]) + (4-int(x.group(1)))
    else:
        working[0][1] = int(working[0][1])+4
    
    x = re.search(r"([0-9]+)S$",working[1][7])
    if x is not None:
        if int(x.group(1))<5:
            working[1][2] = int(working[1][2]) - (5-int(x.group(1)))
    else:
        working[1][2] = int(working[1][2])-5
    
    return "\t".join(map(str,[working[1][0],working[0][1],working[1][2],working[1][12][5:],working[1][3],int((int(working[0][4])+int(working[1][4]))/2)]))+'\n'


"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"


arg_obj = ArgumentParser(description='make fragment file from deduplicated BAM/SAM file')
arg_obj.add_argument('bam', metavar = 'bam', help='deduplicated bam (or sam) file; QNAME_ct!=2 discarded')
arg_obj.add_argument('outFile', metavar = 'outFile', help='output fragment file')
arg_obj.add_argument('-exclude_dupFrag', action='store_true', help="Exclude duplicate fragments if given. Default=False - include duplicate fragments")
arg_obj.add_argument('-qname_error', action='store_true', help="if given, will error out at the end if there are any IGNORNING QNAME messages because of a qname in the bam file that does not have exactly 2 instances.")
args = arg_obj.parse_args()


"""Validate arguments"""
if not isfile(args.bam):
    raise Exception("BAM file must exist.")


"""Convert BAM->BED"""
print(time.strftime("%c"))
temp = GenerateTempFilename(base_string = basename(args.bam).rsplit('.',1)[0])
temp2 = GenerateTempFilename(base_string = basename(args.bam).rsplit('.',1)[0])
toClean_files = [temp,temp2]

cmd = "samtools view -b %s | bam2bed - -d > %s"%(args.bam,temp2)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0:
    raise Exception("bam2bed failed.\nERROR: %s"%(ret[1]))
if not isfile(temp2) or os.stat(temp2).st_size == 0:
    raise Exception("temp2 file does not exist or is empty: %s"%(temp2))

cmd = "sort -k4,4 -k1,1 -k2,2n %s > %s"%(temp2,temp)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0:
    raise Exception("sort failed.\nERROR: %s"%(ret[1]))
if not isfile(temp) or os.stat(temp).st_size == 0:
    raise Exception("temp file does not exist or is empty: %s"%(temp))
print(time.strftime("%c"))


mated = None
intermediate_outFile = args.outFile.rsplit('.',1)[0]+"_withDup.txt"
toClean_files.append(intermediate_outFile)
OUT = open(intermediate_outFile,'w')
IN = open(temp,'r')
for i,line in enumerate(IN):
    line = line.strip()
    tabs = line.split("\t")
    
    CB = [p for p in tabs if len(p)>3 and p[0:4]=='CB:Z'][0]
    
    if i==0:
        mated = [tabs[0:12]+[CB]]
    elif tabs[3]==prevQN:
        mated.append(tabs[0:12]+[CB])
    else:
        toPrint = makeFragment(mated)
        if toPrint is not None:
            OUT.write(toPrint)
        else:
            if args.qname_error:
                OUT.close()
                IN.close()
                getOneShellOutput("rm %s"%(" ".join(toClean_files)))
                raise Exception('This file contained at least one IGNORING QNAME warning, so erroring.')
            
        #refresh
        mated = [tabs[0:12]+[CB]]
    
    prevQN=tabs[3]
IN.close()

if mated is None:
    raise Exception("last mated is None.")
toPrint = makeFragment(mated)
if toPrint is not None:
    OUT.write(toPrint)
else:
    if args.qname_error:
        OUT.close()
        getOneShellOutput("rm %s"%(" ".join(toClean_files)))
        raise Exception('This file contained at least one IGNORING QNAME warning, so erroring.')

OUT.close()
print(time.strftime("%c"))


fragment_file = args.outFile
if not args.exclude_dupFrag:
    cmd='cut -f1-4 %s | sort | uniq -c | awk \'BEGIN{OFS="%s"} {print $2,$3,$4,$5,$1}\' > %s'%(intermediate_outFile,r"\t",fragment_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0:
        raise Exception("Getting unique fragment/CB file.\nERROR: %s"%(ret[1]))
else:
    cmd = 'sort -k1,1 -k2,2n -k3,3n -k4,4 %s | awk \'BEGIN{OFS=FS="%s"} {key=$1"-"$2"-"$3"-"$4; if(NR==1) {bestScore=$6;bestLine=$0;prevKey=key}; if(key==prevKey) {if(score>bestScore) {bestScore=$6;bestLine=$0} } else {print bestLine;bestScore=$6;bestLine=$0}; prevKey=key} END{print bestLine}\' | sort -k1,1 -k2,2n > %s'%(intermediate_outFile, r"\t",fragment_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(fragment_file):
        raise Exception("ERROR: getting fragment file:\n%s"%(ret[1]))
    print(time.strftime("%c"))
    
getOneShellOutput("rm %s"%(" ".join(toClean_files)))

