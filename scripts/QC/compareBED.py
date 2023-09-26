import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname

def getOneShellOutput(cmd):
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    if len(err_stream)!=0:
        print("ERROR:\n%s"%(err_stream.decode("UTF-8")))
    return [out_stream.decode("UTF-8").strip(),err_stream.decode("UTF-8")]

def isInt(x):
    try:
        int(x)
        float(x)
    except ValueError:
        return False
        
    if float(x)!=int(x):
        return False
    
    return True

def mkdir(d):
    if not isdir(d):
        cmd = "mkdir %s"%(d)
        ret = getOneShellOutput(cmd)
        if len(ret[1])!=0 or not isdir(d):
            raise Exception("ERROR: making output directory. %s"%(d))

def getCount(fl):
    cmd = "wc -l %s | awk '{print $1}'"%(fl)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isInt(ret[0]):
        raise Exception("ERROR: count of merge file: %s"%(fl))
    
    return int(ret[0])

def merge(x,outDir):
    outFile = outDir+basename(x).rsplit('.',1)[0]+'_merged.bed'
    cmd = "sort -k1,1 -k2,2n %s | mergeBed -i - > %s"%(x,outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(outFile):
        raise Exception("ERROR: merging file failed: %s"%(x))

    ct = getCount(outFile)

    print("%s all: %s"%(basename(x),ct))
    return ct, outFile

def only(f1,f2,outDir):
    #f1 and f2 should be merged already
    outFile = outDir+basename(f1).rsplit('.',1)[0]+'_only.bed'
    cmd = "intersectBed -a %s -b %s -wao | grep -P '%s0$' | cut -f1,2,3 > %s"%(f1,f2,r"\t",outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(outFile):
        raise Exception("ERROR: only file failed: %s"%(ret[1]))

    ct = getCount(outFile)

    print("%s only: %s"%(basename(f1),ct))
    return ct, outFile

def both(f1,f2,outDir):
    #f1 and f2 should be merged already
    outFile = outDir+basename(f1).rsplit('.',1)[0]+'_both.bed'
    cmd = "intersectBed -a %s -b %s -wa -u > %s"%(f1,f2,outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(outFile):
        raise Exception("ERROR: both file failed: %s"%(ret[1]))

    ct = getCount(outFile)

    print("%s both: %s"%(basename(f1),ct))
    return ct, outFile

def bp_both(f1,f2):
    cmd = 'intersectBed -a %s -b %s -wao | rev | cut -f1 | rev | awk \'BEGIN{OFS=FS="%s";s=0} {s+=$1} END{print s}\''%(f1,f2,r"\t")
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isInt(ret[0]):
        raise Exception("ERROR: bp_both failed: %s"%(ret[1]))
    
    print("bp both: %s"%(ret[0]))
    return int(ret[0])

def bp_ct(f1):
    cmd = 'awk \'BEGIN{OFS=FS="%s";s=0} {s+=($3-$2)} END{print s}\' %s'%(r"\t",f1)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isInt(ret[0]):
        raise Exception("ERROR: bp_both failed: %s"%(ret[1]))
    
    print("%s bp: %s"%(basename(f1),ret[0]))
    return int(ret[0])


arg_obj = ArgumentParser(description='compare overlap between two bed files in both base pairs and peaks')
arg_obj.add_argument('bed1', metavar = 'bed1', help='bed file 1')
arg_obj.add_argument('bed2', metavar = 'bed2', help = 'bed file 2')
arg_obj.add_argument('outDir', metavar = 'outDir', help = 'output directory')
arg_obj.add_argument('-outFile', metavar = 'outFile', help = 'output file - stats')
args = arg_obj.parse_args()

if not isfile(args.bed1) or not isfile(args.bed2):
    raise Exception("both input files need to exist.")

if args.outDir[-1]!='/':
    args.outDir+='/'
mkdir(args.outDir)

f1_merge_ct, f1_merge_file  = merge(args.bed1,args.outDir)
f2_merge_ct, f2_merge_file = merge(args.bed2,args.outDir)
f1_only_ct, f1_only_file = only(f1_merge_file,f2_merge_file,args.outDir)
f2_only_ct, f2_only_file = only(f2_merge_file,f1_merge_file,args.outDir)
f1_both_ct, f1_both_file = both(f1_merge_file,f2_merge_file,args.outDir)
f2_both_ct, f2_both_file = both(f2_merge_file,f1_merge_file,args.outDir)

bp_both_f1 = bp_both(f1_merge_file,f2_merge_file)
bp_both_f2 = bp_both(f2_merge_file,f1_merge_file)
if bp_both_f1 != bp_both_f2:
    print("WARNING: The shared number of bp is not symmetric.")

f1_bp_ct = bp_ct(f1_merge_file)
f2_bp_ct = bp_ct(f2_merge_file)


if args.outFile is None:
    args.outFile = args.outDir+basename(args.outDir[:-1])+'_stats.txt'
    
OUT = open(args.outFile,'w')
OUT.write("\t".join(map(str,["peak",basename(args.bed1),f1_both_ct,f1_merge_ct,int(round(f1_both_ct/f1_merge_ct*100,0)),f1_both_ct/f1_merge_ct*100]))+'\n')
OUT.write("\t".join(map(str,["peak",basename(args.bed2),f2_both_ct,f2_merge_ct,int(round(f2_both_ct/f2_merge_ct*100,0)),f2_both_ct/f2_merge_ct*100]))+'\n')
OUT.write("\t".join(map(str,["bp",basename(args.bed1),bp_both_f1,f1_bp_ct,int(round(bp_both_f1/f1_bp_ct*100,0)),bp_both_f1/f1_bp_ct*100]))+'\n')
OUT.write("\t".join(map(str,["bp",basename(args.bed2),bp_both_f2,f2_bp_ct,int(round(bp_both_f2/f2_bp_ct*100,0)),bp_both_f2/f2_bp_ct*100]))+'\n')
OUT.close()




