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

def isInt(x):
    try:
        int(x)
        float(x)
    except ValueError:
        return False
        
    if float(x)!=int(x):
        return False
    
    return True


"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"

"""Required files"""
blacklist_file = homeDir+'downloadedInputs/hg38-blacklist.v2.bed'
if not isfile(blacklist_file):
    raise Exception("Blacklist file must exist: %s"%(blacklist_file))
chrom_bed_file = homeDir+'downloadedInputs/hg38_normal_chromosomes.bed'
if not isfile(chrom_bed_file):
    raise Exception("Chromosome bed file must exist: %s"%(chrom_bed_file))
    

arg_obj = ArgumentParser(description='statistics comparing pre- and post-deduplication BAMs')
arg_obj.add_argument('bam', metavar = 'bam', help='bam (or sam) file')
arg_obj.add_argument('dedupBam', metavar = 'dedupBam', help='deduplicated bam (or sam) file')
arg_obj.add_argument('outFile', metavar = 'outFile', help='file to be outputted; table of stats')
arg_obj.add_argument('-MAPQ_cutoff', metavar = 'MAPQ_cutoff', type=int, default = 60, help="MAPQ cutoff score. Default=60.")
arg_obj.add_argument('-blacklist', metavar = 'blacklist', default=blacklist_file, help="Blacklist BED file to exclude from resulting BAM file. Default = Boyle Lab hg38 v2")
arg_obj.add_argument('-chrom_bed', metavar = 'chrom_bed', default=chrom_bed_file, help="chr1-22XY BED file to subset from BAM files.")
args = arg_obj.parse_args()

if not isfile(args.bam) or not isfile(args.dedupBam):
    raise Exception('Both BAM (or SAM) files must exist. Please try again.')

if args.chrom_bed is not None and not isfile(args.chrom_bed):
    raise Exception("If given, chrom_bed must exist.")

"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n\n')

descriptions = []
unfilt = []
filt = []
for i,f in enumerate([args.bam,args.dedupBam]):
    
    print("\nStarting: %s"%(f))
    
    print(time.strftime("%c"))
    #convert to BAM file if not already
    if f.rsplit('.',1)[1]!='bam':
        input_BAM = f.rsplit('.',1)[0]+'.bam'
        cmd = "samtools view -b %s > %s"%(f,input_BAM)
        ret = getOneShellOutput(cmd)
        if len(ret[1])!=0 or not isfile(input_BAM):
            raise Exception("ERROR: converting SAM to BAM.\n%s"%(ret[1]))
    else:
        input_BAM = f
    
    print(time.strftime("%c"))
    #flagstat - sam or bam
    cmd = "samtools flagstat %s"%(input_BAM)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or len(ret[0])==0:
        raise Exception("ERROR: flagstat %s:\n%s"%(input_BAM,ret[1]))
    
    for i_for,line in enumerate(ret[0].split('\n')):
        line = line.strip()
        
        x = re.match(r"([0-9]+) \+ ([0-9]+) (.*)",line)
        if x is None:
            raise Exception("ERROR: parsing flagstat.")
        pass_ct,fail_ct,description = x.groups()
        
        if fail_ct!="0":
            print("fail count for this row is 0!:\n%s"%(line))
        
        if i==0:
            descriptions.append(description)
            unfilt.append(int(pass_ct))
        else:
            filt.append(int(pass_ct))
    
    
    print(time.strftime("%c"))
    #MAPQ - sam or bam
    cmd = "samtools view -q %s -c %s"%(args.MAPQ_cutoff,input_BAM)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isInt(ret[0]):
        raise Exception("ERROR: MAPQ cutoff %s:\n%s"%(input_BAM,ret[1]))
    
    if i==0:
        descriptions.append("MAPQ>=%s"%(args.MAPQ_cutoff))
        unfilt.append(int(ret[0]))
    else:
        filt.append(int(ret[0]))
    
    
    print(time.strftime("%c"))
    #CB - if sam, don't need the samtools view, though be aware of header
    cmd = "samtools view %s | grep -vc 'CB:Z'"%(input_BAM)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isInt(ret[0]):
        raise Exception("ERROR: CB presence %s:\n%s"%(input_BAM,ret[1]))
    
    if i==0:
        descriptions.append("NO cell barcode")
        unfilt.append(int(ret[0]))
    else:
        filt.append(int(ret[0]))
    
    
    print(time.strftime("%c"))
    #blacklist - reads that overlap - BAM, not SAM.
    if args.chrom_bed is None:
        cmd = "intersectBed -abam %s -b %s | samtools view -c"%(input_BAM,args.blacklist) 
    else:
        cmd = "samtools view -b %s -L %s | intersectBed -abam stdin -b %s | samtools view -c"%(input_BAM,args.chrom_bed,args.blacklist)
    ret = getOneShellOutput(cmd)
    
    if len(ret[1])!=0:
        print("ERROR stream:\n%s"%(ret[1]))
    if not isInt(ret[0]):
        raise Exception("ERROR: blacklist filtering %s"%(input_BAM))
    
    if i==0:
        descriptions.append("Overlapping blacklist")
        unfilt.append(int(ret[0]))
    else:
        filt.append(int(ret[0]))
    
    print(time.strftime("%c"))


#result formatting
OUT = open(args.outFile,'w')
delim=","
header_1 = ["category",basename(args.bam),"",basename(args.dedupBam),"",""]
header_2 = ["","QC-passed reads","percent","QC-passed reads","percent","percent of unfiltered"]
OUT.write(delim.join(header_1)+'\n')
OUT.write(delim.join(header_2)+'\n')

for idx in range(len(descriptions)):
    OUT.write(delim.join(map(str,[descriptions[idx],unfilt[idx],unfilt[idx]/unfilt[0],filt[idx],filt[idx]/filt[0],filt[idx]/unfilt[0]]))+"\n")

OUT.close()

