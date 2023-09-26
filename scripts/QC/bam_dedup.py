import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import re

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

def isInt(x):
    try:
        int(x)
        float(x)
    except ValueError:
        return False
        
    if float(x)!=int(x):
        return False
    
    return True

def keepLine(flag):
    badIndices = [2,3,8,9,11] #SAM flag shouldn't have any of these bits set
    goodIndices = [0,1] #SAM flag should have all of these bits set
    notEqualIndices = [(4,5),(6,7)] #SAM flag should not have both bits in the pair set or both bits unset
    
    ignore = 0
    for b in badIndices:
        ignore |= (flag >> b & 1)
    for b in notEqualIndices:
        if (flag >> b[0] & 1) == (flag >> b[1] & 1):
            ignore |= 1
        else:
            ignore |= 0
    
    keep = 1
    for b in goodIndices:
        keep &= (flag >> b & 1)
    
    final = ~ignore & keep
    return final


"""Current script directory"""
homeDir = dirname(realpath(__file__))+"/"

"""Required files"""
blacklist_file = homeDir+'downloadedInputs/hg38-blacklist.v2.bed'
if not isfile(blacklist_file):
    raise Exception("Blacklist file must exist: %s"%(blacklist_file))


arg_obj = ArgumentParser(description='keep duplicates between cells, not within cells while maintaining read/mate pair status using read/mate start positions')
arg_obj.add_argument('bam', metavar = 'bam', help='bam file')
arg_obj.add_argument('region', metavar = 'region', help='format: chr:start-stop')
arg_obj.add_argument('outFile', metavar = 'outFile', help='output file; deduplicated bam file')
arg_obj.add_argument('-MAPQ_cutoff', metavar = 'MAPQ_cutoff', type=int, default = 60, help="MAPQ cutoff score. Default=60.")
arg_obj.add_argument('-include_unpaired', action='store_true', help="Include unpaired reads if given. Default=False - exclude unpaired reads")
arg_obj.add_argument('-blacklist', metavar = 'blacklist', default=blacklist_file, help="Blacklist BED file to exclude from resulting BAM file. Default = Boyle Lab hg38 v2")
args = arg_obj.parse_args()


"""Validate arguments"""
if not isfile(args.bam):
    raise Exception("BAM file must exist.")

MAPQ_cutoff = int(args.MAPQ_cutoff)
print("Using MAPQ cutoff of: %s"%(MAPQ_cutoff))

if args.include_unpaired:
    print("Including any unpaired QNAMEs.")
else:
    print("Excluding any unpaired QNAMEs.")

if not isfile(args.blacklist):
    raise Exception("Blacklist file must exist.")

x = re.match(r"(chr.*?):([0-9]+)\-([0-9]+)",args.region)
if x is not None:
    chrom,start,stop = x.groups()

    if int(stop)<int(start):
        raise Exception("The region stop must be greater than the region start.")
    if int(start)<0:
        raise Exception("The region start must be 0 or greater.")
    
    spacer = ' '
    regions = args.region
    tmp_base = chrom
elif isfile(args.region):
    regions_bed_file = GenerateTempFilename(base_string=basename(args.region).rsplit('.',1)[0]+'_bed')
    cmd='samtools view -H %s | awk \'NR==FNR{if(substr($2,1,3)=="SN:") {split($2,b,":"); a[b[2]]=b[2]}; next} {split($1,c,":"); if(c[1] in a) print $0}\' - %s | sed "s/:/%s/" | sed "s/-/%s/" > %s'%(args.bam,args.region,r"\t",r"\t",regions_bed_file)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0 or not isfile(regions_bed_file):
        raise Exception("Regions bed file failed.\nERROR: %s"%(ret[1]))
    
    spacer = ' -L '
    regions = regions_bed_file
    tmp_base = basename(args.region).rsplit('.',1)[0]
else:
    raise Exception("Region must be in the format of chr:start-stop or a BED file")


"""Convert BAM->SAM for given region
I don't know of a way to parse the BAM file directly."""
temp = GenerateTempFilename(base_string=tmp_base)
cmd = "samtools view -b %s%s%s | intersectBed -abam stdin -b %s -v | samtools view > %s"%(args.bam,spacer,regions,args.blacklist,temp)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0:
    raise Exception("Samtools view failed.\nERROR: %s"%(ret[1]))

if isfile(args.region):
    getOneShellOutput("rm %s"%(regions_bed_file))


"""If keeping a read, I should also keep its mate, as determined by QNAME and SAM flag
If the read is written to the output file, then its QNAME is read_mate_pairs"""
read_mate_pairs = set()

"""This keeps track of the SAM lines for the chosen reads"""
dic_qname = {}


"""Iterate through SAM file"""
firstLine_flag = True
dic = {} #dictionary flushed for each read/mate combo regardless of cell barcode
intermediate_outFile = args.outFile.rsplit('.',1)[0]+"_withUnpairedQNAMEs.sam"
OUT = open(intermediate_outFile,'w')
IN = open(temp,'r')
for i,line in enumerate(IN):
    line = line.strip()
    tabs = line.split("\t")
    score = int(tabs[4])
    
    
    """Initial determination if the SAM line is valid
    (i.e., no problematic SAM flags, not chimeric, has a cell barcode, and a good MAPQ score)"""
    # parsing SAM Flag
    if not keepLine(int(tabs[1])):
        continue
    if tabs[6]!="=":
        continue
    if 'CB:Z' not in line:
        continue
    #MAPQ == 255 means N/A
    if score < MAPQ_cutoff or score == 255:
        continue
    
    """Remaking score if read MAPQ score available"""
    if 'MQ:i' in line:
        MQ = [p[5:] for p in tabs if len(p)>5 and p[0:4]=='MQ:i'][0]
        if isInt(MQ):
            score = int((score+int(MQ))/2)
        
    
    """Handling QNAMEs
    I assume that the max count of QNAMEs after SAM flag filtering is 2.
    Therefore, I can just use a binary T/F for if the QNAME has been written to the output file."""
    qname = tabs[0]
    qname_flag = qname in read_mate_pairs
    
    
    """Initializing dictionary values"""
    CB = [p for p in tabs if len(p)>3 and p[0:4]=='CB:Z'][0]
    key = "-".join([tabs[2],tabs[3],tabs[7],CB])
    
    
    """curr is the chrom/start variable
    It determines when to flush the dictionary to the output file and refresh.
    I position sort the BAM file earlier in the pipeline to allow for this."""
    curr = "-".join([tabs[2],tabs[3]])
    if firstLine_flag:
        prev = curr
        firstLine_flag = False
    
    
    """Main loop of the entire thing."""
    if curr==prev:
        """If within same chrom/start group, then still in the process of choosing reads/mates"""
        
        """Only update scoring dictionary if:
        1) the qname_flag is set, meaning this is the mate of a read already written to the output file, OR
        2) this key has already been seen, but that previous read wasn't a mate (so not prioritized) and has a lower MAPQ score, OR
        3) this key hasn't been seen yet, so it is our de facto choice for the moment"""
        if qname_flag or (key in dic and not dic.get(key)[2] and dic.get(key)[0] < score) or key not in dic:
            dic[key] = [score,qname,qname_flag]
        
        """Down here for the read start==mate start"""
        if qname in dic_qname.keys():
            tmp = dic_qname.get(qname)
            tmp.append(line)
            dic_qname[qname]=tmp
        else:
            dic_qname[qname]=[line]
    else:
        """Once we've moved on to a new chrom/start group, we have to flush the dictionarys to the output file and refresh with the new group"""
        
        for k in dic.keys():
            
            """Print SAM lines to output through QNAMEs"""
            val = dic.get(k)
            chosen_qnames = dic_qname.get(val[1])
            OUT.write('\n'.join(chosen_qnames)+'\n')
            
            """Updating read_mate_pairs"""
            #Read was already written to output, now mate was also added, so we can disregard the QNAME from further consideration
            if val[1] in read_mate_pairs:
                read_mate_pairs.discard(val[1])
            #Only the read was written to output, the mate still needs to be added, so we need to keep an eye out for this QNAME
            elif len(chosen_qnames)==1:
                read_mate_pairs.add(val[1])
        
        """Refresh dictionarys with new chrom/start group values"""
        qname_flag = qname in read_mate_pairs #since we just updated read_mate_pairs
        
        dic = {}
        dic[key] = [score,qname,qname_flag]
        
        dic_qname = {}
        dic_qname[qname]=[line]
    
    """Update chrom/start variable for next line"""
    prev = curr
IN.close()

"""Close out the last chrom/start group to output"""
for k in dic:
    
    """Print SAM lines to output through QNAMEs"""
    val = dic.get(k)
    chosen_qnames = dic_qname.get(val[1])
    OUT.write('\n'.join(chosen_qnames)+'\n')
    
    """Updating read_mate_pairs"""
    #Read was already written to output, now mate was also added, so we can disregard the QNAME from further consideration
    if val[1] in read_mate_pairs:
        read_mate_pairs.discard(val[1])
    #Only the read was written to output, the mate still needs to be added, so we need to keep an eye out for this QNAME
    elif len(chosen_qnames)==1:
        read_mate_pairs.add(val[1])
    
OUT.close()

"""Remove SAM intermediate file"""
getOneShellOutput("rm %s"%(temp))

"""Are there any unpaired QNAMEs left at the end?"""
if len(read_mate_pairs) != 0:
    print("read_mate_pairs still has %s QNAMEs left."%(len(read_mate_pairs)))

if args.include_unpaired:
    cmd="mv %s %s"%(intermediate_outFile,args.outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0:
        raise Exception("Moving intermediate outFile to actual outFile.\nERROR: %s"%(ret[1]))
elif len(read_mate_pairs)!=0:
    remain_QNAME_file = args.outFile.rsplit('.',1)[0]+'_remainingQNAMEs.txt'
    with open(remain_QNAME_file, 'w') as f:
        f.write("\n".join(read_mate_pairs)+"\n")
        
    cmd = 'awk \'BEGIN{OFS=FS="%s"} FNR==NR{a[$1]=$1;next} {if(!($1 in a)) print $0}\' %s %s > %s'%(r"\t",remain_QNAME_file,intermediate_outFile,args.outFile)
    ret = getOneShellOutput(cmd)
    if len(ret[1])!=0:
        raise Exception("Excluding unpaired QNAMEs from intermediate outFile.\nERROR: %s"%(ret[1]))
    
    getOneShellOutput("rm %s %s"%(remain_QNAME_file,intermediate_outFile))


#make a sorted BAM file
bam_outFile = args.outFile.rsplit('.',1)[0]+'.bam'
cmd = "samtools view -H %s | cat - %s | samtools sort -O bam -o %s"%(args.bam,args.outFile,bam_outFile)
ret = getOneShellOutput(cmd)
if not isfile(bam_outFile):
    raise Exception("Making final bam file.\nERROR: %s"%(ret[1]))
getOneShellOutput("rm %s"%(args.outFile))
getOneShellOutput("samtools index %s"%(bam_outFile))


