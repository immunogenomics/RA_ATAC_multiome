import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename

def getOneShellOutput(cmd):
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]


arg_obj = ArgumentParser(description='output highest scoring peak for each set of overlapping peaks')
arg_obj.add_argument('peak', metavar = 'peak', help='peak file; must be sorted -k1,1 -k2,2n; score in 5th column')
arg_obj.add_argument('-outFile', metavar = 'outFile', help='outFile file of non-overlapping peaks; peak with highest score kept')
args = arg_obj.parse_args()


if not isfile(args.peak):
    raise Exception("Peak file must exist. Please try again.")
    
ret = getOneShellOutput('sort -c -k1,1 -k2,2n %s'%(args.peak))
if len(ret[1])!=0 and 'disorder' in ret[1]:
    raise Exception("Peak file must be sorted -k1,1 -k2,2n")

if args.outFile is not None:
    outFile = args.outFile
else:
    outFile = args.peak.rsplit('.',1)[0]+'_iterative.'+args.peak.rsplit('.',1)[1]

IN = open(args.peak,'r')
OUT = open(outFile,'w')
for i,line in enumerate(IN):
    line = line.strip()
    tabs = line.split('\t')
    
    chrom = tabs[0]
    start = int(tabs[1])
    stop = int(tabs[2])
    score = int(tabs[4])
    
    if i==0:
        prevChrom = chrom
        prevStop = stop
        
        bestScore = score
        bestLine = line
        
        continue
    
    #= includes bookends
    if chrom==prevChrom and start<=prevStop:
        if score>bestScore:
            bestScore=score
            bestLine=line
    else:
        OUT.write(bestLine+'\n')
        
        #reset
        bestScore = score
        bestLine = line
    
    
    prevChrom = chrom
    prevStop = stop

IN.close()

#last one
OUT.write(bestLine+'\n')

OUT.close()


