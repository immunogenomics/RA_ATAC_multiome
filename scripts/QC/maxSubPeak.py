import argparse
from argparse import ArgumentParser
from os.path import isfile
import re
import decimal


arg_obj = ArgumentParser(description='determine MACS2 subpeak with the best peak signal')
arg_obj.add_argument('summits', metavar = 'summits', help='summits or peak file from MACS2')
arg_obj.add_argument('-peaks', action='store_true', help='narrowPeak file from MACS2')
arg_obj.add_argument('-outFile', metavar = 'outFile', help='outFile file')
args = arg_obj.parse_args()

if not isfile(args.summits):
    raise Exception("Input file need to exist.")

summits_file = args.summits

if args.outFile is not None:
    bestSummit_file = args.outFile
else:
    bestSummit_file = summits_file.rsplit('.',1)[0]+'_bestSubPeak.'+summits_file.rsplit('.',1)[1]
OUT = open(bestSummit_file,'w')

prevNum=-1
IN = open(summits_file,'r')
for i,line in enumerate(IN):
    tabs = line.strip().split("\t")

    if args.peaks:
        peakSignal = float(decimal.Decimal(tabs[8]))
    else:
        peakSignal = float(decimal.Decimal(tabs[4]))

    x = re.search(r"peak_([0-9]+)",tabs[3])
    peakNum = x.groups(1)[0]

    if peakNum == prevNum:
        if peakSignal > currSignal:
            currSignal = peakSignal
            currPrint = line.strip()
    else:
        if i!=0:
            OUT.write(currPrint+'\n')
        currSignal = peakSignal
        currPrint = line.strip()

    prevNum = peakNum
OUT.write(currPrint+'\n')
IN.close()
OUT.close()


