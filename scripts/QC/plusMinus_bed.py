import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename

def isInt(x):
    try:
        int(x)
        float(x)
    except ValueError:
        return False

    if float(x)!=int(x):
        return False

    return True


arg_obj = ArgumentParser(description='manipulate peak lengths')
arg_obj.add_argument('inFile', metavar = 'inFile', help='BED file')
arg_obj.add_argument('outFile', metavar = 'outFile', help='BED with plus/minus X amount')
arg_obj.add_argument('flank', metavar = 'flank', type=float, help='flank size')
arg_obj.add_argument('-perc', action='store_true', help='flank is percentage of region length')
arg_obj.add_argument('type', metavar = 'type', help='start, end, both, center (or summit if given in column 10)')
arg_obj.add_argument('-keepBody', action='store_true', help='keep gene body if type is start or end.')
arg_obj.add_argument('bounds', metavar = 'bounds', help='genome bounds - chrom.sizes file')
args = arg_obj.parse_args()

if not isfile(args.inFile) or not isfile(args.bounds) or args.type not in ['start','end','both','center']:
    raise Exception('Input errors. Please try again.')

if args.keepBody and args.type not in ['start','end']:
    raise Exception("keep body is only meaningful when type is start or end.")

bounds = {}
IN = open(args.bounds,'r')
for i,line in enumerate(IN):
    tabs=line.strip().split("\t")
    
    if not isInt(tabs[1]):
        raise Exception("Bounds file second column not int.")
    
    bounds[tabs[0]]=int(tabs[1])
IN.close()

IN = open(args.inFile,'r')
OUT = open(args.outFile,'w')
for i,line in enumerate(IN):
    tabs = line.strip().split('\t')
    
    if not isInt(tabs[1]) or not isInt(tabs[2]):
        raise Exception("inFile file second or third column not int.")
    
    if args.perc:
        flank = int(args.flank*(int(tabs[2])-int(tabs[1])))
    else:
        flank = int(args.flank)
    
    if args.type=="start":
        if (len(tabs)<=5) or (len(tabs)>5 and (tabs[5]=="+" or tabs[5]=='.')):
            x = int(tabs[1])
            newStart = x-flank if x-flank >= 0 else 0
            if args.keepBody:
                newEnd = int(tabs[2])
            else:
                newEnd = x
        elif len(tabs)>5 and tabs[5]=="-":
            x=int(tabs[2])
            newEnd = x+flank if x+flank <= bounds[tabs[0]] else bounds[tabs[0]]
            if args.keepBody:
                newStart = int(tabs[1])
            else:
                newStart = x
        else:
            raise Exception("Shouldn't have gotten here.")
    elif args.type=="end":
        if (len(tabs)<=5) or (len(tabs)>5 and (tabs[5]=="+" or tabs[5]=='.')):
            x = int(tabs[2])
            newEnd = x+flank if x+flank <= bounds[tabs[0]] else bounds[tabs[0]]
            if args.keepBody:
                newStart = int(tabs[1])
            else:
                newStart = x
        elif len(tabs)>5 and tabs[5]=="-":
            x=int(tabs[1])
            newStart = x-flank if x-flank >= 0 else 0
            if args.keepBody:
                newEnd = int(tabs[2])
            else:
                newEnd = x
        else:
            raise Exception("Shouldn't have gotten here.")
    elif args.type=="both":
        newStart = int(tabs[1])-flank if int(tabs[1])-flank >= 0 else 0
        newEnd = int(tabs[2])+flank if int(tabs[2])+flank <= bounds[tabs[0]] else bounds[tabs[0]]
    elif args.type=="center":
        if len(tabs)==10 and isInt(tabs[9]):
            #summit
            if tabs[5]=="+" or tabs[5]=='.':
                x = tabs[1]
            if tabs[5]=="-":
                x=tabs[2]
            center = int(x)+int(tabs[9])
        else:
            #center
            center = (int(tabs[1])+int(tabs[2]))/2
        newStart = center-flank if center-flank >=0 else 0
        newEnd = center+flank if center+flank <= bounds[tabs[0]] else bounds[tabs[0]]
    else:
        raise Exception("Shouldn't have gotten here. type is start, end, both, or center.")
    
    OUT.write("\t".join(map(str,[tabs[0],newStart,newEnd]+tabs[3:]))+'\n')
OUT.close()
IN.close()

