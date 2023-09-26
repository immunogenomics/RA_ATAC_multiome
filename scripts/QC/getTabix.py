import argparse
from argparse import ArgumentParser
from os.path import isfile, isdir, basename, dirname, realpath
import time


def getOneShellOutput(cmd):
    print(cmd)
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    return [out_stream.strip().decode("UTF-8"),err_stream.decode("UTF-8")]


def mkdir(x):
    if not isdir(x):
        cmd = "mkdir %s"%(x)
        ret = getOneShellOutput(cmd)
        print(ret)
        
        if not isdir(x):
            return False
        else:
            return True
    else:
        return True


"""Input Arguments"""
arg_obj = ArgumentParser(description='sort, bgzip, tabix a fragment file')
arg_obj.add_argument('inFile',metavar = 'inFile', help='input fragment file')
arg_obj.add_argument('outFile', metavar = 'outFile', help='output sorted fragment file - .gz and .gz.tbi will be added within')
args = arg_obj.parse_args()


"""Verifying Input Arguments
A lot of this is done in the ArgumentParser"""
if not isfile(args.inFile):
    raise Exception("Input file must exist.")


"""Printing Arguments"""
print("Argument List:")
args_print = vars(args)
for k in args_print:
    print("%s = %s"%(k,args_print[k]))
print('\n')

print(time.strftime("%c"))
cmd = "sort -k1,1 -k2,2n %s > %s"%(args.inFile, args.outFile)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(args.outFile):
    raise Exception("sort failed.\nERROR: %s"%(ret[1]))

print(time.strftime("%c"))
bgzip_outFile = args.outFile+'.gz'
cmd = "bgzip %s"%(args.outFile)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(bgzip_outFile):
    raise Exception("bgzip failed.\nERROR: %s"%(ret[1]))

print(time.strftime("%c"))
tabix_outFile = args.outFile+'.gz.tbi'
cmd = "tabix -p bed %s"%(bgzip_outFile)
ret = getOneShellOutput(cmd)
if len(ret[1])!=0 or not isfile(tabix_outFile):
    raise Exception("tabix failed.\nERROR: %s"%(ret[1]))

print(time.strftime("%c"))
