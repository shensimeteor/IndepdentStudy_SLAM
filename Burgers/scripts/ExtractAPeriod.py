#!/usr/bin/env python

import sys
import numpy as np

def help():
    print('''
usage:  python ExtractAPeriod.py [input.bin] [output.bin] [start t, 0-based] [length t, default 1] [nx,default 1000] 
    ''')
    exit(0)


if(len(sys.argv[1:]) in [3,4,5] ):
    filebin = sys.argv[1]
    outbin = sys.argv[2]
    start_t = int(sys.argv[3])
    if(len(sys.argv[1:]) >= 4):
        len_t = int(sys.argv[4])
    else:
        len_t = 1
    if(len(sys.argv[1:]) == 5):
        nx = int(sys.argv[5])
    else:
        nx = 1000
else:
    help()

Xs1=np.fromfile(filebin, dtype="double", count=-1)
Nt = int(len(Xs1) / nx)
Xs1=Xs1.reshape((Nt, nx))


end_t = start_t + len_t - 1
if (end_t > Nt):
    print("Warning: end_t = %d, exceeds the length of actual Nt")
    end_t = Nt

outputXs = Xs1[start_t:end_t+1,:]
print("start_t = %d, length_t = %d, end_t = %d" %(start_t, end_t + 1 -start_t, end_t))
outputXs.tofile(outbin)



