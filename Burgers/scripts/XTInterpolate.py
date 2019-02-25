#!/usr/bin/env python

import sys
import numpy as np

def help():
    print('''
usage:  python XTInterpolate.py [input.bin] [output.bin] [input NX] [input NT] [output Nx] [output Nt]
interpolate the [Nt,Nx] 2d domain of input, to [out_Nt, out_Nx] 2d domain (same domain, but different grid)
    ''')
    exit(0)


if(len(sys.argv[1:]) == 6):
    in_bin = sys.argv[1]
    out_bin = sys.argv[2]
    in_nx = int(sys.argv[3])
    in_nt = int(sys.argv[4])
    out_nx = int(sys.argv[5])
    out_nt = int(sys.argv[6])
else:
    help()

inXs=np.fromfile(in_bin, dtype="double", count=-1).reshape((in_nt,in_nx))
# to do

in_x_len = len(in_nx)


end_t = start_t + len_t - 1
if (end_t > Nt):
    print("Warning: end_t = %d, exceeds the length of actual Nt")
    end_t = Nt

outputXs = Xs1[start_t:end_t+1,:]
print("start_t = %d, length_t = %d, end_t = %d" %(start_t, end_t + 1 -start_t, end_t))
outputXs.tofile(outbin)



