#!/usr/bin/env python
import numpy as np
import sys
def help():
    print('''
    AddGaussianNoise.py  <input bin file> <output name> <noise mean> <noise stdv>
    ''')
    exit(0)



if(len(sys.argv[1:]) == 4):
    filebin = sys.argv[1]
    outbin = sys.argv[2]
    mean = float(sys.argv[3])
    stdv = float(sys.argv[4])
else:
    help()

xs=np.fromfile(filebin, dtype="double", count=-1)
n=xs.shape[0]
noise = np.random.normal(mean, stdv, n)
xs += noise
xs.tofile(outbin)

