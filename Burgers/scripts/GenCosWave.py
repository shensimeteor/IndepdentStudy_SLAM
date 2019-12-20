#!/usr/bin/env python
import numpy as np
import sys
def help():
    print '''
<usage> Python GenCosWave.py <output file> <nx> <amplitude> <n-wave> <offset, in degree>
- offset: e.g. to gen sin wave, offset=-90, no matter <n_wave>
'''
    exit(0)

if __name__ == "__main__":
    if (len(sys.argv) < 6):
        help()
    outfile = sys.argv[1]
    nx = int(sys.argv[2])
    amp = float(sys.argv[3])
    nwave = int(sys.argv[4])
    offset = float(sys.argv[5])

    x=np.arange(nx)
    w=2*np.pi 
