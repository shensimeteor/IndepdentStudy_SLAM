#!/usr/bin/env python

import sys
import numpy as np

def help():
    print('''
usage:  python ShortSummaryDiff.py [*.bin 1] [*.bin 2] [nx,default 1000]
    ''')
    exit(0)


if(len(sys.argv[1:]) in [2,3] ):
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    if(len(sys.argv[1:]) == 3):
        nx = int(sys.argv[3])
    else:
        nx = 1000
else:
    help()

Xs1=np.fromfile(file1, dtype="double", count=-1)
Nt = int(len(Xs1) / nx)
Xs1=Xs1.reshape((Nt, nx))

Xs2=np.fromfile(file2, dtype="double", count=-1).reshape((Nt,nx))

diff = Xs1 - Xs2  # Nt, nx
rmse = np.sqrt(np.mean(diff*diff, axis = 1))

print("Time AVG rmse: "+str(np.mean(rmse)))
print("Time MAX rmse: "+str(np.max(rmse)))

npart = 10
part_size = int(np.ceil(Nt / npart))
for i in range(npart):
    idx_start = i*part_size
    idx_end = (i+1)*part_size if i < npart-1 else Nt
    print(str(i)+"th 10% part, AVG rmse:" + str(np.mean(rmse[idx_start:idx_end])))



