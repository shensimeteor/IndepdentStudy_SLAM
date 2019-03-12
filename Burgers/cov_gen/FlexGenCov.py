#!/usr/bin/env python
import numpy as np
import pprint as pp
import csv
import json
import sys
import matplotlib.pyplot as plt


def help():
    print( '''
usage: FlexGenCov.py [fgcov_conf.in]
''')
    exit(0)

def genGaussianCorrMatrix(nx, is_periodic, ref_len):
    B = np.zeros((nx,nx), np.float64)
    for i in range(nx):
        for j in range(i, nx):
            dis  = (j-i)*1.0
            if(is_periodic):
                dis2 = (i+nx-j)*1.0
                if(dis2 < dis):
                    dis = dis2
            dis /= ref_len
            B[i,j] = np.exp(- (dis**2) / 2)
            B[j,i] = B[i,j]
    return B


#A = E E^T, by EigDecomposition, A = U S U^T, so E = U*S^0.5
#ratio_threshold, is for selecting modes
def MatrixDecomp(A, ratio_threshold=1.0):
    eigvalue,eigVectors = np.linalg.eig(A)
    eigvalue = np.real(eigvalue)  #only keep the real part
    eigVectors = np.real(eigVectors)
    sortIdx = np.flip(np.argsort(eigvalue))  #sort by eigvalue
    maxIdx = np.count_nonzero(eigvalue > 0)  # only keep the positive eigvalue
    n_positive_eigvalues = maxIdx
    eigvalue_positive = eigvalue[sortIdx[0:maxIdx]] 
    eigVectors_positive = eigVectors[:, sortIdx[0:maxIdx]]
    
    sum_eigvalue = np.sum(eigvalue_positive)
    if(ratio_threshold < 1.0):
        maxIdx = np.nonzero(np.cumsum(eigvalue_positive) >= ratio_threshold*sum_eigvalue)[0][0]
    print(maxIdx)
    n_selected_eigvalues = maxIdx
    eigvalue_selected = eigvalue_positive[0:maxIdx]
    eigVectors_selected = eigVectors_positive[:, 0:maxIdx]
    print("total: %d, positive: %d, selected: %d" %(A.shape[0], n_positive_eigvalues, n_selected_eigvalues))
    
    E = np.matmul(eigVectors_selected, np.diag(np.sqrt(eigvalue_selected)))
    # test
    print("decomposition error (max) = %f" %np.max(np.abs( A - E.dot(E.T))))
    return E


if(len(sys.argv[1:]) == 1):
    fileconf = sys.argv[1]
else:
    help()

res_obs=[]
with open(fileconf, "r") as f:
    dct = json.load(f)
    cnt=0
    for cov_set in dct["cov_config_array"]:
        cnt+=1
        print("cov %d --------------" %(cnt))
        if(cov_set["correlation_model"] == "gaussian"):
            Cov = genGaussianCorrMatrix(cov_set["nx"], cov_set["nx_periodic"], cov_set["gaussian_reflen_grid"]) * cov_set["stdv"]**2
            E = MatrixDecomp(Cov, cov_set["modes_eigvalue_ratio_threshold"])
        else:
            print("Only Guassian model is supported")
        #output
        if(cov_set["output_Cov_matrix"]):
            Cov.tofile(cov_set["output_Cov_matrix"])
        if(cov_set["output_Cov_modes"]):
            E.T.tofile(cov_set["output_Cov_modes"])
        #img
        if(cov_set["out_image_cov"]):
            plt.figure()
            plt.imshow(Cov)
            plt.title("Cov Model")
            plt.colorbar()
            plt.savefig(cov_set["out_image_cov"])
            plt.close()
        if(cov_set["out_image_cov_reconstructed"]):
            plt.figure()
            Cov_rec = E.dot(E.T)
            plt.imshow(Cov_rec)
            plt.title("Cov Model (Reconstructed)")
            plt.colorbar()
            plt.savefig(cov_set["out_image_cov_reconstructed"])
            plt.close()

        


