#!/usr/bin/env python
import numpy as np
import sys

base = 0.2
nx=1000
xs = np.ones((nx,), np.double)*base
xs[600:900] = 0.3 

xs.tofile("xb0_err_stdv_array.bin")

