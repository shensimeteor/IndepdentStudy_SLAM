#!/usr/bin/env python

import os
import numpy as np
import glob

# concat run_dir/*/output.bin's istep to concat_file_path
def concat_ensemble_astep(run_dir, istep, output_file_name, concat_file_path, nx):
    files=glob.glob(run_dir+"ens/ens_*/"+output_file_name)
    files.sort()
    print(istep)
    nfile = len(files)
    concat_data = np.zeros((nfile, nx), dtype="double")
    for i, filename in enumerate(files):
        Xs=np.fromfile(filename, dtype="double", count=-1)
        nt = int(len(Xs) / nx)
        Xs=Xs.reshape((nt, nx))
        concat_data[i,:] = Xs[istep,:]
    concat_data.tofile(concat_file_path)


steps=288
n_ens=200
cwd=os.getcwd()
template_burgers_dir=cwd+"/burgers_template/"
command_burgers="python "+cwd+"/../../PyODE/burgers.py dopri5"
run_dir=cwd+"/B288_dopri5/"
x_init_file=cwd+"/x_init/x_init_000.bin"
concat_file_path=run_dir+"concat_ens_step288.bin"

for iens in range(n_ens):
    print(iens)
    mem_dir=run_dir+"/ens_%0.5d" %iens
    if(os.path.isdir(mem_dir)):
        continue
    os.makedirs(mem_dir)
    os.chdir(mem_dir)
    os.symlink(x_init_file,"x_init_base.bin")
    os.system("AddGaussianNoise.py x_init_base.bin x_init.bin 0 0.05")
    os.remove("x_init_base.bin")
    os.system(command_burgers + " x_init.bin %d" %steps)
#concat
concat_ensemble_astep(run_dir, 0, "output.bin", run_dir+"concat_ens_step000.bin", 1000)
concat_ensemble_astep(run_dir, 1, "output.bin", run_dir+"concat_ens_step001.bin", 1000)
concat_ensemble_astep(run_dir, 3, "output.bin", run_dir+"concat_ens_step003.bin", 1000)
concat_ensemble_astep(run_dir, 6, "output.bin", run_dir+"concat_ens_step006.bin", 1000)
concat_ensemble_astep(run_dir, 36, "output.bin", run_dir+"concat_ens_step036.bin", 1000)
concat_ensemble_astep(run_dir, 144, "output.bin", run_dir+"concat_ens_step144.bin", 1000)
concat_ensemble_astep(run_dir, 288, "output.bin", run_dir+"concat_ens_step288.bin", 1000)
