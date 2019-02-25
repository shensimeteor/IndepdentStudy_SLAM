Burgers Model tests:

normal_run: normal setting
fwd_run: same with normal, but change frogleap to fwd numeric solution, diff rmse: 0.002
hires_run: same with normal, but increase x/t resolution by 5 times, diff rmse: 0.03
perturbed_run: same with normal, but add N(0,0.05^2) noise in each step, diff rmse:
perturbed_run_replicate: to test if perturbed_run can be replicated, (by recording the seed from burgers_conf.out)
test_restart: to test of restart/start a model with a given initial condition file. 

Burgers SLAM DA tests:
slam_test/: w constraints xb0, e0 = 0.05
slam_test2: w/o xb0 constraints
slam_test3: w constraints xb0, e0 = 0.2
slam_test4: w constraints xb0, e0 = 1
