1. pre-requirement:
 - install ceres-solver and its dependencies
 - install libconfig (for C++)
  
2. dir/code structure: 
 - src/: C++ code
 --   Burgers.h/cpp:  Burgers class, for Burgers model setting, running, data input, output
 --   run_burgers.cpp: main code to run burgers model (run_burgers.x executable)
 --   Observation.h/cpp: for Observation read
 --   burgers_slam.h/cpp:  classes/structs to run burgers slam
 --   run_burgers_slam.cpp: main code to run burgers slam DA (run_burgers_slam.x executable)
 --   burgers_4dvar.h/cpp: classes/structs to run burgers 4dvar & weak-constraint 4dvar
 --   run_burgers_4dvar.cpp: main code to run burgers 4dvar/weak-constraint 4dvar DA (run_burgers_4dvar.x executable)
 --   template.conf:  the config file template (for burgers model & slam & 4dvar)
 --   test*.cpp:  ignore these files, these are for simple test
 
 - obs_gen/: Python code
 --   FlexGenObs.py: code to generate obs (from model state)
 --   fgobs_conf.json:  obs config file

 - cov_gen/: Python code
 --   FlexGenCov.py:  code to generate Covariance and its square-root (i.e. Modes) (for 4dvar run) 
 --   fgcov_conf.json: config file

 - fake_obs/: ignore that, used for generate obs from model state, not used now
 

 - scripts/: Python code, some useful tools. To see how to use it, just type <name_of_script> without any argument, it will printout help info
 --   ShortSummaryDiff.py: summarize diff (rmse) between 2 model state bin files
 --   ExtractAPeriod.py: extract 1 or more time slice/slices from a bin file
 --   QuickPlots.py: run inside a 4dvar or slam running directory (after 4dvar or slam finished), to generate "images/" directory
 --   LinkAllQt.py: link all needed Qt files to a Cov Modes file 
 --   AddGaussianNoise.py: add a gaussian noise (independent with each other) to a model state file (*.bin), and output
 --   QuickPlotATimeStep.py: Just plot one time step

 - build/: to build burgers model
 
 - build_ceres/: to build burgers slam

 - test/: simple test cases. ignore

 - expt/: all experiments, check its README.txt to see more

3. how to compile & run:
   see: HOW_TO_COMPILE_RUN.txt
