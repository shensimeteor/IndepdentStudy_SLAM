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
 --   template.conf:  the config file template (for burgers model & slam)
 
 - obs_gen/: Python code
 --   FlexGenObs.py: code to generate obs (from model state)
 --   fgobs_conf.json:  obs config file

 - scripts/: Python code, some useful tools
 --   ShortSummaryDiff.py: summarize diff (rmse) between 2 model state bin files
 --   ExtractAPeriod.py: extract 1 or more time slice/slices from a bin file

 - build/: to build burgers model
 
 - build_ceres/: to build burgers slam

 - test/: common test cases

3. how to run:
 - compile: 
 --  cd build,  make clean,  make
 --  cd build_ceres, rm CMakeCache.txt,  cmake .,  make

 - run burgers model:
 --  cd test/normal_run, cp ../../src/template.conf burgers_conf.in,  <modify burgers_conf.in>,  ./run_burgers.x

 - run burgers slam:
 --  prepare obs:  cd test/obs_gen_test/,  cp ../../obs_gen/fgobs_conf.json ., <modify fgobs_conf.json>,  FlexGenObs.py <parameters>
 --  extract the window size from normal_run: use ExtractAPeriod.py
 --  run slam:  modify burgers_conf.in, ./run_burgers_slam.x
