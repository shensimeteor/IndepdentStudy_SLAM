#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include <fstream>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include "burgers_4dvar.h"
#include <string>
#include <vector>

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::ResidualBlockId;
using namespace std;



int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    
    Burgers bg("burgers_conf.in", "model.burgers");
    cout<< "[log] Read config from burgers_conf.in" << endl;
    bg.outputConfig("burgers_model.out");
    
    Burgers_4Dvar bg4dvar("burgers_conf.in", "da_4dvar_experiment");
    cout<< "[log] Read config from burgers_conf.in" << endl;

    double** xb0 = bg.readXs(bg4dvar.xb0_file.c_str(), 1, bg.getNx()); 
    cout<< "[log] Read Xb0 from " << bg4dvar.xb0_file << endl;
//    double** Xs = &Xs_all[bgslam.spinup_steps];

    Observations obss;
    obss.load(bgslam.obs_file.c_str());
    cout<< "[log] Read Obs from " << bgslam.obs_file << ";  n_obs = " << obss.nobs << endl;
   
    ObsTimeGrouper obstg(bg4dvar.window_steps, bg4dvar.da_steps_unit, 0, 0, 0);
    obstg.group(obss);

    if(!bg4dvar.weak_constraint){
        CovModel* covB0 = bg4dvar.inflate(bg4dvar.readB0(bg.getNx()), bg4dvar.B0_inflate_factor);
        int nw = covB0->n_mode; //dimension of w
        double* w = new double[nw];
        
        Problem problem;
        CostFunction* cost_function_bckg = CostFunctorWb0::createAutoDiffCostFunction(nw);
        problem.AddResidualBlock(cost_function_bckg, NULL, w);
    }

    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = bg4dvar.max_num_iterations;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
//    bg.writeXs(bgslam.xas_file.c_str(), Xs, bgslam.window_steps, bg.getNx());
    cout<< "[log] Finish Burgers 4DVar" << endl;
    return 0;
}
