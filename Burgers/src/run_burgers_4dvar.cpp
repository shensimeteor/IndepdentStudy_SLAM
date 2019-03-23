#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include <stdio.h>
#include <fstream>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include "burgers_4dvar.h"
#include <string>
#include <vector>

using ceres::Jet;
using ceres::DynamicAutoDiffCostFunction;
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
    cout<< "[log] Read config from burgers_conf.in -- model.burgers" << endl;
    bg.outputConfig("burgers_model.out");
    
    Burgers_4DVar bg4dvar("burgers_conf.in", "da_4dvar_experiment");
    cout<< "[log] Read config from burgers_conf.in -- da_4dvar_experiment" << endl;

    double** xb0 = bg.readXs(bg4dvar.xb0_file.c_str(), 1, bg.getNx()); 
    cout<< "[log] Read Xb0 from " << bg4dvar.xb0_file << endl;
//    double** Xs = &Xs_all[bgsl.spinup_steps];

    Observations obss;
    obss.load(bg4dvar.obs_file.c_str());
    cout<< "[log] Read Obs from " << bg4dvar.obs_file << ";  n_obs = " << obss.nobs << endl;
   
    ObsTimeGrouper obstg(bg4dvar.window_steps, bg4dvar.da_steps_unit, 0, 0, 0);
    obstg.group(obss);
    cout<< "[log] Group Obs into "<<obstg.nt_obs<<" groups; n_obs = "<<obstg.n_obs<<endl;

    Problem problem;
    CovModel* covB0 = bg4dvar.inflate(bg4dvar.readB0(bg.getNx()), bg4dvar.B0_inflate_factor);
    cout<< "[log] Read CovB0" << endl;
    int nw = covB0->n_mode; //dimension of w
    double* w = new double[nw];
    for(int i=0; i<nw; i++){
        w[i] = 0.0;
    }

    //background term
    CostFunction* cost_function_bckg = CostFunctorWb0::createDynamicAutoDiffCostFunction(nw);
    vector<double*> parameter_blocks;
    parameter_blocks.clear();
    parameter_blocks.push_back(w);
    problem.AddResidualBlock(cost_function_bckg, NULL, parameter_blocks);


    //if(!bg4dvar.weak_constraint){
        //DynamicAutoDiffCostFunction<CostFunctorWb0, 4>* cost_function_bckg = 
        //    new DynamicAutoDiffCostFunction<CostFunctorWb0, 4>(new CostFunctorWb0(nw));
        //CostFunction* cost_function_bckg = CostFunctorWb0::createAutoDiffCostFunction(nw);
        //cost_function_bckg->AddParameterBlock(nw);
        //cost_function_bckg->SetNumResiduals(nw);
        CostFunction* cost_function_obs = CostFunctor_4DVar_FullObs::createDynamicAutoDiffCostFunction(covB0, &obstg, &bg, xb0[0]);
        problem.AddResidualBlock(cost_function_obs, NULL, parameter_blocks);
//    }

    cout<< "[log] Problem Defined" << endl; 
    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = bg4dvar.max_num_iterations;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
//    bg.writeXs(bgslam.xas_file.c_str(), Xs, bgslam.window_steps, bg.getNx());
    cout<< "[log] Finish Optimization and get W" << endl;

    //postprocess: w->xa, M(xa) -> xas
    covB0->updateXb0byW(w, xb0[0]);
    bg.init(xb0[0]);
    bg.advanceNStepsAndOutputBin(bg4dvar.window_steps, bg4dvar.xas_file.c_str(),1,1);
    cout<< "[log] Finish Re-run model and output xas" <<endl;
    return 0;
}
