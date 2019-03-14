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

        


    // Build the problem.
    Problem problem;
    vector<ResidualBlockId> resblock_ids;
    
    // -1. Process Contraints
    int nx = bg.getNx();
    int cnt_resiblock_proc=0;
    double proc_err_var = bgslam.proc_err_stdv * bgslam.proc_err_stdv; 
    for( int i=1; i<=bgslam.window_steps; i++){
        for(int j=0; j<nx; j++){
            if(i==1 || bg.getNumericOption() == NUMERIC_OPTION_FORWARD){
                CostFunction* cost_functionProc = CostFunctorProc_fwd::createAutoDiffCostFunction(proc_err_var, i, &bg);
                ResidualBlockId id = problem.AddResidualBlock(cost_functionProc, NULL, &Xs[i-1][(j-1+nx)%nx], &Xs[i-1][j], &Xs[i-1][(j+1)%nx], &Xs[i][j]);
                resblock_ids.push_back(id);
            } else {
                CostFunction* cost_functionProc = CostFunctorProc_frogleap::createAutoDiffCostFunction(proc_err_var, i, &bg);
                ResidualBlockId id = problem.AddResidualBlock(cost_functionProc, NULL, &Xs[i-2][j], &Xs[i-1][(j-1+nx)%nx], &Xs[i-1][j], &Xs[i-1][(j+1)%nx], &Xs[i][j]);
                resblock_ids.push_back(id); 
            }
            cnt_resiblock_proc++;
        }
    }
    
    // -2. Observation Contraints
    int cnt_resiblock_obs = 0;
    for(int i=0; i< obss.nobs; i++){
        int ix=obss.obs[i].ix;
        int it=obss.obs[i].it;
        double obs_value=obss.obs[i].value;
        double obs_stdv =obss.obs[i].error;
        CostFunction* cost_functionObs = CostFunctorObs::createAutoDiffCostFunction(ix, obs_value, obs_stdv*obs_stdv); 
        ResidualBlockId id = problem.AddResidualBlock(cost_functionObs, NULL, &Xs[it][ix]);
        resblock_ids.push_back(id);
        cnt_resiblock_obs++;
    }
    
    int cnt_resiblock_x0 = 0;
    if(bgslam.xb0_constraint){
        double xb0_err_var = bgslam.xb0_err_stdv * bgslam.xb0_err_stdv; 
        for(int i=0; i<nx; i++){
            CostFunction* cost_functionX0 = CostFunctorX0::createAutoDiffCostFunction(xb0_err_var, Xs[0][i]);
            ResidualBlockId id = problem.AddResidualBlock(cost_functionX0, NULL, &Xs[0][i]);
            resblock_ids.push_back(id);
            cnt_resiblock_x0++;
        }
    }
    cout<< "[log] Problem Defined" << endl; 

    //evaluate before solve
    ResiBlockGroup rbgroup;
    if(rbgroup_conf.output_cost_groups){
        //set groups & subgroups
        rbgroup.clear_groups();
        rbgroup.add_groups_byComponent(cnt_resiblock_proc, cnt_resiblock_obs, cnt_resiblock_x0);
        if(rbgroup_conf.output_proc_subgroups){
            rbgroup.add_groups_proc_subgroup(cnt_resiblock_proc, rbgroup_conf.proc_subgroup_size, 0);
        }
        if(rbgroup_conf.output_obs_subgroups){
            rbgroup.add_groups_obs_subgroup(cnt_resiblock_obs, rbgroup_conf.obs_subgroup_size, cnt_resiblock_proc);
        }
        //evaluate all 
        rbgroup.evaluate_all(problem, resblock_ids);
        //sum by groups
        rbgroup.sum_metrics_bygroup();
        rbgroup.output_cost_bygroup(rbgroup_conf.xb_cost_filename.c_str());
        cout<< "[log] Cost Groups for xb is generated: " << rbgroup_conf.xb_cost_filename << endl;
    }
    
    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = bgslam.max_num_iterations;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    //evaluate after solve
    if(rbgroup_conf.output_cost_groups){
        rbgroup.evaluate_all(problem, resblock_ids);
        rbgroup.sum_metrics_bygroup();
        rbgroup.output_cost_bygroup(rbgroup_conf.xa_cost_filename.c_str());
        cout<< "[log] Cost Groups for xa is generated: " << rbgroup_conf.xa_cost_filename << endl;
    }

    std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    bg.writeXs(bgslam.xas_file.c_str(), Xs, bgslam.window_steps, bg.getNx());
    cout<< "[log] Finish SLAM Burgers" << endl;
    return 0;
}
