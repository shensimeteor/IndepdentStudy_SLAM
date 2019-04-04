#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include <fstream>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include "burgers_slam.h"
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
    
    Burgers_SLAM bgslam("burgers_conf.in", "da_experiment");
    cout<< "[log] Read config from burgers_conf.in" << endl;

    ResiBlockGroup_Config rbgroup_conf("burgers_conf.in", "da_experiment_cost_groups");
    
    double** Xs_all = bg.readXs(bgslam.xbs_file.c_str(), bgslam.window_steps + bgslam.spinup_steps, bg.getNx());
    cout<< "[log] Read Xbs from " << bgslam.xbs_file << endl;
    double** Xs = &Xs_all[bgslam.spinup_steps];

    Observations obss;
    obss.load(bgslam.obs_file.c_str());
    cout<< "[log] Read Obs from " << bgslam.obs_file << ";  n_obs = " << obss.nobs << endl;
    
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
        if(it > bgslam.window_steps) 
            break;
        double obs_value=obss.obs[i].value;
        double obs_stdv =obss.obs[i].error;
        CostFunction* cost_functionObs = CostFunctorObs::createAutoDiffCostFunction(ix, obs_value, obs_stdv*obs_stdv); 
        ResidualBlockId id = problem.AddResidualBlock(cost_functionObs, NULL, &Xs[it][ix]);
        resblock_ids.push_back(id);
        cnt_resiblock_obs++;
    }
    
    // -3. initial background constraint (optional)
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
    // -4. neighbor constraint (optional)
    int cnt_resiblock_xneighbor = 0;
    if(bgslam.xneighbor_constraint){
        double xneighbor_diff_variance = bgslam.xneighbor_diff_stdv * bgslam.xneighbor_diff_stdv;
        int tend = bgslam.window_steps;
        if(bgslam.xneighbor_only_x0) {
            tend = 0;
        }
        for(int t=0; t<=tend; t++){
            for(int i=0; i<nx; i++){
                CostFunction* cost_functionXnbr;
                if(bgslam.xneighbor_constraint_xa_dxa == "xa"){
                    cost_functionXnbr = CostFunctorXneighbor::createAutoDiffCostFunction(xneighbor_diff_variance);
                }else{
                    cost_functionXnbr = CostFunctorXneighbor::createAutoDiffCostFunction(xneighbor_diff_variance, Xs[t][(i-1+nx)%nx], Xs[t][i], Xs[t][(i+1)%nx]);
                }
                ResidualBlockId id = problem.AddResidualBlock(cost_functionXnbr, NULL, &Xs[t][i], &Xs[t][(i-1+nx)%nx], &Xs[t][(i+1)%nx]);
                resblock_ids.push_back(id);
                cnt_resiblock_xneighbor++;
            }
        }
    }
    cout<< "[log] Problem Defined" << endl; 

    //evaluate before solve
    ResiBlockGroup rbgroup;
    if(rbgroup_conf.output_cost_groups){
        //set groups & subgroups
        rbgroup.clear_groups();
        rbgroup.add_groups_byComponent(cnt_resiblock_proc, cnt_resiblock_obs, cnt_resiblock_x0, cnt_resiblock_xneighbor);
        if(rbgroup_conf.output_proc_subgroups){
            rbgroup.add_groups_proc_subgroup(cnt_resiblock_proc, rbgroup_conf.proc_subgroup_size, 0);
        }
        if(rbgroup_conf.output_obs_subgroups){
            rbgroup.add_groups_obs_subgroup(cnt_resiblock_obs, rbgroup_conf.obs_subgroup_size, cnt_resiblock_proc);
        }
        if(rbgroup_conf.output_xneighbor_subgroups){
            rbgroup.add_groups_xneigbor_subgroup(cnt_resiblock_xneighbor, rbgroup_conf.xneighbor_subgroup_size, cnt_resiblock_proc+cnt_resiblock_obs);
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
