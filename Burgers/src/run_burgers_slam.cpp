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

    double* xb0; 
    if (!bgslam.xb0_file.empty()) {
        xb0 = bg.readX(bgslam.xb0_file.c_str(), bg.getNx());
        cout<< "[log] Read xb0 (separately) from "<<bgslam.xb0_file << endl;
        cout<< " -- xbs provides initial guess, xb0 provides initial background constraint" <<endl;
    } else {
        xb0 = Xs[0];
    }
    
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
    cout<< "[log] model-error term is added" << endl;
    
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
    cout<< "[log] obs-error term is added" << endl;
    
    // -3. initial background constraint (optional)
    int cnt_resiblock_x0 = 0; 
    int cnt_resiblocK_hybrid_weak = 0; 
    if (bgslam.option_hybrid_4dvar == "") { // no hybrid
        if(bgslam.xb0_constraint){
            double* xb0_err_stdv_array; // xb0 error stdv array [Nx]
            if (!bgslam.xb0_err_stdv_file.empty()) {
                xb0_err_stdv_array = bg.readX(bgslam.xb0_err_stdv_file.c_str(), nx);
            } else {
                xb0_err_stdv_array = new double[nx];
                for (int i=0; i<nx; i++) {
                    xb0_err_stdv_array[i] = bgslam.xb0_err_stdv;
                }
            }
            for(int i=0; i<nx; i++){
                double xb0_err_var = pow(xb0_err_stdv_array[i],2);
                CostFunction* cost_functionX0 = CostFunctorX0::createAutoDiffCostFunction(xb0_err_var, xb0[i]);
                ResidualBlockId id = problem.AddResidualBlock(cost_functionX0, NULL, &Xs[0][i]);
                resblock_ids.push_back(id);
                cnt_resiblock_x0++;
            }
            cout <<"[log] xb0 constraint term is added -- diagonal only" << endl;
        }
    } else { // hybrid 4dvar+slamda
        CovModel* covB0 = bgslam.inflate(bgslam.readB0(bg.getNx()), bgslam.B0_inflate_factor);
        int nw = covB0->n_mode;
        // ||w0||^2 term
        CostFunction* cost_functionW0 = CostFunctorWb0::createDynamicAutoDiffCostFunction(nw);
        double* w0 = new double[nw]; // w0 
        for (int i=0; i<nw; i++) {
            w0[i] = 0.0;
        }
        vector<double*> paras1;
        paras1.clear();
        paras1.push_back(w0);
        ResidualBlockId id = problem.AddResidualBlock(cost_functionW0, NULL, paras1);
        resblock_ids.push_back(id);
        cnt_resiblock_x0 = 1;
        cout <<"[log] xb0 constraint term is added -- w0 from hybrid-4dvar" << endl;
        if (bgslam.option_hybrid_4dvar == "weak") {
            // weak constraint term, ||x0-xb0-E0w0||^2
            for (int i=0; i<nx; i++) {
                CostFunction* cost_functionW0X0_weak = CostFunctorX0W0_weak::createDynamicAutoDiffCostFunction(bgslam.weak_hybrid_4dvar_wstd, covB0, xb0[i], i);
                vector<double*> paras2;
                paras2.clear();
                paras2.push_back(&Xs[0][i]);
                paras2.push_back(w0);
                ResidualBlockId id = problem.AddResidualBlock(cost_functionW0X0_weak, NULL, paras2);
                resblock_ids.push_back(id);
                cnt_resiblocK_hybrid_weak = 1;
            }
            cout <<"[log] hybrid-4dvar weak constraint term is added" << endl;
        } else {
            cout << "option_hybrid_4dvar = "<<bgslam.option_hybrid_4dvar<<" unsupported!" <<endl;
            exit(1);
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
        rbgroup.add_groups_byComponent(cnt_resiblock_proc, cnt_resiblock_obs, cnt_resiblock_x0, cnt_resiblock_xneighbor, cnt_resiblocK_hybrid_weak);
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
