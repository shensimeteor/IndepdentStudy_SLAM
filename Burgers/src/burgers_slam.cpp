#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <fstream>
#include <libconfig.h++>
#include "burgers_slam.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::ResidualBlockId;
using namespace std;
using namespace libconfig;

//Burgers_SLAM
Burgers_SLAM::Burgers_SLAM(const char* config_file, const char* config_path){
    Config cfg;
    try{
        cfg.readFile(config_file);
        Setting & root = cfg.getRoot();
        Setting& slamda = root.lookup(config_path);
        if(slamda.exists("spinup_steps"))
            this->spinup_steps = slamda["spinup_steps"];
        if(slamda.exists("window_steps"))
            this->window_steps = slamda["window_steps"];
        if(slamda.exists("xb0_constraint"))
            this->xb0_constraint = slamda["xb0_constraint"];
        if(slamda.exists("xb0_err_stdv"))
            this->xb0_err_stdv = slamda["xb0_err_stdv"];
        if(slamda.exists("proc_err_stdv"))
            this->proc_err_stdv = slamda["proc_err_stdv"];
        if(slamda.exists("max_num_iterations"))
            this->max_num_iterations = slamda["max_num_iterations"];
        slamda.lookupValue("xbs_file", this->xbs_file);
        slamda.lookupValue("obs_file", this->obs_file);
        slamda.lookupValue("xas_file", this->xas_file);
    }catch(ParseException &pex){
        cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                      << " - " << pex.getError() << endl;
        exit(1);
    }catch(FileIOException &fioex){
        cerr << "File I/O Error " << config_file << endl;
        exit(1);
    }catch(SettingNotFoundException &nfex){
        cerr << "Setting Not Found "<<nfex.getPath() << endl;
        exit(1);
    }
}


//CostFunctorObs
CostFunctorObs::CostFunctorObs(int iloc, double obs, double obs_error_variance){
    this->iloc = iloc;
    this->obs = obs;
    this->obs_error_variance = obs_error_variance;
    this->obs_error_sqrtinv  = sqrt(1./obs_error_variance);
}

template <typename T> bool CostFunctorObs::operator()(const T* const x, T* residual) const {
    residual[0] = obs_error_sqrtinv * (x[0] - obs);
    return true;
}

CostFunction* CostFunctorObs::createAutoDiffCostFunction(int iloc, double obs, double obs_error_variance){
    return new AutoDiffCostFunction<CostFunctorObs, 1, 1>(new CostFunctorObs(iloc, obs, obs_error_variance));
}


//CostFunctorProc: baseClass
CostFunctorProc::CostFunctorProc(double proc_error_variance, int t, Burgers* bg){
    this->proc_error_variance = proc_error_variance;
    this->proc_error_sqrtinv = sqrt(1./proc_error_variance);
    this->t = t;
    this->bg = bg;
    bg->getConvenient(&dtdx, &dtdx2, &c0, &c1);
}


//CostFunctorProc_fwd
template <typename T> bool CostFunctorProc_fwd::operator()(const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const {
    residual[0] = proc_error_sqrtinv*(xp1[0] - 0.5*dtdx*xp1[0]*(xp1p1[0]-xp1m1[0]) + 0.5*c1*(xp1m1[0]+xp1p1[0] - 2.*xp1[0]) - x[0]);
    return true;
}

CostFunction* CostFunctorProc_fwd::createAutoDiffCostFunction(double proc_error_variance, int t, Burgers* bg){
    return new AutoDiffCostFunction<CostFunctorProc_fwd, 1, 1, 1, 1, 1>(new CostFunctorProc_fwd(proc_error_variance, t, bg));
}

//CostFunctorProc_frogleap
template <typename T> bool CostFunctorProc_frogleap::operator()(const T* const xp2, const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const {
    residual[0] = proc_error_sqrtinv*(c0*(xp2[0] + c1*(xp1p1[0] + xp1m1[0] - xp2[0]) - dtdx*(xp1[0]*(xp1p1[0]-xp1m1[0]))) - x[0]);
    return true;
}

CostFunction* CostFunctorProc_frogleap::createAutoDiffCostFunction(double proc_error_variance, int t, Burgers* bg){
    return new AutoDiffCostFunction<CostFunctorProc_frogleap, 1,1,1,1,1,1>(new CostFunctorProc_frogleap(proc_error_variance, t, bg));
}


//CostFunctorX0
CostFunctorX0::CostFunctorX0(double init_error_variance, double init_priori){
    this->init_error_variance = init_error_variance;
    this->init_error_sqrtinv = sqrt(1./init_error_variance);
    this->init_priori = init_priori;
}

template <typename T> bool CostFunctorX0::operator()(const T* const x, T* residual) const{
    residual[0] = init_error_sqrtinv * (x[0] - init_priori);
    return true;
}

CostFunction* CostFunctorX0::createAutoDiffCostFunction(double init_error_variance, double init_priori){
    return new AutoDiffCostFunction<CostFunctorX0, 1, 1>(new CostFunctorX0(init_error_variance, init_priori));
}


//CostFunction Evaluator

void ResiBlockGroup::add_resi_block_group(string gname, int start_idx, int end_idx){
    group_name.push_back(gname);
    start_ids.push_back(start_idx);
    end_ids.push_back(end_idx);
}

void ResiBlockGroup::evaluate_all(Problem & problem, vector<ResidualBlockId> & resblock_ids){
    Problem::EvaluateOptions options;
    options.residual_blocks = resblock_ids;
    each_costs.clear();
    total_cost = 0;
    problem.Evaluate(options, &total_cost, &each_costs, NULL, NULL);
    for(int i=0; i<each_costs.size(); i++){
        each_costs[i] *= each_costs[i] * 0.5;
    }
    evaluated = true;
}

void ResiBlockGroup::sum_metrics_bygroup(){
    group_costs.clear();
    if(!evaluated){
        cout<<"Error, need to call ResiBlockGroup::evaluate_all first!"<<endl;
        exit(1);
    }
    for(int i=0; i<group_name.size(); i++){
        int start_idx = start_ids[i];
        int end_idx = end_ids[i];
        double sum=0;
        for(int j=start_idx; j<=end_idx; j++){
            sum+= each_costs[j];
        }
        group_costs.push_back(sum);
    }
}

void ResiBlockGroup::output_cost_bygroup(const char* csvname){
    ofstream f(csvname, ios::out);
    if(f){
        f<<"group_name, start_idx, end_idx, cnt_resiblock, cost \n";
        for(int i=0; i<group_name.size(); i++){
            f<<group_name[i]<<","<<start_ids[i]<<","<<end_ids[i]<<","<< (end_ids[i] - start_ids[i]+1) << "," << group_costs[i] <<"\n";
        }
    }else{
        cout<<"Error opening file: "<<csvname<<endl;
    }
}

void ResiBlockGroup::clear_groups(){
    group_name.clear();
    start_ids.clear(); 
    end_ids.clear();
    group_costs.clear();
}

//into groups: all_proc, all_obs, all_x0, 
void ResiBlockGroup::add_groups_byComponent(int cnt_resiblock_proc, int cnt_resiblock_obs, int cnt_resiblock_x0){
    add_resi_block_group("all_proc", 0, cnt_resiblock_proc-1);
    add_resi_block_group("all_obs", cnt_resiblock_proc, cnt_resiblock_proc+cnt_resiblock_obs-1);
    add_resi_block_group("all_X0", cnt_resiblock_obs+cnt_resiblock_proc, cnt_resiblock_proc+cnt_resiblock_obs+cnt_resiblock_x0-1);
}

void ResiBlockGroup::add_groups_proc_subgroup(int cnt_resiblock_proc, int subgroup_size, int offset){
    int ngroup = ceil(cnt_resiblock_proc *1.0/ subgroup_size);
    for(int i=0; i< ngroup; i++){
        int istart = i*subgroup_size + offset;
        int iend = offset + (i<ngroup-1? (i+1)*subgroup_size-1 : cnt_resiblock_proc-1);
        string name = string("proc_")+to_string(i+1);
        add_resi_block_group(name, istart, iend);
    }
}

//for the normal order of residual blocks: proc, obs, X0, here offset should be cnt_resiblock_proc
void ResiBlockGroup::add_groups_obs_subgroup(int cnt_resiblock_obs, int subgroup_size, int offset){
    int ngroup = ceil(cnt_resiblock_obs *1.0/ subgroup_size);
    for(int i=0; i< ngroup; i++){
        int istart = i*subgroup_size + offset;
        int iend = offset+ (i<ngroup-1? (i+1)*subgroup_size-1 : cnt_resiblock_obs-1) ;
        string name = string("obs_")+to_string(i+1);
        add_resi_block_group(name, istart, iend);
    }
}

ResiBlockGroup_Config::ResiBlockGroup_Config(const char* config_file, const char* config_path){
    Config cfg;
    try{
        cfg.readFile(config_file);
        Setting & root = cfg.getRoot();
        Setting& slamda = root.lookup(config_path);
        if(slamda.exists("output_cost_groups"))
            this->output_cost_groups = slamda["output_cost_groups"];
        if(this->output_cost_groups){
            if(slamda.exists("output_obs_subgroups"))
                this->output_obs_subgroups = slamda["output_obs_subgroups"];
            if(this->output_obs_subgroups){
                this->obs_subgroup_size = slamda["obs_subgroup_size"];
            }
            if(slamda.exists("output_proc_subgroups"))
                this->output_proc_subgroups = slamda["output_proc_subgroups"];
            if(this->output_proc_subgroups){
                this->proc_subgroup_size = slamda["proc_subgroup_size"];
            }
            slamda.lookupValue("xb_cost_filename", this->xb_cost_filename);
            slamda.lookupValue("xa_cost_filename", this->xa_cost_filename);
        }
    }catch(ParseException &pex){
        cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                      << " - " << pex.getError() << endl;
        exit(1);
    }catch(FileIOException &fioex){
        cerr << "File I/O Error " << config_file << endl;
        exit(1);
    }catch(SettingNotFoundException &nfex){
        cerr << "Setting Not Found "<<nfex.getPath() << endl;
        exit(1);
    }
}
