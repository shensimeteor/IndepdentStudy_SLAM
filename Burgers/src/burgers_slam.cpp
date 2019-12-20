#include "ceres/ceres.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <fstream>
#include <libconfig.h++>
#include "burgers_slam.h"

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
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
        if(slamda.exists("xb0_err_stdv_file"))
            slamda.lookupValue("xb0_err_stdv_file", this->xb0_err_stdv_file);
        if(slamda.exists("proc_err_stdv"))
            this->proc_err_stdv = slamda["proc_err_stdv"];
        if(slamda.exists("max_num_iterations"))
            this->max_num_iterations = slamda["max_num_iterations"];
        if(slamda.exists("xneighbor_constraint"))
            this->xneighbor_constraint = slamda["xneighbor_constraint"];
        if(slamda.exists("xneighbor_only_x0"))
            this->xneighbor_only_x0 = slamda["xneighbor_only_x0"];
        if(slamda.exists("xneighbor_constraint_xa_dxa"))
            slamda.lookupValue("xneighbor_constraint_xa_dxa", this->xneighbor_constraint_xa_dxa);
        if(slamda.exists("xneighbor_diff_stdv"))
            this->xneighbor_diff_stdv = slamda["xneighbor_diff_stdv"];
        if(slamda.exists("xb0_file"))
            slamda.lookupValue("xb0_file", this->xb0_file);
        if(slamda.exists("option_hybrid_4dvar"))
            slamda.lookupValue("option_hybrid_4dvar", this->option_hybrid_4dvar);
        if(slamda.exists("weak_hybrid_4dvar_wstd"))
            this->weak_hybrid_4dvar_wstd = slamda["weak_hybrid_4dvar_wstd"];
        if(slamda.exists("cov_dir"))
            slamda.lookupValue("cov_dir", this->cov_dir);
        if(slamda.exists("B0_inflate_factor"))
            this->B0_inflate_factor = slamda["B0_inflate_factor"];
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

CovModel* Burgers_SLAM::readB0(int nx){
    string B0file = this->cov_dir + string("/B0_modes.bin");
    CovModel* pcov = new CovModel(B0file.c_str(), nx);
    return pcov;
}

CovModel* Burgers_SLAM::inflate(CovModel* B0, double inflator){
    B0->inflate(inflator);
    return B0;
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


//CostFunctors for hybrid-4dvar-slamda
CostFunctorX0W0_weak::CostFunctorX0W0_weak(double wstd, CovModel* B0, double xb0, int xidx) {
    this->weak_stdv = wstd; 
    this->B0 = B0;
    this->xidx = xidx;
    this->xb0 = xb0;
}

template <typename T> bool CostFunctorX0W0_weak::operator()(T const* const* paras, T* residual) const {
    T const* x0 = paras[0]; // size 1
    T const* w0 = paras[1]; // size nw
    // sum = xb0 + E0w0
    T sum = (T)0.0;
    for (int j=0; j<this->B0->n_mode; j++) {
        sum = sum + w0[j] * B0->modes[j][xidx];
    }
    //  || (x0 - xz0)/wstd ||^2 
    residual[0] = (x0[0] - sum) / weak_stdv;
    return true;
}

CostFunction* CostFunctorX0W0_weak::createDynamicAutoDiffCostFunction(double wstd, CovModel* B0, double xb0, int xidx) {
    int w_size = B0->n_mode;
    CostFunctorX0W0_weak* functor = new CostFunctorX0W0_weak(wstd, B0, xb0, xidx);
    DynamicAutoDiffCostFunction<CostFunctorX0W0_weak, 4>* cost_function = new 
        DynamicAutoDiffCostFunction<CostFunctorX0W0_weak, 4>(functor);
    cost_function->AddParameterBlock(1);
    cost_function->AddParameterBlock(w_size);
    cost_function->SetNumResiduals(1);
    return cost_function;
}



//CostFunctorXneighbor
CostFunctorXneighbor::CostFunctorXneighbor(double neighbor_diff_variance, double xb0, double xb1, double xb2){
    this->neighbor_diff_variance = neighbor_diff_variance;
    this->neighbor_diff_sqrtinv = sqrt(1./neighbor_diff_variance);
    this->xb0 = xb0;
    this->xb1 = xb1;
    this->xb2 = xb2;
}

template <typename T> bool CostFunctorXneighbor::operator()(const T* const x1, const T* const x0, const T* const x2, T* residual) const{
    residual[0] = neighbor_diff_sqrtinv * ( (x1[0]-xb1) - (x0[0]-xb0+x2[0]-xb2)/2.0);
    return true;
}

CostFunction* CostFunctorXneighbor::createAutoDiffCostFunction(double neighbor_diff_variance, double xb0, double xb1, double xb2){
    return new AutoDiffCostFunction<CostFunctorXneighbor, 1, 1, 1, 1>(new CostFunctorXneighbor(neighbor_diff_variance, xb0, xb1, xb2));
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
void ResiBlockGroup::add_groups_byComponent(int cnt_resiblock_proc, int cnt_resiblock_obs, int cnt_resiblock_x0, int cnt_resiblock_xneighbor, int cnt_resiblock_hybrid_weak) {
    int sum=0;
    add_resi_block_group("all_proc", sum, sum+cnt_resiblock_proc-1);
    sum += cnt_resiblock_proc;
    add_resi_block_group("all_obs", sum, sum+cnt_resiblock_obs-1);
    sum += cnt_resiblock_obs;
    add_resi_block_group("all_X0", sum, sum+cnt_resiblock_x0-1);
    sum += cnt_resiblock_x0;
    add_resi_block_group("all_hybridw", sum, sum+cnt_resiblock_hybrid_weak-1);
    sum += cnt_resiblock_hybrid_weak;
    add_resi_block_group("all_xneighbor", sum, sum+cnt_resiblock_xneighbor-1);
    sum += cnt_resiblock_xneighbor;
    add_resi_block_group("all_all", 0, sum-1);
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

//for the normal order of residual blocks: proc, obs, X0, xneighbor, here offset should be cnt_resiblock_proc
void ResiBlockGroup::add_groups_obs_subgroup(int cnt_resiblock_obs, int subgroup_size, int offset){
    int ngroup = ceil(cnt_resiblock_obs *1.0/ subgroup_size);
    for(int i=0; i< ngroup; i++){
        int istart = i*subgroup_size + offset;
        int iend = offset+ (i<ngroup-1? (i+1)*subgroup_size-1 : cnt_resiblock_obs-1) ;
        string name = string("obs_")+to_string(i+1);
        add_resi_block_group(name, istart, iend);
    }
}


void ResiBlockGroup::add_groups_xneigbor_subgroup(int cnt_resiblock_xneighbor, int subgroup_size, int offset){
    int ngroup = ceil(cnt_resiblock_xneighbor *1.0/ subgroup_size);
    for(int i=0; i< ngroup; i++){
        int istart = i*subgroup_size + offset;
        int iend = offset+ (i<ngroup-1? (i+1)*subgroup_size-1 : cnt_resiblock_xneighbor-1) ;
        string name = string("xneighbor_")+to_string(i+1);
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
            if(slamda.exists("output_xneighbor_subgroups"))
                this->output_xneighbor_subgroups = slamda["output_xneighbor_subgroups"];
            if(this->output_xneighbor_subgroups){
                this->xneighbor_subgroup_size = slamda["xneighbor_subgroup_size"];
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
