#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <libconfig.h++>
#include "burgers_slam.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
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
