#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <libconfig.h++>
#include <string>

using namespace std;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

//for config read 
class Burgers_SLAM {
public:
    int window_steps=288;
    string xbs_file;
    int spinup_steps=0;
    string obs_file;
    string xas_file;
    bool xb0_constraint = false;
    double xb0_err_stdv = 0.05;
    double proc_err_stdv = 0.05;
    int max_num_iterations = 10;

    Burgers_SLAM() {}
    Burgers_SLAM(const char* config_file, const char* config_path);
};



//must be on a grid (iloc)
struct CostFunctorObs {
    double obs, obs_error_variance, obs_error_sqrtinv;
    int iloc;
    CostFunctorObs(int iloc, double obs, double obs_error_variance);
    template <typename T> bool operator()(const T* const x, T* residual) const;
    static CostFunction* createAutoDiffCostFunction(int iloc, double obs, double obs_error_variance);
};

//process
struct CostFunctorProc {
    double proc_error_variance, proc_error_sqrtinv;
    int t; //0 means init condition, start from 1,..,Nt
    Burgers* bg;
    double c0, c1, dtdx, dtdx2;
    CostFunctorProc(double proc_error_variance, int t, Burgers* bg);
};


struct CostFunctorProc_fwd : CostFunctorProc {
    CostFunctorProc_fwd(double proc_error_variance, int t, Burgers* bg): CostFunctorProc(proc_error_variance, t, bg) {}
    template <typename T> bool operator()(const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const;
    static CostFunction* createAutoDiffCostFunction(double proc_error_variance, int t, Burgers* bg);
};

struct CostFunctorProc_frogleap : CostFunctorProc {
    CostFunctorProc_frogleap(double proc_error_variance, int t, Burgers* bg): CostFunctorProc(proc_error_variance, t, bg) {}
    //xp2,(t-2) time slice, just one location (same loc),  xp1: (t-1) time slice, , x: t time slice, i location
    template <typename T> bool operator()(const T* const xp2, const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const;
    static CostFunction* createAutoDiffCostFunction(double proc_error_variance, int t, Burgers* bg);
};

struct CostFunctorX0 {
    double init_error_variance, init_error_sqrtinv; 
    double init_priori;
    CostFunctorX0(double init_error_variance, double init_priori);
    template <typename T> bool operator()(const T* const x, T* residual) const;
    static CostFunction* createAutoDiffCostFunction(double init_error_variance, double init_priori);
};

