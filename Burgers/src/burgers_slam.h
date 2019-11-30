#ifndef burgers_slam_h
#define burgers_slam_h

#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <libconfig.h++>
#include <string>
#include "burgers_4dvar.h" // use CovModel, CostFunctorWb0

using namespace std;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::ResidualBlockId;
using ceres::Solve;

//for config read 
class Burgers_SLAM {
public:
    int window_steps=288;
    string xbs_file;
    int spinup_steps=0;
    string obs_file;
    string xas_file;
    string xb0_file="";
    bool xb0_constraint = false; // only works if option_hybrid_4dvar = false
    double xb0_err_stdv = 0.05;
    string xb0_err_stdv_file="";
    double proc_err_stdv = 0.05;
    bool xneighbor_constraint = false;
    bool xneighbor_only_x0 = false;  //only x0 not work well
    string xneighbor_constraint_xa_dxa = "xa";
    double xneighbor_diff_stdv = 0.02;
    int max_num_iterations = 10;

    // hybrid 4dvar slamda options:
    string option_hybrid_4dvar = ""; // "" for no hybrid, "weak" for "weak constraint" 
    double weak_hybrid_4dvar_wstd = 0.001; // weak hybrid 4dvar wstd (ew)
    string cov_dir = ""; // dir of B0_Modes.bin
    double B0_inflate_factor = 1.0;

    Burgers_SLAM() {}
    Burgers_SLAM(const char* config_file, const char* config_path);
    CovModel* readB0(int nx); 
    CovModel* inflate(CovModel* cov, double inflator);
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

// hybrid-4dvar-slamda

// weak constraint: || (x0 - xb0 - E0w0)/ew ||^2
struct CostFunctorX0W0_weak {
    double weak_stdv; //ew
    double* xb0; // will make a copy here
    CovModel* B0; // won't copy here, so caller should not change B0 hereafter
    CostFunctorX0W0_weak(double wstd, CovModel* B0, double* xb0);
    // paras: [x0, w0] 
    template <typename T> bool operator()(T const* const* paras, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(double wstd, CovModel* B0, double* xb0);
};

// CostFunctorProc_fwd_w0, for the 1st time-step, as w0 is needed
/*
struct CostFunctorProc_fwd_w0: CostFunctorProc_fwd {
    CovModel* B0;
    double* xb0;
    CostFunctorProc_fwd_w0(double proc_error_variance, Burgers* bg, CovModel* B0);
    // paras: [0]: x0, [1]: w0, [2]: lambda, residual: size of x0
    template <typename T> bool operator()(T const* const* paras, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(CovModel* B0, double* xb0);
}; */


/* lagrange can't work here, because ceres assumes cost function must be in ||..||^2, form
// new constraint term: lambda (x0-E0w0-xb), (to implement constraints: x0 = xb0 + E0w0)
struct CostFunctorLagrange {
    CovModel* B0;
    double* xb0; // will make a copy, so caller is free to change xb0 later
    CostFunctorLagrange(CovModel* B0, double* xb0);  
    // paras: [0]: x0, [1]: w0, [2]: lambda, residual: size of x0
    template <typename T> bool operator()(T const* const* paras, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(CovModel* B0, double* xb0);
};
*/

//x1 - (x0+x2)/2
struct CostFunctorXneighbor {
    double neighbor_diff_variance, neighbor_diff_sqrtinv;
    double xb1, xb0, xb2;
    CostFunctorXneighbor(double neighbor_diff_variance, double xb0=0.0, double xb1=0.0, double xb2=0.0);
    template <typename T> bool operator()(const T* const x1, const T* const x0, const T* const x2, T* residual) const;
    static CostFunction* createAutoDiffCostFunction(double neighbor_diff_variance, double xb0=0.0, double xb1=0.0, double xb2=0.0);
};


struct ResiBlockGroup{
    vector<string> group_name;
    vector<int> start_ids, end_ids;
    vector<double> each_costs; //all of residual blocks
    double total_cost;
    vector<double> group_costs; //costs sumed by group
    bool evaluated = false;

    void add_resi_block_group(string gname, int start_idx, int end_idx);  //add a group manually 

    void evaluate_all(Problem & problem, vector<ResidualBlockId> & resblock_ids); // do evaluation (calculate cost function value for each residual block) 

    void sum_metrics_bygroup(); //after evaluate_all, sum costs for every group (group_costs)

    void output_cost_bygroup(const char* csvname);  //output group_costs to a file

    void clear_groups(); //clear the groups

    void add_groups_byComponent(int cnt_resiblock_proc, int cnt_resiblock_obs, int cnt_resiblock_x0, int cnt_resiblock_xneighbor, int cnt_resiblock_hybrid_weak); //create automatic groups (by 4 components & 1 whole)

    void add_groups_proc_subgroup(int cnt_resiblock_proc, int subgroup_size, int offset=0); //create automatic groups (for proc subgroups, each subgroup has fixed size)

    //for the normal order of residual blocks: proc, obs, X0, here offset should be cnt_resiblock_proc
    void add_groups_obs_subgroup(int cnt_resiblock_obs, int subgroup_size, int offset=0); //create automatic groups (for obs subgroups, each subgroup has fixed size)

    void add_groups_xneigbor_subgroup(int cnt_resiblock_xneighbor, int subgroup_size, int offset=0);

};

struct ResiBlockGroup_Config{
    bool output_cost_groups = false;
    bool output_obs_subgroups = false;
    bool output_proc_subgroups = false;
    bool output_xneighbor_subgroups = false;
    int obs_subgroup_size = 1;
    int proc_subgroup_size = 1;
    int xneighbor_subgroup_size = 1;
    string xb_cost_filename; 
    string xa_cost_filename;

    ResiBlockGroup_Config(const char* config_file, const char* config_path);  //create groups from config file
};

#endif
