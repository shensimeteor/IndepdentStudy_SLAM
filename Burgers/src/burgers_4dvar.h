#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <libconfig.h++>
#include <vector>
#include <string>

using namespace std;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::ResidualBlockId;
using ceres::Solve;


struct CovModel{
    int nx = 1000;
    int n_mode = 1000;
    double ** modes; //[n_mode][nx]
    CovModel() { }
    CovModel(const char* cov_mode_file, int nx); 
    void load(const char* cov_mode_file, int nx); 
    void updateXb0byW(double* w, double* xb0);
};


//to group Observations into multiple Observations objects, based on time
class ObsTimeGrouper{
public:
    int steps_unit=6; //should be same with da_steps_unit, how many steps will be grouped 
    int start_step=0; //0
    int left_range_steps=0, right_range_steps=0; 
    int window_steps=288;
    int nt_obs; //number of timegroups of obs
    int last_step; //step of last obs, i.e. the last step needed to be integrated to
    Observations* obss; //size of group-size, i.e. each member is obs for each grouped-time
    int* obs_steps;  // same size with obss, obs_steps[idx] = model-step of that idx
    int n_obs; //total number of obs (~nt*nx)
    ObsTimeGrouper() {}
    ObsTimeGrouper(int window_steps, int steps_unit, int start_step=0, int left_range_steps=0, int right_range_steps=0);
    void group(Observations& allobs);//allobs must already be sorted by time
};


//for config read 
class Burgers_4DVar {
public:
    int window_steps=288;
    int da_steps_unit=6;
    int da_steps;
    bool weak_constraint = true;
    string xb0_file;
    string obs_file;
    string xas_file;
    string cov_dir; //B0 file, B0_modes.bin, Qt file: Q<t in XXX>_modes.bin
    double B0_inflate_factor=1.0;
    double Qt_inflate_factor=1.0;
    int max_num_iterations;
    Burgers_4DVar() {}
    Burgers_4DVar(const char* config_file, const char* config_path);
    CovModel* readB0(int nx); 
    CovModel* readQt(int nx, int t); //t in range [1,window_steps]
    CovModel* inflate(CovModel* cov, double inflator);
    CovModel* check_read_inflate_allQts(ObsTimeGrouper* obstg, int nx, int& nt_Qts); 
    // a high-level API, to check, read, inflate all Qts (return array of Qt, and its size nt_Qts), by time info from obstg
};


class SingleTimeObsOperator{
public:
    int nx_obs;
    int* obs_xidx;
    SingleTimeObsOperator() {nx_obs = -1;}
    SingleTimeObsOperator(const Observations& obss); //must make sure, obss only contains same time step, or canbe treated as a single time step (for da)
    void setByObservations(const Observations& obss); 
};


template <class T>
class Burgers_T: public Burgers{
    T *curX, *preX, *preX2; 
public: 
    Burgers_T() {}
    Burgers_T(const Burgers& bg); //copy constructor
    void init(T* initX); 
    void advanceStep();
    T* getCurrentX();
};


//for w^T * w pattern residual (used in 4DVar B term and Weak-Constraint 4Dvar's B&Q terms)
struct CostFunctorWb0 {
    int w_size;
    CostFunctorWb0(int w_size); 
    template <typename T> bool operator()(T const* const* ptr_w, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(int w_size); 
};

//for strong constraint 4DVar
struct CostFunctor_4DVar_FullObs{
    //Covariance term
    CovModel* B0; 
    //Obs & ObsOperator 
    ObsTimeGrouper* obstg; //grouped obs
    SingleTimeObsOperator* obsop; //obsop[nt_obs]
    //Model
    Burgers* bg; 
    //background xb0
    double* xb0;

    CostFunctor_4DVar_FullObs(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0); 
    template <typename T> bool operator()(T const* const* ptr_w, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0);
}; 

struct  CostFunctor_Weak4DVar_FullObs: public CostFunctor_4DVar_FullObs{
    int nt_Qts; //size of Qt array, should be nt_obs-1 
    CovModel* Qts;

    CostFunctor_Weak4DVar_FullObs(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0, CovModel* Qt);
    template <typename T> bool operator()(T const* const* ptr_w, T* residual) const;
    static CostFunction* createDynamicAutoDiffCostFunction(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0, CovModel* Qt);
};


