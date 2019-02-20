// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: keir@google.com (Keir Mierle)
//
// A simple example of using the Ceres minimizer.
//
// Minimize 0.5 (10 - x)^2 using jacobian matrix computed using
// automatic differentiation.

#include "ceres/ceres.h"
#include "glog/logging.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

// A templated cost functor that implements the residual r = 10 -
// x. The method operator() is templated so that we can then use an
// automatic differentiation wrapper around it to generate its
// derivatives.
// (x-5)^2 + (8-x)^2

//must be on a grid (iloc)
struct CostFunctorR {
    double obs, obs_error_variance, obs_error_sqrtinv;
    int iloc;
    CostFunctorR(int iloc, double obs, double obs_error_variance){
        this->iloc = iloc;
        this->obs = obs;
        this->obs_error_variance = obs_error_variance;
        this->obs_error_sqrtinv  = sqrt(1./obs_error_variance);
    }
    template <typename T> bool operator()(const T* const x, T* residual) const {
        residual[0] = obs_error_sqrtinv * (x[0] - obs);
        return true;
    }
};

//process
struct CostFunctorX {
    double proc_error_variance, proc_error_sqrtinv;
    int t; //0 means init condition, start from 1,..,Nt
    Burgers* bg;
    double c0, c1, dtdx, dtdx2;
    CostFunctorX(double proc_error_variance, int t, Burgers* bg){
        this->proc_error_variance = proc_error_variance;
        this->proc_error_sqrtinv = sqrt(1./proc_error_variance);
        this->t = t;
        this->bg = bg;
        bg->getConvenient(&dtdx, &dtdx2, &c0, &c1);
    }
};

struct CostFunctorX_fwd : CostFunctorX {
    CostFunctorX_fwd(double proc_error_variance, int t, Burgers* bg): CostFunctorX(proc_error_variance, t, bg) {}
    template <typename T> bool operator()(const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const {
        residual[0] = proc_error_sqrtinv*(xp1[0] - 0.5*dtdx*xp1[0]*(xp1p1[0]-xp1m1[0]) + 0.5*c1*(xp1m1[0]+xp1p1[0] - 2.*xp1[0]) - x[0]);
        return true;
    }
};

struct CostFunctorX_frogleap : CostFunctorX {
    CostFunctorX_frogleap(double proc_error_variance, int t, Burgers* bg): CostFunctorX(proc_error_variance, t, bg) {}
    //xp2,(t-2) time slice, just one location (same loc),  xp1: (t-1) time slice, , x: t time slice, i location
    template <typename T> bool operator()(const T* const xp2, const T* const xp1m1, const T* const xp1, const T* const xp1p1, const T* const x, T* residual) const {
        residual[0] = proc_error_sqrtinv*(c0*(xp2[0] + c1*(xp1p1[0] + xp1m1[0] - xp2[0]) - dtdx*(xp1[0]*(xp1p1[0]-xp1m1[0]))) - x[0]);
        return true;
    }
};

struct CostFunctorX0 {
    double init_error_variance, init_error_sqrtinv; 
    int nx;
    double* init_priori;
    CostFunctorX0(double init_error_variance, double* init_priori, int nx){
        this->init_error_variance = init_error_variance;
        this->init_error_sqrtinv = sqrt(1./init_error_sqrtinv);
        this->nx = nx;
        this->init_priori = new double[nx];
        for(int i=0; i<nx; i++)
            this->init_priori[i] = init_priori[i];
    }
    //x, all x at init condition [nx], residual: (x-xb)*sqrtinv, [nx]
    template <typename T> bool operator()(const T* const x, T* residual) const{
        for(int i=0; i<nx; i++){
            residual[i] = init_error_sqrtinv * (x[i] - init_priori[i]);
        }
        return true;
    }
};

struct CostFunctorX0_single {
    double init_error_variance, init_error_sqrtinv; 
    double init_priori;
    CostFunctorX0_single(double init_error_variance, double init_priori){
        this->init_error_variance = init_error_variance;
        this->init_error_sqrtinv = sqrt(1./init_error_variance);
        this->init_priori = init_priori;
    }
    //single point
    template <typename T> bool operator()(const T* const x, T* residual) const{
        residual[0] = init_error_sqrtinv * (x[0] - init_priori);
        return true;
    }
};



int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
    
    //read first guess Xs
    const int nx=1000;
    int nStep=288;
    int spinup_steps=200;
    double R=1e-6;
    double dx=40000.0; //2pi*Rearth ~ 40K, Km
    double dt=600;
    Burgers bg(nx, dx, R, dt, BC_OPTION_PERIOD, 0., LINEAR_OPTION_NONLINEAR, 0.); //note, dt<dx/U
    double** Xs_all=bg.readXs("output/burgers_1e-6.bin", nStep+spinup_steps, nx);
    double** Xs = &Xs_all[spinup_steps];

    //read observations
    Observations obss;
//    obss.load("obs_gen/exptB_1sttimefullobs/burgers_1e-6_noise_obs_sample.csv");
//    obss.load("obs_gen/exptA_modelobserror/burgers_1e-6_noise_obs_sample.csv");
    obss.load("obs_gen/exptC_spinup200/burgers_1e-6_noise_obs_sample.csv");
    std::cout<<"N obs read: "<<obss.nobs<<"\n";
    
    // Build the problem.
    Problem problem;

    // -1. Process Contraints
    for( int i=1; i<nStep+1; i++){
        for(int j=0; j<nx; j++){
            if(i==1 || bg.getNumericOption() == NUMERIC_OPTION_FORWARD){
                CostFunction* cost_functionX = 
                  new AutoDiffCostFunction<CostFunctorX_fwd, 1, 1, 1, 1, 1>(new CostFunctorX_fwd(0.05*0.05, i, &bg));
                problem.AddResidualBlock(cost_functionX, NULL, &Xs[i-1][(j-1+nx)%nx], &Xs[i-1][j], &Xs[i-1][(j+1)%nx], &Xs[i][j]);
            } else {
                CostFunction* cost_functionX = 
                  new AutoDiffCostFunction<CostFunctorX_frogleap, 1, 1, 1, 1, 1, 1>(new CostFunctorX_frogleap(0.05*0.05, i, &bg));
                problem.AddResidualBlock(cost_functionX, NULL, &Xs[i-2][j], &Xs[i-1][(j-1+nx)%nx], &Xs[i-1][j], &Xs[i-1][(j+1)%nx], &Xs[i][j]);
            }
        }
    }
    
    // -2. Observation Contraints
    for(int i=0; i< obss.nobs; i++){
        int ix=obss.obs[i].ix;
        int it=obss.obs[i].it;
        double obs_value=obss.obs[i].value;
        double obs_stdv =obss.obs[i].error;
        CostFunction* cost_functionR = 
          new AutoDiffCostFunction<CostFunctorR, 1, 1>(new CostFunctorR(ix, obs_value, obs_stdv*obs_stdv));
        problem.AddResidualBlock(cost_functionR, NULL, &Xs[it][ix]);
    }
    
    bool useInitPriori = false;
    if(useInitPriori){
//        CostFunction* cost_functionX0 = 
 //         new AutoDiffCostFunction<CostFunctorX0, nx, nx>(new CostFunctorX0(0.05*0.05, Xs[0], nx));
//        problem.AddResidualBlock(cost_functionX0, NULL, Xs[0]);
        for(int i=0; i<nx; i++){
            CostFunction* cost_functionX0 = 
              new AutoDiffCostFunction<CostFunctorX0_single, 1,1>(new CostFunctorX0_single(0.20*0.20, Xs[0][i]));
            problem.AddResidualBlock(cost_functionX0, NULL, &Xs[0][i]);
        }
    }
  

  // Run the solver!
  Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  options.max_num_iterations = 10;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
 // std::cout << "x : " << initial_x0[0] << "," << initial_x0[1]
 //           << " -> " << xs[0][0] << "," << xs[0][1] << "\n";
  // get the 0->6 position
//  for(int i=0; i<=6; i++){
//      printf("t=%d, position = %f, %f \n", i, xs[i][0], xs[i][1]);
//  }
  bg.writeXs("xa_burgers_1e-6.bin", Xs, nStep, nx);
  return 0;
}
