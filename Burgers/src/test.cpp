#include "Burgers.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

int main(){
    int nx=1000;
    double dx=40000.0; //2pi*Rearth ~ 40K, Km
    const double pi=3.14159;
    //init
    double *x = new double[nx];
    double U=-20;
    for(int i=0; i<nx; i++){
        x[i] = U * cos(i*dx * 2*pi/ (nx*dx));
    }
    double R=1e-6;
    double dt=600;
    //nonlinear, with period BC
    int Nstep = 1440;
    string file_base_name="output/burgers_1e-6"; 
    Burgers bg(nx, dx, R, dt, BC_OPTION_PERIOD, 0., LINEAR_OPTION_NONLINEAR, 0.); //note, dt<dx/U
    bg.init(x);
    bg.advanceNStepsAndOutputBin(Nstep, (file_base_name+".bin").c_str());
    
    //high res
    int ratio=5;
    double *x_hi = new double[nx*ratio];
    for(int i=0; i<nx*ratio; i++){
        x_hi[i] = U*cos(i*dx/ratio * 2*pi/(nx*dx));
    }
    Burgers bg1(nx*ratio, dx/ratio, R, dt/ratio, BC_OPTION_PERIOD, 0., LINEAR_OPTION_NONLINEAR, 0.); //note, dt<dx/U
    bg1.init(x_hi);
    bg1.advanceNStepsAndOutputBin(Nstep*5, (file_base_name+"_high.bin").c_str(), 5, 5);
    
   
    //lores , with forward 
    Burgers bg2(nx, dx, R, dt, BC_OPTION_PERIOD, 0., LINEAR_OPTION_NONLINEAR, 0., NUMERIC_OPTION_FORWARD); //note, dt<dx/U
    bg2.init(x);
    bg2.advanceNStepsAndOutputBin(Nstep, (file_base_name+"_fwd.bin").c_str());

    //lores, frog leap, with noise
    Burgers bg3(nx, dx, R, dt, BC_OPTION_PERIOD, 0., LINEAR_OPTION_NONLINEAR, 0., NUMERIC_OPTION_FORWARD, STOCHASTIC_OPTION_GAUSSIAN, 0.0, 0.05); //note, dt<dx/U
    bg3.init(x);
    bg3.advanceNStepsAndOutputBin(Nstep, (file_base_name+"_noise.bin").c_str());

//    //linear, with period BC
//    Burgers bg2(nx, dx, 1./1000., 200., BC_OPTION_PERIOD, 0., LINEAR_OPTION_LINEAR, -U); //note, dt<dx/U
//    bg2.init(x);
//    bg2.advanceNStepsAndOutputBin(Nstep, "output/burgers_run_periodBC_linear.bin");
//
//    //linear, with fixed BC
//    Burgers bg3(nx, dx, 1./1000., 200., BC_OPTION_FIXED, 0., LINEAR_OPTION_LINEAR, -U); //note, dt<dx/U
//    bg3.init(x);
//    bg3.advanceNStepsAndOutputBin(Nstep, "output/burgers_run_fixedBC_linear.bin");

    return 0;

}
