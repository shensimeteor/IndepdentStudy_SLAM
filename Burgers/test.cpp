#include "Burgers.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

int main(){
    int nx=1000;
    double dx=40000.0; //2pi*Rearth ~ 40K, Km
    const double pi=3.14159;
    Burgers bg(nx, dx, 1000., 500., BC_OPTION_PERIOD); //note, dt<dx/U
    //init
    double *x = new double[nx];
    double U=-20;
    for(int i=0; i<nx; i++){
        x[i] = U * sin(i*dx * 2*pi/ (nx*dx));
    }
    bg.init(x);
    //run
    int Nstep = 1000;
//    double **Xs= bg.advanceNSteps(Nstep);
//    bg.outputBinX2d("burgers_run.bin", Xs, Nstep);
//    bg.advanceNStepsAndOutputBin(Nstep, "burgers_run.bin");
    bg.advanceNStepsAndOutputBin(Nstep, "burgers_run_periodBC.bin");
    return 0;

}
