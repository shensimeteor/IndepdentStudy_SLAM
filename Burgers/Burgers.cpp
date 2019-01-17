#include "Burgers.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

void Burgers::setConvenient(){
    dtdx = dt/dx;
    dtdx2 = dt/(dx*dx);
    c1 = (2./R) * dtdx2;
    c0 = 1. / (1.+c1);
}

Burgers::Burgers(){
    nx=100;
    istep=0;
    dx=1;
    R=1000;
    dt=1;
    this->setBC(0,1);
    curX = new double[nx];
    preX = new double[nx];
    preX2 = new double[nx];
    this->setConvenient();
}

Burgers::Burgers(int nx, double dx, double R, double dt, int bc_option, double bc_value){
    this->nx=nx;
    this->dx=dx;
    this->R=R;
    this->dt=dt;
    this->istep=0;
    this->setBC(bc_option, bc_value);
    curX = new double[nx];
    preX = new double[nx];
    preX2 = new double[nx];
    this->setConvenient();
}


Burgers::~Burgers(){
    delete curX;
    delete preX;
    delete preX2;
}

void Burgers::setBC(int bc_option, double bc_value){
    if(bc_option == BC_OPTION_FIXED){
        this->bc_option = bc_option;
        this->bc_value = bc_value;
    }else if(bc_option == BC_OPTION_PERIOD){
        this->bc_option = bc_option;
    }else{
        printf("Error, bc_option unknown\n");
        exit(1);
    }
}

void Burgers::init(double initX[]){
    istep=0;
    for(int i=0; i<nx; i++){
        curX[i] = initX[i];
    }
}

//for constant boundary condition =0 
void Burgers::advanceStep(){
    //copy curX -> preX
    istep++; 
    //forward step
    if(istep==1){
        for(int i=0; i<nx; i++){
            preX[i] = curX[i];
        }
        if(bc_option == BC_OPTION_FIXED){
            curX[0] = bc_value;
            curX[nx-1] = bc_value;
            for(int i=1; i<nx-1; i++){
                curX[i] = preX[i] - 0.5*dtdx*preX[i]*(preX[i+1] - preX[i-1]) + 0.5*c1*(preX[i+1]+preX[i-1] - 2*preX[i]);
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                curX[i] = preX[i] - 0.5*dtdx*preX[i]*(preX[ip1] - preX[im1]) + 0.5*c1*(preX[ip1]+preX[im1] - 2*preX[i]);
            }
        }
    }else{ //frog leap 
        for(int i=0; i<nx; i++){
            preX2[i] = preX[i];
            preX[i] = curX[i];
        }
        if(bc_option == BC_OPTION_FIXED){
            curX[0] = bc_value;
            curX[nx-1] = bc_value;
            for(int i=1; i<nx-1; i++){
                curX[i] = c0*(preX2[i] + c1*(preX[i+1] + preX[i-1] - preX2[i]) - dtdx*(preX[i]*(preX[i+1] - preX[i-1])));
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                curX[i] = c0*(preX2[i] + c1*(preX[ip1] + preX[im1] - preX2[i]) - dtdx*(preX[i]*(preX[ip1] - preX[im1])));
            }
        }
    }
}

void Burgers::advanceStep(double x[]){
    this->advanceStep();
    for(int i=0; i<nx; i++){
        x[i] = curX[i];
    }
}

double** Burgers::advanceNSteps(int N){
    //allocate
    double** x = new double*[N+1];
    for(int i=0; i<=N; i++)
        x[i]=new double[nx];
    //copy and run
    for(int i=0; i<nx; i++){
        x[0][i] = curX[i];
    }
    for(int t=0; t<N; t++){
        this->advanceStep();
        for(int i=0; i<nx; i++){
            x[t+1][i] = curX[i];
        }
    }
    return x;
}

void Burgers::outputBinX2d(const char* file, double **x, int N){
    std::ofstream f(file, std::ios::out|std::ios::binary);
    if(f){
        for(int i=0; i<=N; i++)
            f.write((char*)x[i], nx*sizeof(double));
        f.close();
    }

}

void Burgers::advanceNStepsAndOutputBin(int N, const char* file){
    std::ofstream f(file, std::ios::out|std::ios::binary);
    if(f){
        f.write((char*)curX, nx*sizeof(double));
        for(int i=1; i<=N; i++){
            this->advanceStep();
            f.write((char*)curX, nx*sizeof(double));
        }
        f.close();
    }


}

