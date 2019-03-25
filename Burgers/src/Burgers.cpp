#include "Burgers.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h++>

using namespace std;
using namespace libconfig; 

void Burgers::setConvenient(){
    dtdx = dt/dx;
    dtdx2 = dt/(dx*dx);
    c1 = (2./R) * dtdx2;
    c0 = 1. / (1.+c1);
}

void Burgers::getConvenient(double* dtdx, double *dtdx2, double *c0, double *c1){
    *dtdx = this->dtdx;
    *dtdx2 = this->dtdx2;
    *c0 = this->c0;
    *c1 = this->c1;
}

int Burgers::getNumericOption(){
    return numeric_option;
}

void Burgers::setupRandomGenerator(unsigned seed, double noise_mean, double noise_stdv){
    this->generator = new std::default_random_engine(seed);
    this->distribution = new std::normal_distribution<double>(noise_mean, noise_stdv);
    this->noise = new double[nx];
}

Burgers::Burgers(){
    nx=100;
    istep=0;
    dx=1;
    R=1000;
    dt=1;
    this->setBC(BC_OPTION_PERIOD,1);
    this->setLinearOption(LINEAR_OPTION_LINEAR, 10);
    curX = new double[nx];
    preX = new double[nx];
    preX2 = new double[nx];
    this->setConvenient();
}

Burgers::Burgers(int nx, double dx, double R, double dt, int bc_option, double bc_value, int linear_option, double linear_velocity, int numeric_option, int stochastic_option, double noise_mean, double noise_stdv, bool preset_seed, unsigned seed){
    this->nx=nx;
    this->dx=dx;
    this->R=R;
    this->dt=dt;
    this->istep=0;
    this->setBC(bc_option, bc_value);
    this->setLinearOption(linear_option, linear_velocity);
    curX = new double[nx];
    preX = new double[nx];
    preX2 = new double[nx];
    this->setConvenient();
    this->numeric_option = numeric_option;
    this->stochastic_option = stochastic_option;
    if(stochastic_option == STOCHASTIC_OPTION_GAUSSIAN){
        this->preset_seed = preset_seed;
        if(preset_seed){
            this->seed = seed;
        }else{
            this->seed = std::chrono::system_clock::now().time_since_epoch().count(); 
        }
        this->setupRandomGenerator(this->seed, noise_mean, noise_stdv);
    }
}

Burgers::Burgers(const char* config_file, const char* config_path){
    Config cfg; 
    try{
        cfg.readFile(config_file);
    }catch(const FileIOException &fioex){
        cerr << "I/O Error while reading file: "<<config_file<<endl;
        exit(1);
    }catch(const ParseException &pex){
        cerr << "Parse Error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << endl;
        exit(1);
    }

    const Setting& root = cfg.getRoot();

    try{
        Setting& burgers = cfg.lookup(config_path);
        int nx = burgers["nx"];
        double R = burgers["R"];
        double dx = burgers["dx"];
        double dt = burgers["dt"];
        int bc_option = burgers["bc_option"];
        double bc_fixed_value = 0.;
        if(bc_option == 1){
            bc_fixed_value = burgers["bc_fixed_value"];
        }
        int linear_option = burgers["linear_option"];
        double linear_velocity = 0.;
        if(linear_option == LINEAR_OPTION_LINEAR){
            linear_velocity = burgers["linear_velocity"];
        }
        int numeric_option = burgers["numeric_option"];
        int stochastic_option = burgers["stochastic_option"];
        double stochastic_proc_err_mean=0., stochastic_proc_err_stdv=0.;
        bool preset_seed = false;
        unsigned seed = 0;
        if(stochastic_option == 2){
            stochastic_proc_err_mean = burgers["stochastic_proc_err_mean"];
            stochastic_proc_err_stdv = burgers["stochastic_proc_err_stdv"];
            preset_seed = burgers["preset_seed"];
            if(preset_seed){
                long long lseed;
                burgers.lookupValue("seed", lseed);
                seed = (unsigned) lseed;
            }
        }
        this->nx = nx;
        new (this)Burgers(nx, dx, R, dt, bc_option, bc_fixed_value, linear_option, linear_velocity, numeric_option, stochastic_option, stochastic_proc_err_mean, stochastic_proc_err_stdv, preset_seed, seed);
    }catch(const SettingNotFoundException &nfex){
        cerr << "SettingNotFound at " << nfex.getPath() << endl;
        exit(1);
    }catch(const SettingTypeException &stex){
        cerr << "SettingTypeException at " << stex.getPath() << endl;
        exit(1);
    }
}


Burgers::~Burgers(){
    delete curX;
    delete preX;
    delete preX2;
}


void Burgers::addNoiseToCurX(){
    for(int i=0; i<nx; i++){
        curX[i] += (*this->distribution)((*(this->generator)));
    }
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

void Burgers::setStochasticOption(int stochastic_option, double noise_mean, double noise_stdv){
    if(stochastic_option == STOCHASTIC_OPTION_GAUSSIAN){
        this->stochastic_option = stochastic_option;
        this->noise_mean = noise_mean;
        this->noise_stdv = noise_stdv;
    }else if(stochastic_option == STOCHASTIC_OPTION_NONE){
        this->stochastic_option = stochastic_option;
    }else{
        printf("Error, stochastic_option unknown\n");
        exit(1);
    }
}

void Burgers::setLinearOption(int linear_option, double linear_velocity){
    if(linear_option == LINEAR_OPTION_LINEAR){
        this->linear_option = linear_option;
        this->linear_velocity = linear_velocity;
    }else if(linear_option == LINEAR_OPTION_NONLINEAR){
        this->linear_option = linear_option;
    }else{
        printf("Error, linear_option unknown\n");
        exit(1);
    }
}

void Burgers::init(double initX[]){
    istep=0;
    for(int i=0; i<nx; i++){
        curX[i] = initX[i];
        preX[i] = 0.0;
        preX2[i] = 0.0;
    }
}

void Burgers::advanceStep(){
    //copy curX -> preX
    istep++; 
    //forward step
    if(istep==1 || numeric_option == NUMERIC_OPTION_FORWARD){
        for(int i=0; i<nx; i++){
            preX[i] = curX[i];
        }
        if(bc_option == BC_OPTION_FIXED){
            curX[0] = bc_value;
            curX[nx-1] = bc_value;
            for(int i=1; i<nx-1; i++){
                double velo = linear_option==LINEAR_OPTION_LINEAR? linear_velocity : preX[i];
                curX[i] = preX[i] - 0.5*dtdx*velo*(preX[i+1] - preX[i-1]) + 0.5*c1*(preX[i+1]+preX[i-1] - 2*preX[i]);
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                double velo = linear_option==LINEAR_OPTION_LINEAR? linear_velocity : preX[i];
                curX[i] = preX[i] - 0.5*dtdx*velo*(preX[ip1] - preX[im1]) + 0.5*c1*(preX[ip1]+preX[im1] - 2*preX[i]);
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
                double velo = linear_option==LINEAR_OPTION_LINEAR? linear_velocity : preX[i];
                curX[i] = c0*(preX2[i] + c1*(preX[i+1] + preX[i-1] - preX2[i]) - dtdx*(velo*(preX[i+1] - preX[i-1])));
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                double velo = linear_option==LINEAR_OPTION_LINEAR? linear_velocity : preX[i];
                curX[i] = c0*(preX2[i] + c1*(preX[ip1] + preX[im1] - preX2[i]) - dtdx*(velo*(preX[ip1] - preX[im1])));
            }
        }
    }
    if(stochastic_option == STOCHASTIC_OPTION_GAUSSIAN){
        this->addNoiseToCurX();
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

void Burgers::advanceNStepsAndOutputBin(int N, const char* file, int output_step_t, int output_step_x){
    std::ofstream f(file, std::ios::out|std::ios::binary);
    if(f){
        int n_outX=(nx - 1)/output_step_x + 1;
        double *outputX = new double[n_outX];
        for(int i=0; i<n_outX; i++){
            outputX[i] = curX[i*output_step_x]; 
        }
        f.write((char*)outputX, n_outX*sizeof(double));

        for(int i=1; i<=N; i++){
            this->advanceStep();
            if( (i-1) % output_step_t == 0){
                for(int i=0; i<n_outX; i++){
                    outputX[i] = curX[i*output_step_x];
                }
                f.write((char*)outputX, n_outX*sizeof(double));
            }

        }
        f.close();
    }
}


//real Xs in file, [Nx+1][Nt]
double** Burgers::readXs(const char* file, int Nt, int Nx){
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if(!f){
        printf("file not found!\n");
        exit(1);
    }else{
        double** x = new double*[Nt+1];
        for(int i=0; i<=Nt; i++){
            x[i]=new double[Nx];
            f.read((char*)x[i], Nx*sizeof(double));
        }
        return x;
    }
}
        


void Burgers::writeXs(const char* file, double** Xs, int Nt, int Nx){
    std::ofstream f(file, std::ios::out|std::ios::binary);
    if(f){
        for(int i=0; i<=Nt; i++){
            f.write((char*)Xs[i], Nx*sizeof(double));
        }
        return;
    }
    printf("fail to create file %s\n", file);
}

void Burgers::outputConfig(const char* config_file){
    Config cfg;

    Setting& root= cfg.getRoot();
    Setting& model = root.add("model", Setting::TypeGroup);
    Setting& burgers = model.add("burgers", Setting::TypeGroup);
    burgers.add("nx", Setting::TypeInt) = this->nx;
    burgers.add("dx", Setting::TypeFloat) = this->dx;
    burgers.add("R", Setting::TypeFloat) = this->R;
    burgers.add("dt", Setting::TypeFloat) = this->dt;
    burgers.add("bc_option", Setting::TypeInt) = this->bc_option;
    if(this->bc_option == 1){
        burgers.add("bc_fixed_value", Setting::TypeFloat) = this->bc_value;
    }
    burgers.add("linear_option", Setting::TypeInt) = this->linear_option; 
    if(this->linear_option == 1) {
        burgers.add("linear_velocity", Setting::TypeFloat) = this->linear_velocity;
    }
    burgers.add("stochastic_option", Setting::TypeInt) = this->stochastic_option;
    if(this->stochastic_option == 2){
        burgers.add("preset_seed", Setting::TypeBoolean) = this->preset_seed;
        burgers.add("seed", Setting::TypeInt64) = (long long)this->seed;
        burgers.add("stochastic_proc_err_mean", Setting::TypeFloat) = this->noise_mean;
        burgers.add("stochastic_proc_err_stdv", Setting::TypeFloat) = this->noise_stdv;
    }
    try{
        cfg.writeFile(config_file);
    }catch(FileIOException &fio){
        cerr << "I/O Error when writng file " << config_file << endl;
        exit(1);
    }
}

int Burgers::getNx(){
    return nx;
}

int Burgers::getIstep(){
    return istep;
}

double* Burgers::getCurrentX(){
    return curX;
}
