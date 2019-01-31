#include "Observation.h"
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

SingleObs::SingleObs(int it, int ix, double value, double error): it(it), ix(ix), value(value), error(error){
}

void SingleObs::setAll(int it, int ix, double value, double error){
    this->it = it;
    this->ix = ix;
    this->value = value;
    this->error = error;
}

Observations::Observations(int nt_obs, int nx_obs){
    this->nt_obs = nt_obs;
    this->nx_obs = nx_obs;
    this->obs = new SingleObs*[nt_obs];
    for(int i=0; i<nt_obs; i++){
        this->obs[i] = new SingleObs[nx_obs];
    }
}

Observations::~Observations(){
    for(int i=0; i<nt_obs; i++){
        delete this->obs[i];
    }
}

void Observations::load(const char* obsfile){
    std::ifstream inFile(obsfile, std::ios::in);
    std::string lineStr;
    std::getline(inFile, lineStr);  //skip header
    int prevT=-1, prevX=-1, it=-1, ix=-1;
    while(std::getline(inFile, lineStr)){
        std::stringstream ss(lineStr);
        std::string strT, strX, strObs, strErr;
        std::getline(ss, strT, ',');
        std::getline(ss, strX, ',');
        std::getline(ss, strObs, ',');
        std::getline(ss, strErr, ',');
        int t=atoi(strT.c_str());
        int x=atoi(strX.c_str());
        double val = atof(strObs.c_str());
        double err = atof(strErr.c_str());
        if(t != prevT)  { it++; ix=-1; }
        ix++;
        printf("%d, %d, %d, %d, %f, %f\n", it,ix,t,x,val,err);
        obs[it][ix].setAll(it, ix, val, err);
        prevT=t;
    }
}
