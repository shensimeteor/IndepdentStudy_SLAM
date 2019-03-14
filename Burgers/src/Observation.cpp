#include "Observation.h"
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector> 

SingleObs::SingleObs(int it, int ix, double value, double error): it(it), ix(ix), value(value), error(error){
}

void SingleObs::setAll(int it, int ix, double value, double error){
    this->it = it;
    this->ix = ix;
    this->value = value;
    this->error = error;
}

Observations::Observations(){
    this->nobs=-1;
    this->obs = NULL;
}

Observations::Observations(int nobs){
    this->nobs = nobs;
    this->obs = new SingleObs[nobs];
}
/*

Observations::Observations(int nobs, SingleObs* obs){
    this->nobs = nobs;
    this->obs = new SingleObs[nobs];
    for(int i=0; i<nobs; i++){
        this->obs[i] = obs[i];
    }
}
*/
Observations::~Observations(){
    delete this->obs;
}


void Observations::setByCopy(int nobs, SingleObs* obs){
    this->nobs = nobs;
    this->obs = new SingleObs[nobs];
    for(int i=0; i<nobs; i++){
        this->obs[i] = obs[i];
    }
}

void Observations::load(const char* obsfile){
    std::ifstream inFile(obsfile, std::ios::in);
    std::string lineStr;
    std::getline(inFile, lineStr);  //skip header
    int iobs=0;
    std::vector<std::string> v;
    while(std::getline(inFile, lineStr)){
        v.push_back(lineStr);
    }
    if(this->obs == NULL){ //test number of lines first
        this->nobs = v.size();
        this->obs = new SingleObs[this->nobs];
    }
    for(int i=0; i< v.size(); i++){
        std::stringstream ss(v[i]);
        std::string strT, strX, strObs, strErr;
        std::getline(ss, strT, ',');
        std::getline(ss, strX, ',');
        std::getline(ss, strObs, ',');
        std::getline(ss, strErr, ',');
        int t=atoi(strT.c_str());
        int x=atoi(strX.c_str());
        double val = atof(strObs.c_str());
        double err = atof(strErr.c_str());
        obs[iobs++].setAll(t, x, val, err);
    }
}
