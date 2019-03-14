
#ifndef OBSERVATION_H
#define OBSERVATION_H

class SingleObs{
public:
    int it, ix;
    double value, error;
    SingleObs(int it=-1, int ix=-1, double value=0., double error=9999.);
    void setAll(int it, int tx, double value, double error);
};


class Observations{
public:
    SingleObs* obs;
    int nobs;

    Observations(int nobs);
    Observations();
    //Observations(int nobs, SingleObs* obs);
    void setByCopy(int nobs, SingleObs* obs);
    void load(const char* obsfile);
    ~Observations();

};

#endif
