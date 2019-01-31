#include "Observation.h"

int main(){
    int nt_obs = 25, nx_obs = 10;
    Observations obss(nt_obs, nx_obs);
    obss.load("obs_gen/burgers_1e-6_noise_obs_sample.csv");
}
