#include "Observation.h"

int main(){
    int nt_obs = 25;
    Observations obss(nt_obs);
    obss.load("obs_gen/burgers_1e-6_noise_obs_sample.csv");
}
