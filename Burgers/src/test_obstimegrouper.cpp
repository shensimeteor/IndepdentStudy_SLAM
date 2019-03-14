#include "Observation.h"
#include "burgers_4dvar.h"
#include <iostream> 
using namespace std;

int main(){
    Observations obss;
    obss.load("test_obs.csv");
    ObsTimeGrouper obstg(288, 6, 0, 0, 0);
    obstg.group(obss);

    cout<<obstg.nt_obs<<endl;
    for(int i=0; i<obstg.nt_obs; i++){
            cout<<i<<" "<<obstg.obss[i].nobs<< 
            " "<<obstg.obss[i].obs[0].ix<<" "<<obstg.obss[i].obs[0].it<<" "<<
            obstg.obss[i].obs[obstg.obss[i].nobs-1].ix<<" "<<
            obstg.obss[i].obs[obstg.obss[i].nobs-1].it<<endl;
    }
}

