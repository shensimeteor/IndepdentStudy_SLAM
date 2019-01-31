
class SingleObs{
public:
    int it, ix;
    double value, error;
    SingleObs(int it=-1, int ix=-1, double value=0., double error=9999.);
    void setAll(int it, int tx, double value, double error);
};


class Observations{
public:
    SingleObs** obs;
    int nt_obs, nx_obs;

    Observations(int nt_obs, int nx_obs);
    void load(const char* obsfile);
    ~Observations();

};

