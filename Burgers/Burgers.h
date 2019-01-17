#define BC_OPTION_FIXED 1
#define BC_OPTION_PERIOD 2

class Burgers{
private:
    int nx, istep;
    double dx,R,dt;
    double *curX, *preX, *preX2;
    double dtdx, dtdx2, c1, c0;  //set for convenient
    int bc_option; //boundary condition
    double bc_value; //bc value (for fixed bc_option)
    void setConvenient();

public:
    Burgers();
    Burgers(int nx, double dx, double R, double dt, int bc_option=BC_OPTION_FIXED, double bc_value=0.);
    ~Burgers();
//    void init();
    void init(double initX[]);
    void setBC(int bc_option, double bc_value);
    void advanceStep();
    void advanceStep(double x[]);
    double** advanceNSteps(int N);  // return X[N+1][nx]
    void outputBinX2d(const char* file, double **x, int N); //x should be [N+1][nx]
    void advanceNStepsAndOutputBin(int N, const char* file); //advance and output, so that no need to save every step in memory
};
