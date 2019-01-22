#define BC_OPTION_FIXED 1
#define BC_OPTION_PERIOD 2
#define LINEAR_OPTION_LINEAR 1
#define LINEAR_OPTION_NONLINEAR 2

class Burgers{
private:
    int nx, istep;
    double dx,R,dt;
    double *curX, *preX, *preX2;
    double dtdx, dtdx2, c1, c0;  //set for convenient
    int bc_option; //boundary condition
    double bc_value; //bc value (for fixed bc_option)
    void setConvenient();
    int linear_option; //linear or nonlinear option
    double linear_velocity;

public:
    Burgers();
    Burgers(int nx, double dx, double R, double dt, int bc_option=BC_OPTION_FIXED, double bc_value=0., int linear_option=LINEAR_OPTION_NONLINEAR, double linear_velocity=10.);
    ~Burgers();
//    void init();
    void init(double initX[]);
    void setBC(int bc_option, double bc_value);
    void setLinearOption(int linear_option, double linear_velocity);
    void advanceStep();
    void advanceStep(double x[]);
    double** advanceNSteps(int N);  // return X[N+1][nx]
    void outputBinX2d(const char* file, double **x, int N); //x should be [N+1][nx]
    void advanceNStepsAndOutputBin(int N, const char* file, int output_step_t=1, int output_step_x=1); //advance and output, so that no need to save every step in memory
};
