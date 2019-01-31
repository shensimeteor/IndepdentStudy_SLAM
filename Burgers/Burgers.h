#define BC_OPTION_FIXED 1
#define BC_OPTION_PERIOD 2
#define LINEAR_OPTION_LINEAR 1
#define LINEAR_OPTION_NONLINEAR 2
#define NUMERIC_OPTION_FROGLEAP 1
#define NUMERIC_OPTION_FORWARD 2
#define STOCHASTIC_OPTION_NONE 1
#define STOCHASTIC_OPTION_GAUSSIAN 2
#include <random>
#include <chrono>


class Burgers{
private:
    int nx, istep;
    double dx,R,dt;
    double *curX, *preX, *preX2;
    double dtdx, dtdx2, c1, c0;  //set for convenient
    int bc_option; //boundary condition
    double bc_value; //bc value (for fixed bc_option)
    int linear_option; //linear or nonlinear option
    double linear_velocity;
    int numeric_option; //
    int stochastic_option; // not implemented yet
    double noise_mean, noise_stdv; 
    unsigned seed;
    std::default_random_engine *generator;
    std::normal_distribution<double> *distribution; 
    
    void setConvenient();
    double *noise; 
    void addNoiseToCurX();


public:
    Burgers();
    Burgers(int nx, double dx, double R, double dt, int bc_option=BC_OPTION_FIXED, double bc_value=0., int linear_option=LINEAR_OPTION_NONLINEAR, double linear_velocity=10., int numeric_option=NUMERIC_OPTION_FROGLEAP, int stochastic_option=STOCHASTIC_OPTION_NONE, double noise_mean = 0., double noise_stdv = 0.);
    ~Burgers();
//    void init();
    void init(double initX[]);
    void setBC(int bc_option, double bc_value);
    void setupRandomGenerator(unsigned seed, double noise_mean, double noise_stdv);
    void setLinearOption(int linear_option, double linear_velocity);
    void setStochasticOption(int stochastic_option, double noise_mean, double noise_stdv);
    void getConvenient(double* dtdx, double* dtdx2, double* c0, double* c1);
    int getNumericOption();
    void advanceStep();
    void advanceStep(double x[]);
    double** advanceNSteps(int N);  // return X[N+1][nx]
    void outputBinX2d(const char* file, double **x, int N); //x should be [N+1][nx]
    void advanceNStepsAndOutputBin(int N, const char* file, int output_step_t=1, int output_step_x=1); //advance and output, so that no need to save every step in memory

    static double** readXs(const char* file, int Nt, int Nx); //Xs[Nt+1][Nx]
};
