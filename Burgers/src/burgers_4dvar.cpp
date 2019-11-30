#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include <stdio.h>
#include "Burgers.h"
#include "Observation.h"
#include "math.h"
#include <fstream>
#include <libconfig.h++>
#include "burgers_4dvar.h"
#include <iostream>

using ceres::Jet;
using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::ResidualBlockId; 
using namespace std;
using namespace libconfig;

void CovModel::load(const char* cov_mode_file, int nx){
    std::ifstream f(cov_mode_file, std::ios::in | std::ios::binary);
    if(!f){
        printf("file not found -- %s!\n", cov_mode_file);
        exit(1);
    }else{
        /*
        this->nx = nx;
        f.seekg(0, f.end);
        int length = f.tellg();
        f.seekg(0, f.beg);
        this->n_mode = (length/(8*nx));
        this->modes = new double*[n_mode];
        printf("Cov Modes number = %d\n", this->n_mode);
        for(int i=0; i<n_mode; i++){
            this->modes[i]=new double[nx];
            f.read((char*)modes[i], nx*sizeof(double));
        }
        f.close(); */
        this->nx = nx;
        f.seekg(0, f.beg);
        std::vector<double*> vec; 
        while(f.peek()!=EOF){
            double* x = new double[nx];
            f.read((char*)x, nx*sizeof(double));
            vec.push_back(x);
        }
        //printf("vector_size = %d\n", vec.size());
        this->n_mode = vec.size();
        this->modes = new double*[n_mode];
        for(int i=0; i<this->n_mode; i++){
            this->modes[i] = vec[i];
        }
    }
}

//cov_mode_file, matrix of [nx][n_mode]
CovModel::CovModel(const char* cov_mode_file, int nx){
    this->load(cov_mode_file, nx);
}

void CovModel::inflate(double inflator) {
    double sqrt_inflator = sqrt(inflator);
    for(int i=0; i<n_mode; i++){
        for(int j=0; j<nx; j++){
            modes[i][j] *= sqrt_inflator;
        }
    }
}


// xb0 <- xb0 + Modes*w
void CovModel::updateXb0byW(double* w, double* xb0){
    int nw = this->n_mode;
    for(int i=0; i<nx; i++){
        for(int j=0; j<nw; j++){
            xb0[i]+= w[j] * modes[j][i];
        }
    }
}

//Burgers_SLAM
Burgers_4DVar::Burgers_4DVar(const char* config_file, const char* config_path){
    Config cfg;
    try{
        cfg.readFile(config_file);
        Setting & root = cfg.getRoot();
        Setting& varda = root.lookup(config_path);
        if(varda.exists("window_steps"))
            this->window_steps = varda["window_steps"];
        if(varda.exists("da_steps_unit"))
            this->da_steps_unit = varda["da_steps_unit"];
        this->da_steps = this->window_steps / this->da_steps_unit;
        if(varda.exists("weak_constraint"))
            this->weak_constraint = varda["weak_constraint"];
        if(varda.exists("B0_inflate_factor"))
            this->B0_inflate_factor = varda["B0_inflate_factor"];
        if(varda.exists("Qt_inflate_factor"))
            this->Qt_inflate_factor = varda["Qt_inflate_factor"];
        if(varda.exists("max_num_iterations"))
            this->max_num_iterations = varda["max_num_iterations"];
        varda.lookupValue("xb0_file", this->xb0_file);
        varda.lookupValue("obs_file", this->obs_file);
        varda.lookupValue("xas_file", this->xas_file);
        varda.lookupValue("cov_dir", this->cov_dir);
    }catch(ParseException &pex){
        cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                      << " - " << pex.getError() << endl;
        exit(1);
    }catch(FileIOException &fioex){
        cerr << "File I/O Error " << config_file << endl;
        exit(1);
    }catch(SettingNotFoundException &nfex){
        cerr << "Setting Not Found "<<nfex.getPath() << endl;
        exit(1);
    }
}

CovModel* Burgers_4DVar::readB0(int nx){
    string B0file = this->cov_dir + string("/B0_modes.bin");
    CovModel* pcov = new CovModel(B0file.c_str(), nx);
    return pcov;
}

CovModel* Burgers_4DVar::readQt(int nx, int t){
    char buffer[20];
    sprintf(buffer, "/Q%03d_modes.bin", t);
    string Qtfile = this->cov_dir + string(buffer);
    CovModel* pcov = new CovModel(Qtfile.c_str(), nx);
    return pcov;
}

CovModel* Burgers_4DVar::inflate(CovModel* B0, double inflator){
    B0->inflate(inflator);
    return B0;
}

CovModel* Burgers_4DVar::check_read_inflate_allQts(ObsTimeGrouper* obstg, int nx, int& nt_Qts){
    //get nt_Qts from obstg
    nt_Qts=0;
    if(obstg->start_step == 0){
        nt_Qts = obstg->nt_obs - 1;
    }else{
        nt_Qts = obstg->nt_obs;
    }
    CovModel* ptr_qts = new CovModel[nt_Qts];
    //read & inflate Qts files
    for(int i=0; i<nt_Qts; i++){
        int qt_tidx;
        if(obstg->start_step == 0){
            qt_tidx = obstg->obs_steps[i+1];
        }else{
            qt_tidx = obstg->obs_steps[i];
        }
        char buffer[20];
        sprintf(buffer, "/Q%03d_modes.bin", qt_tidx);
        string Qtfile = this->cov_dir + string(buffer);
        ptr_qts[i].load(Qtfile.c_str(), nx);
    }
    return ptr_qts;    
}


ObsTimeGrouper::ObsTimeGrouper(int window_steps, int steps_unit, int start_step, int left_range_steps, int right_range_steps){
    this->window_steps = window_steps;
    this->steps_unit = steps_unit;
    this->start_step = start_step;
    this->left_range_steps = left_range_steps; 
    this->right_range_steps = right_range_steps; 
    this->nt_obs = (window_steps - right_range_steps - start_step)/steps_unit + 1;
    this->last_step = start_step + (this->nt_obs-1) * steps_unit;
    this->obss = new Observations[this->nt_obs];
    this->obs_steps = new int[this->nt_obs];
    this->obs_steps[0] = this->start_step;
    for(int i=1; i<this->nt_obs; i++){
        this->obs_steps[i] = this->obs_steps[i-1] + steps_unit;
    }
}

void ObsTimeGrouper::group(Observations& allobs){
    int start_idx=-1, end_idx=0, igroup=0;
    int istep=start_step;
    int i=0;
    bool flag=false;
    //cout<<allobs.nobs<<endl;
    while(i<=allobs.nobs){
        if(i==allobs.nobs){
            if(flag){
                end_idx = i-1;
               // cout<<start_idx<<" "<<end_idx<<" "<<igroup<<endl;
                obss[igroup++].setByCopy(end_idx-start_idx+1, &allobs.obs[start_idx]);
            }
            i++;
        }
        else if( allobs.obs[i].it <= istep + right_range_steps && allobs.obs[i].it >= istep - left_range_steps){
            if(!flag)  {start_idx = i; flag=true;}
            i++;
        }else if( allobs.obs[i].it > istep + right_range_steps && flag){
            end_idx = i-1;
            //cout<<start_idx<<" "<<end_idx<<" "<<igroup<<endl;
            obss[igroup++].setByCopy(end_idx-start_idx+1, &allobs.obs[start_idx]);
            if(igroup==nt_obs) break;
            istep+= steps_unit; 
            flag=false;
        }
        else{
            i++;
        }
    }
    nt_obs = igroup;
    //get n_obs;
    n_obs = 0;
    for(int t=0; t<nt_obs; t++){
        n_obs += obss[t].nobs;
    }
}


SingleTimeObsOperator::SingleTimeObsOperator(const Observations& obss){
    this->nx_obs = obss.nobs;
    this->obs_xidx = new int[nx_obs];
    for(int i=0; i<nx_obs; i++){
        obs_xidx[i] = obss.obs[i].ix;
    }
}

void SingleTimeObsOperator::setByObservations(const Observations& obss){
    this->nx_obs = obss.nobs;
    this->obs_xidx = new int[nx_obs];
    for(int i=0; i<nx_obs; i++){
        obs_xidx[i] = obss.obs[i].ix;
    }
}
    

//Burgers_T
template <class T> Burgers_T<T>::Burgers_T(const Burgers& bg): Burgers(bg) {
    curX=new T[nx];
    preX=new T[nx];
    preX2=new T[nx];
}

template<class T> void Burgers_T<T>::init(T* initX){
    istep=0;
    for(int i=0; i<nx; i++){
        curX[i] = initX[i];
        preX[i] = (T) 0.;
        preX2[i] = (T) 0.;

    }
}
template <class T> T* Burgers_T<T>::getCurrentX(){
    return curX;
}

template <class T> void Burgers_T<T>::advanceStep() {
    //copy curX -> preX
    istep++; 
    //forward step
    if(istep==1 || numeric_option == NUMERIC_OPTION_FORWARD){
        for(int i=0; i<nx; i++){
            preX[i] = curX[i];
        }
        if(bc_option == BC_OPTION_FIXED){
            curX[0] = (T) bc_value;
            curX[nx-1] = (T) bc_value;
            for(int i=1; i<nx-1; i++){
                T velo = linear_option==LINEAR_OPTION_LINEAR? (T)linear_velocity : preX[i];
                curX[i] = preX[i] - 0.5*dtdx*velo*(preX[i+1] - preX[i-1]) + 0.5*c1*(preX[i+1]+preX[i-1] - 2.0*preX[i]);
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                T velo = linear_option==LINEAR_OPTION_LINEAR? (T)linear_velocity : preX[i];
                curX[i] = preX[i] - 0.5*dtdx*velo*(preX[ip1] - preX[im1]) + 0.5*c1*(preX[ip1]+preX[im1] - 2.0*preX[i]);
            }
        }
    }else{ //frog leap 
        for(int i=0; i<nx; i++){
            preX2[i] = preX[i];
            preX[i] = curX[i];
        }
        if(bc_option == BC_OPTION_FIXED){
            curX[0] = (T) bc_value;
            curX[nx-1] = (T) bc_value;
            for(int i=1; i<nx-1; i++){
                T velo = linear_option==LINEAR_OPTION_LINEAR? (T)linear_velocity : preX[i];
                curX[i] = c0*(preX2[i] + c1*(preX[i+1] + preX[i-1] - preX2[i]) - dtdx*(velo*(preX[i+1] - preX[i-1])));
            }
        }else if(bc_option == BC_OPTION_PERIOD){
            for(int i=0; i<nx; i++){
                int im1 = (i-1+nx) % nx;
                int ip1 = (i+1) %nx;
                T velo = linear_option==LINEAR_OPTION_LINEAR? (T)linear_velocity : preX[i];
                curX[i] = c0*(preX2[i] + c1*(preX[ip1] + preX[im1] - preX2[i]) - dtdx*(velo*(preX[ip1] - preX[im1])));
            }
        }
    }
    //if(stochastic_option == STOCHASTIC_OPTION_GAUSSIAN){
    //    this->addNoiseToCurX();
   // }
}





//CostFunctorWb0
CostFunctorWb0::CostFunctorWb0(int w_size): w_size(w_size) {}

template <typename T> bool CostFunctorWb0::operator()(T const* const* ptr_w, T* residual) const{
    for(int i=0; i<this->w_size; i++){
        residual[i] = (*ptr_w)[i];
    }
    return true;
}

CostFunction* CostFunctorWb0::createDynamicAutoDiffCostFunction(int w_size){
    CostFunctorWb0* functor = new CostFunctorWb0(w_size);
    DynamicAutoDiffCostFunction<CostFunctorWb0, 4>* cost_function = new 
        DynamicAutoDiffCostFunction<CostFunctorWb0, 4>(functor);
    cost_function->AddParameterBlock(w_size);
    cost_function->SetNumResiduals(w_size);
    return cost_function;
}


//CostFunctor_4DVar_FullObs
CostFunctor_4DVar_FullObs::CostFunctor_4DVar_FullObs(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double *xb0){
    //Covariance
    this->B0 = B0;
    //Obs & ObsOperator
    this->obstg = obstg;
    int nt_obs = obstg->nt_obs;
    this->obsop = new SingleTimeObsOperator[nt_obs];
    for(int i=0; i<nt_obs; i++){
        this->obsop[i].setByObservations(obstg->obss[i]);
    }
    //Model
    this->bg = bg;
    //background
    int nx = bg->getNx();
    this->xb0 = new double[nx];
    for(int i=0; i<nx; i++){
        this->xb0[i] = xb0[i];
    }
}

//w: input size, residual: n_obs
//template <typename T> bool CostFunctor_4DVar_FullObs::operator()(const T* const w, T* residual) const{
template <typename T> bool CostFunctor_4DVar_FullObs::operator()(T const* const* ptr_w, T* residual) const{
    Burgers_T<T>* bgt = new Burgers_T<T>(*bg);
    T* x0 = new T[this->bg->getNx()];
    // x0 = xb0 + E*w
    for(int i=0; i<this->bg->getNx(); i++){
        T sum=(T)0.0;
        for(int j=0; j<B0->n_mode; j++){
            sum = sum + (*ptr_w)[j] * B0->modes[j][i];
        }
        x0[i] = (T) (xb0[i]) + sum;

    }
    // M(x0)
    int last_step = obstg->last_step;
    int istep=0;
    int obs_tidx=0;
    int iobs=0;
    bgt->init(x0);
    while(istep<=last_step){
    //    cout<<"istep="<<istep<<endl;
        if(istep == obstg->obs_steps[obs_tidx]){
            T* curX = bgt->getCurrentX();
            SingleTimeObsOperator& st_obsop = obsop[obs_tidx];
            for(int ix=0; ix<st_obsop.nx_obs; ix++){
                int x = st_obsop.obs_xidx[ix];
      //          cout<<"ix="<<ix<<", x="<<x<<"iobs="<<iobs<<endl;
                residual[iobs++] = (curX[x] - obstg->obss[obs_tidx].obs[ix].value) / sqrt(obstg->obss[obs_tidx].obs[ix].error);
            }
            obs_tidx++;
        }
        bgt->advanceStep();
        istep = bgt->getIstep();
    }
    return true;
}

CostFunction* CostFunctor_4DVar_FullObs::createDynamicAutoDiffCostFunction(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0){
    int w_size = B0->n_mode;
    int n_obs = obstg->n_obs;
    CostFunctor_4DVar_FullObs* functor = new CostFunctor_4DVar_FullObs(B0, obstg, bg, xb0);
    DynamicAutoDiffCostFunction<CostFunctor_4DVar_FullObs, 4>* cost_function = new 
        DynamicAutoDiffCostFunction<CostFunctor_4DVar_FullObs, 4>(functor);
    cost_function->AddParameterBlock(w_size);
    cost_function->SetNumResiduals(n_obs);
    return cost_function;
}


// WeakConstraint 4DVar
CostFunctor_Weak4DVar_FullObs::CostFunctor_Weak4DVar_FullObs(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double *xb0, CovModel* Qts) :
    CostFunctor_4DVar_FullObs(B0, obstg, bg, xb0){
        this->Qts = Qts;
        if(obstg->start_step == 0){
            this->nt_Qts = obstg->nt_obs-1;
        }else{
            this->nt_Qts = obstg->nt_obs;
        }
}

template <typename T> bool CostFunctor_Weak4DVar_FullObs::operator()(T const* const* ws, T* residual) const{
    Burgers_T<T>* bgt = new Burgers_T<T>(*bg);
    T* x0 = new T[this->bg->getNx()];
    // x0 = xb0 + E0*w0
    for(int i=0; i<this->bg->getNx(); i++){
        T sum=(T)0.0;
        for(int j=0; j<B0->n_mode; j++){
            sum = sum + ws[0][j] * B0->modes[j][i];
        }
        x0[i] = (T) (xb0[i]) + sum;

    }
    int obs_tidx=0;
    int iobs=0;
    //if 0th step has obs
    if(obstg->obs_steps[0] == 0){
        SingleTimeObsOperator& st_obsop = obsop[0];
        for(int ix=0; ix<st_obsop.nx_obs; ix++){
            int x = st_obsop.obs_xidx[ix];
            residual[iobs++] = (x0[x] - obstg->obss[0].obs[ix].value) / sqrt(obstg->obss[0].obs[ix].error);
        }
        obs_tidx++;
    }
    //
    int curr_step=0;
    for(int irun = 0; irun< nt_Qts; irun++ ){
        int nsteps = obstg->obs_steps[obs_tidx] - curr_step;  //steps to integrate
        // x(t) = M(x(t-1)) , in curX
        bgt->init(x0);
        while(bgt->getIstep() < nsteps){
            bgt->advanceStep();
        }
        T* curX = bgt->getCurrentX();
        // curX += Et*wt
        for(int i=0; i<bgt->getNx(); i++){
            T sum=(T)0.0;
            for(int j=0; j<Qts[irun].n_mode; j++){
                sum = sum + ws[irun+1][j] * Qts[irun].modes[j][i];
            }
            curX[i] = curX[i] + sum;
        } 
        // residual = H(curX) - obs
        SingleTimeObsOperator& st_obsop = obsop[obs_tidx];
        for(int ix=0; ix<st_obsop.nx_obs; ix++){
            int x = st_obsop.obs_xidx[ix];
            residual[iobs++] = (curX[x] - obstg->obss[obs_tidx].obs[ix].value) / sqrt(obstg->obss[obs_tidx].obs[ix].error);
        }
        // copy curX to x0
        for(int i=0; i<bgt->getNx(); i++){
            x0[i] = curX[i];
        }
        // update obs_tidx, curr_step
        curr_step = obstg->obs_steps[obs_tidx];
        obs_tidx++;
    }
    //evaluate
    T rmse_diff_obs=(T)0.0;
    for(int i=0; i<iobs; i++){
        rmse_diff_obs +=  residual[i]*residual[i];
    }
    T t_iobs = (T) (iobs*1.0);
    rmse_diff_obs = sqrt(rmse_diff_obs/t_iobs);
    cout<<"iobs="<<iobs<<", rms_residual = "<<rmse_diff_obs<<endl;
    //for(int i=0; i<iobs; i++){
    //    cout<<"i="<<i<<", residual="<<residual[i]<<endl;
   // }
//    cout<<"curr_step="<<curr_step<<",  obs_tidx="<<obs_tidx<<endl;
    return true;
}


CostFunction* CostFunctor_Weak4DVar_FullObs::createDynamicAutoDiffCostFunction(CovModel* B0, ObsTimeGrouper* obstg, Burgers* bg, double* xb0, CovModel* Qts){
    int w_size = B0->n_mode;
    int n_obs = obstg->n_obs;
    CostFunctor_Weak4DVar_FullObs* functor = new CostFunctor_Weak4DVar_FullObs(B0, obstg, bg, xb0, Qts);
    DynamicAutoDiffCostFunction<CostFunctor_Weak4DVar_FullObs, 8>* cost_function = new 
        DynamicAutoDiffCostFunction<CostFunctor_Weak4DVar_FullObs, 8>(functor);
    //parameter: w0
    cost_function->AddParameterBlock(B0->n_mode);
//    cout<<B0->n_mode<<endl;
    //parameter: wt (Qt): t=1,..,nt_Qts
    for(int i=0; i<functor->nt_Qts; i++){ 
  //      cout<<Qts[i].n_mode<<endl;
        cost_function->AddParameterBlock(Qts[i].n_mode);
    }
    //residual: n_obs
    cost_function->SetNumResiduals(n_obs);
    return cost_function;
}




