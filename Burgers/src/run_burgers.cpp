#include "Burgers.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <libconfig.h++>
using namespace std;
using namespace libconfig;

struct ModelRun {
    int init_condition_pattern = 1; 
    double init_condition_v0 = -20.0;
    int n_steps = 1440;
    int output_t_ratio=1, output_x_ratio=1;
    string output_filename;
    string init_condition_file;

    ModelRun() {}

    ModelRun(const char* config_file, const char* config_path){
        Config cfg;
        try{
            cfg.readFile(config_file);
            Setting & root = cfg.getRoot();
            Setting& modelrun = root.lookup(config_path);
            if(modelrun.exists("init_condition_pattern"))
                this->init_condition_pattern = modelrun["init_condition_pattern"];
            if(modelrun.exists("init_condition_v0")){
                if(modelrun["init_condition_v0"].getType() == Setting::Type::TypeInt){
                    int tmp = modelrun["init_condition_v0"];
                    this->init_condition_v0 = (double)tmp;
                }else{
                    this->init_condition_v0 = modelrun["init_condition_v0"] ;
                }
            }
            if(modelrun.exists("n_steps"))
                this->n_steps = modelrun["n_steps"];
            if(modelrun.exists("output_t_ratio"))
                this->output_t_ratio = modelrun["output_t_ratio"];
            if(modelrun.exists("output_x_ratio"))
                this->output_x_ratio = modelrun["output_x_ratio"];
            modelrun.lookupValue("output_filename", this->output_filename);
            if(this->init_condition_pattern == 0){ //read from file
                modelrun.lookupValue("init_condition_file", this->init_condition_file);
            }
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
    
    void readX(const char* filename, double *x, int nx){
        ifstream f(filename, ios::in | ios::binary);
        if(!f){
            cout<< "file nout found!\n";
            exit(1);
        }else{
            f.read((char*) x, nx*sizeof(double));
        }
    }

    string getIC(double *x, int nx){
        const double pi = 3.14159;
        if(init_condition_pattern == 1){
            for(int i=0; i<nx; i++){
                x[i] = init_condition_v0 * cos(i * 2*pi/ nx);
            }
            return string("cosine 1-wave, v0=") + to_string(init_condition_v0); 
        }else if(init_condition_pattern == 0){
            readX(this->init_condition_file.c_str(), x, nx);
            return string("read from file") + this->init_condition_file;
        }
    }



};

int main(){
    Burgers bg("burgers_conf.in", "model.burgers");
    cout<< "[log] Read config from burgers_conf.in" << endl;
    bg.outputConfig("burgers_conf.out");
    cout<< "[log] Output config to burgers_conf.out" << endl;
    
    ModelRun modelrun("burgers_conf.in", "model_run");
    cout<<"[log] output_filename = " << modelrun.output_filename <<endl;
    
    double *x = new double[bg.getNx()];
    string ic_name = modelrun.getIC(x, bg.getNx());
    cout<<"[log] set InitCondition: " << ic_name << endl;
    //run model
    bg.init(x);
    bg.advanceNStepsAndOutputBin(modelrun.n_steps, modelrun.output_filename.c_str(),
            modelrun.output_t_ratio, modelrun.output_x_ratio);
    cout<<"[log] model finished! " << endl;
    return 0;
}
