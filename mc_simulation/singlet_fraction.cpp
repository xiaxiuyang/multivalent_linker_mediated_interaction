#include "system.h"

using namespace std;

double getAverage(vector<double>& sequence);

int main(int argc, char* argv[]){

    omp_set_num_threads(8);// use 8 cores
    ofstream singlet,Naa,Nab,Nbb,Nba;
    singlet.open("singlet_fact.txt");
    Naa.open("number_aa_perC.txt");
    Nab.open("number_ab_perC.txt");
    Nbb.open("number_bb_perC.txt");
    Nba.open("number_ba_perC.txt");// to store the numbers of linkers in different states
    ofstream position;
    position.open("configuration.xyz");
    vector<double> singlet_tot,Naa_tot,Nab_tot,Nbb_tot,Nba_tot;
    vector<double> strands={1,2,4,6}; // number of strands
    double temperature=1.;
    double scaleL=stod(argv[1]);// size of the system
    double chemical_potential=stod(argv[2]);// chemical potential
    double binding_energy=stod(argv[3]);// binding energy
    int Ncolloid=stod(argv[4]);// number of colloids
    for(int i=0;i<4;i++){
        singlet_tot={};Naa_tot={};Nab_tot={};Nbb_tot={};Nba_tot={};
        vector<double> recording;
        double ns=strands[i];
        double number_of_iteration=35000000;// number of iterations
        System mySystem(Ncolloid,scaleL,temperature,chemical_potential,binding_energy,ns);
        mySystem.randomSetting(Ncolloid,mySystem.Ncell);// initialize the system and the position of colloids
        for(int q=0;q<number_of_iteration;q++){
            if(q>5000000 && q%100 == 0){// start the recording after 5 million iterations every 100 iterations
                recording=mySystem.statistics();
                singlet_tot.push_back(recording[0]);
                Naa_tot.push_back(recording[1]);
                Nbb_tot.push_back(recording[2]);
                Nab_tot.push_back(recording[3]);
                Nba_tot.push_back(recording[4]);
                if(q%1000000==0){
                    mySystem.takePicture(position);
                }
                
                // for(int i=0;i<Nc;i++){
                //     int coll=mySystem.check_overlap(i);
                //     if(coll ==1){
                //         cout<<"something is wrong"<<endl;
                //     }
                // }
            }
        mySystem.MCmove();// make a MC move of the colloid
        mySystem.adapt_step();// adapt the step size
        }
        singlet<<getAverage(singlet_tot)<<endl;
        Naa<<getAverage(Naa_tot)<<endl;
        Nab<<getAverage(Nab_tot)<<endl;
        Nbb<<getAverage(Nbb_tot)<<endl;
        Nba<<getAverage(Nba_tot)<<endl;
    }
    
    return 0;
}

double getAverage(vector<double>& sequence){// calculate the average of a vector
    int length=sequence.size();
    double sum=0,average;
    for(int i =0; i<length; i++){
        sum=sum+sequence[i];
    }
    average=sum/length;
    return average;
}

//done verified