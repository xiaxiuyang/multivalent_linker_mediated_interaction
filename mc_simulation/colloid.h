#include "array_computation.h"
#include "geometry.h"
#include <omp.h>

using namespace std;

mt19937 randNum(time(nullptr));


const double col_size=100;
const double rec_size=10;

double randCV(void){// generate a random number between 0 and 1
    return randNum()*(1.0/4294967296.0);
}

class Colloid{//define the class of colloids
public:
    double R, sigma,rc,ni,V_tot_r;

    double boxs[3];
    double V_free_r;
    double chi;
    double sum_m,sum_km;
    double bar_n,J;
    
    double backup_V_free_r;
    double backup_chi;
    double backup_sum_m, backup_sum_km;
    double backup_bar_n;

    int index,type;
    vector<double> V_b_ij,sum_q;
    vector<double> distances;
    vector<int> neighbors;

    double backup_V_b_ij,backup_sum_q;
    double backup_distance;
    int backup_neighbor;

    Colloid *previousPtr, *nextPtr;

    void initialize(void);
    void clear_vectors(void);
    void prepare_update(void);
    // void copy_old_to_new(void);
    // void copy_new_to_old(void);
    void random_move(double& ds);
    void cal_chi(double& binding_energy, double& rho0);
    void backup_and_erase(int& i,double& binding_energy, double& rho0);
    void reject_addition(int&i,double& binding_energy, double& rho0);
    void restore_to_backup(void);
    double repulsion(void);
    void random_position(void);
    Colloid(){
        R=col_size; rc=rec_size; ni=100; sigma=2*R;
        V_tot_r=cal_shell(R,rc); 
        V_free_r=V_tot_r;
        backup_V_free_r=V_free_r;
        
        chi=0; sum_m=0; sum_km=0;  bar_n=ni/2.; 
        backup_chi=0; backup_sum_m=0; backup_sum_km=0; backup_bar_n=ni/2.;
        
        V_b_ij={}; sum_q={}; distances={}; neighbors={};
        previousPtr=nullptr;nextPtr=nullptr;
    }
    void label(int indexing,int typeGiven){
        index=indexing;
        type=typeGiven;
    }
};

void Colloid::initialize(void){// initialize a colloid
    initialize_1d_array(boxs,0,3);
    V_free_r=V_tot_r;
    chi=0; sum_m=0; sum_km=0; bar_n=ni/2.;
    V_b_ij={}; sum_q={}; neighbors={}; distances={};
}

void Colloid::clear_vectors(void){ // clear all vectors in the colloid
    V_b_ij={}; sum_q={}; neighbors={}; distances={};
}

void Colloid:: prepare_update(void){
    V_free_r=V_tot_r;
    V_b_ij={}; sum_q={}; neighbors={}; distances={};
}

/*void Colloid::copy_old_to_new(void){
    copy(new_boxs,boxs);
    new_V_free_r=V_free_r;
    new_chi=chi; new_sum_m=sum_m; new_sum_km=sum_km; new_bar_n=bar_n;
    new_V_b_ij=V_b_ij; new_sum_q=sum_q; new_sum_kq=sum_kq; new_neighbors=neighbors;
}

void Colloid::copy_new_to_old(void){
    copy(boxs,new_boxs);
    V_free_r=new_V_free_r;
    chi=new_chi; sum_m=new_sum_m; sum_kq=new_sum_kq; bar_n=new_bar_n;
    V_b_ij=new_V_b_ij; sum_q=new_sum_q; sum_kq=new_sum_kq; neighbors=new_neighbors;
}*/

void Colloid::random_move(double& ds){//make a random move
    
    #pragma omp parallel for schedule(static,3)
    for(int i=0; i<3;i++){
        double abs_pos=(boxs[i]+2*ds*(randCV()-0.5));
        boxs[i]=abs_pos-floor(abs_pos);
    }
}

void Colloid::cal_chi(double& binding_energy, double& rho0){
    chi=exp(-1*binding_energy)/(rho0*V_free_r);
}

void Colloid::backup_and_erase(int& i,double& binding_energy, double& rho0){
    // make backup for trial moves
    vector<int>::iterator pos_i=find(neighbors.begin(),neighbors.end(),i);
    int position_i=pos_i-neighbors.begin();
    
    backup_V_free_r=V_free_r;
    backup_chi=chi;
    backup_sum_m=sum_m; backup_sum_km=sum_km;
    backup_bar_n=bar_n;
    
    backup_V_b_ij=V_b_ij[position_i];
    V_b_ij.erase(V_b_ij.begin()+position_i);
    backup_sum_q=sum_q[position_i];
    sum_q.erase(sum_q.begin()+position_i);
    
    backup_neighbor=i;
    neighbors.erase(neighbors.begin()+position_i);
    backup_distance=distances[position_i];
    distances.erase(distances.begin()+position_i);
    V_free_r=V_free_r+overlap(R+rc,R,backup_distance);
    cal_chi(binding_energy,rho0);

}

void Colloid::reject_addition(int& i,double& binding_energy, double& rho0){
    //reject the addition of a new colloid i to its neighbor info vector
    if(neighbors.back()==i){
        double distance=distances.back();
        V_b_ij.pop_back();
        sum_q.pop_back();
        neighbors.pop_back();
        distances.pop_back();
        V_free_r=V_free_r+overlap(R+rc,R,distance);
        cal_chi(binding_energy,rho0);
    }
}

void Colloid::restore_to_backup(void){// restore the backup info
    V_free_r=backup_V_free_r;
    chi=backup_chi;
    sum_m=backup_sum_m;
    bar_n=backup_bar_n;

    V_b_ij.push_back(backup_V_b_ij);
    sum_q.push_back(backup_sum_q);
    distances.push_back(backup_distance);
    neighbors.push_back(backup_neighbor);
}

double Colloid::repulsion(void){// repulsion from the receptors
    double result=0;
    result=-1*ni*log(V_free_r/V_tot_r);
    return result;
}

void Colloid::random_position(void){
    #pragma omp for schedule(static,3)
    for(int d=0;d<3;d++){
        boxs[d]=randCV();
    }
}