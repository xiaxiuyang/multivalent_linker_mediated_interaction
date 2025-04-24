#include "colloid.h"

using namespace std;

// const int Nc=360;

// const int Nc=540;

class Cell_List {
public:
    Colloid *headPtr; // head pointer that points to the first particle in the cell
    int index; // size of the cell

    void insert(Colloid *pNewPtr) { // a function to insert a particle pointed by pNewPtr in a chosen cell
    // cannot insert a particle in a corresponding cell automatically, but need to pick a cell first
    // the inserted particle is pointed by the head pointer
        if (headPtr != nullptr) {
            pNewPtr->nextPtr = headPtr;
            headPtr->previousPtr=pNewPtr;
        } else {
            pNewPtr->nextPtr = nullptr;
        }
        pNewPtr->previousPtr = nullptr;
        headPtr = pNewPtr;

    }

    void remove(Colloid *pOldPtr) {// a function to remove a particle from a cell
    // need to pick a cell first to make it work
    // to cut the linked list and to link the two neighbor particle around the removed particle
        if (pOldPtr != headPtr){
           (pOldPtr->previousPtr)->nextPtr = pOldPtr->nextPtr;
            
        }else{
            headPtr = pOldPtr->nextPtr;
        }

        if (pOldPtr->nextPtr != nullptr) {
            (pOldPtr->nextPtr)->previousPtr = pOldPtr->previousPtr;
        }

        pOldPtr->nextPtr = nullptr;
        pOldPtr->previousPtr = nullptr;
     

    }

};

class System{// in this system N=2, and the size of the particles are determined
public:
    // double V;
    double T,mu_L;
    double lunit;
    double F;
    double ns,V_tot_L;
    double rho_0;
    double G;
    int Ncell;
    double cellLength;
    double ds, Nacc, Ntot;
    int Nc;
    vector<Colloid> myColloids;
    vector<vector<vector<Cell_List>>> cells;
    
    double backup_F;
    double temp_boxs[3];
    double temp_V_free_r;
    double temp_chi;
    double temp_bar_n;

    vector<double> temp_V_b_ij,temp_distances,temp_sum_q;
    vector<int> temp_neighbors;
    
    void transform_n_to_J(void);
    void transform_J_to_n(void);
    double cal_mutual_volume(int& i, int& j,double& distance);
    double cal_ex_volume(int& i,int& j,double& distance);
    void cal_sum_km(int& i);
    void cal_sum_kq(int& i, int& j);
    vector<double> sum_km_od(int& i);
    vector<double> sum_kq_od(int& i, int& order_j);
    void cal_sum_m(int& i);
    void cal_sum_q(int& i,int& j);
   
    double cal_sum_Q_iL(int& i);
    void remove_from_cell(int& i);
    void extract_colloid(int& i);
    void put_in_cell(int& i);
    void update_mutual_info(int& i,int& j,double& distance);
    int update_neighbors(int& i);
    void extract_colloid_nobackup(int& i);
    void move_colloid_back(int& i);
    void reject_MC(int& i);
    void find_solution(void);
    void energy(void);
    void MCmove(void);
    void adapt_step(void);
    vector<double> function_and_derivative(int& i);
    int update_own_info(int& i);
    int check_overlap(int& i);
    vector<int> neighbors(int& i);
    
    void randomSetting(int number,int Ncell);

    vector<double> statistics(void);
    void takePicture(ofstream& position);
    System(int Ncolloid, double scaleL,double temperature, double chemPt_L,double binding_eng,double strands){
        ds=0.01; // later adapt it
        // V=volume_box;
        Nc=Ncolloid;//set the number of colloids
        vector<Colloid> temporary(Nc);
        myColloids=temporary;// create vector of colloids
        lunit=scaleL;// size of the system
        T=temperature;//
        mu_L=chemPt_L;//chemical potential
        ns=strands;
        rho_0=0.6022;
        G=binding_eng;
        Ncell=floor(lunit/(2.*col_size+2.*rec_size+5)); // number of cells in one direction
        cellLength=1./((double) Ncell);
        F=0;
        Nacc=0;Ntot=0;
    }
    
    // void setColloids(double distance,string type_difference){
    //     double sd=0;
    //     sd=distance/lunit;
    //     myColloids[0].boxs[0]=0.5-sd/2.;
    //     myColloids[1].boxs[0]=0.5+sd/2.;
    //     myColloids[0].boxs[1]=0.5;
    //     myColloids[1].boxs[1]=0.5;
    //     myColloids[0].boxs[2]=0.5;
    //     myColloids[1].boxs[2]=0.5;
    //     int i=0,j=1;
    //     myColloids[i].index=i; myColloids[j].index=j;
    //     if(type_difference=="same"){
    //         myColloids[i].type=0;myColloids[j].type=0;
    //     }else{
    //         myColloids[i].type=0;myColloids[j].type=1;
    //     }
    //     if(distance<220){
    //         myColloids[i].neighbors.push_back(j);
    //         myColloids[j].neighbors.push_back(i);
    //     }
        
    //     update_mutual_info(i,j,distance);
    //     myColloids[i].cal_chi(G,rho_0);
        
    //     myColloids[i].sum_q.push_back(0);
    //     myColloids[j].sum_q.push_back(0);
        
    // }
};

void System::randomSetting(int number,int number_cell){
    Ncell=number_cell;
    cellLength=1./((double) Ncell);
    vector<vector<vector<Cell_List>>> cellsCopy(Ncell,vector<vector<Cell_List>>(Ncell,vector<Cell_List>(Ncell)));
    cells=cellsCopy;
    for(int i =0;i<Ncell;i++){// initialize the cell matrix
        for(int j=0;j<Ncell;j++){
            for(int k=0;k<Ncell;k++){
                cells[i][j][k].headPtr=nullptr;
                cells[i][j][k].index=i*pow(Ncell,2)+j*Ncell+k;
            }
        }
    }// set cells
    for(int i=0;i<Nc;i++){
        myColloids[i].index=i;
        myColloids[i].type = i%2;
        myColloids[i].previousPtr=nullptr;
        myColloids[i].nextPtr=nullptr;
        int coll=1;
        while (coll == 1){
            myColloids[i].random_position();
            put_in_cell(i);
            coll=check_overlap(i);
            if(coll == 1){
                remove_from_cell(i);
            }
        }
    }// initialize the position of the colloids without overlap
    for(int i=0;i<Nc;i++){
        int coll=update_own_info(i);
        if(coll == 1){
            cout<<"something is wrong"<<endl;
        }
    }
    energy();
}

void System::transform_n_to_J(void){
    for(int i=0;i<Nc;i++){
        myColloids[i].J=myColloids[i].bar_n*myColloids[i].chi+1.;
    }
}

void System::transform_J_to_n(void){
    for(int i=0;i<Nc;i++){
        myColloids[i].bar_n=(myColloids[i].J-1)/myColloids[i].chi;
    }
}

double System::cal_mutual_volume(int& i, int& j,double& distance){
    double result=0;
    result=shell_overlap(myColloids[i].R,myColloids[j].R,myColloids[i].rc,myColloids[j].rc,distance);
    return result;
}

double System::cal_ex_volume(int& i, int& j,double& distance){
    double result;
    result=overlap(myColloids[i].R+myColloids[i].rc,myColloids[j].R,distance);
    return result;
}

void System::cal_sum_km(int& i){
    myColloids[i].sum_km=myColloids[i].V_free_r*exp(mu_L/T)*ns*
                    myColloids[i].bar_n*myColloids[i].chi*
                    pow(myColloids[i].bar_n*myColloids[i].chi+1,ns-1);
}

// void System::cal_sum_kq(int& i, int& j){
//     vector<int>::iterator pos_j=find(myColloids[i].neighbors.begin(),myColloids[i].neighbors.end(),j);
//     int position_j=pos_j-myColloids[i].neighbors.begin();
    
//     if(myColloids[i].type == myColloids[j].type){
//         myColloids[i].sum_kq[position_j]=myColloids[i].V_b_ij[position_j]*exp(mu_L/T)*ns*myColloids[i].bar_n*
//                         myColloids[i].chi*
//                         (pow(myColloids[i].bar_n*myColloids[i].chi+myColloids[j].bar_n*myColloids[j].chi,ns-1)
//                         - pow(myColloids[i].bar_n*myColloids[i].chi+1,ns-1));
//     }else{
//         myColloids[i].sum_kq[position_j]=myColloids[i].V_b_ij[position_j]*exp(mu_L/T)*ns*myColloids[i].bar_n*
//                         pow(myColloids[i].bar_n*myColloids[i].chi+1,ns-1)*
//                         (pow(myColloids[j].bar_n*myColloids[j].chi+1,ns)-1);
//     }
    
// }

vector<double> System::sum_km_od(int& i){
    double result1=myColloids[i].V_free_r*exp(mu_L/T)*ns*
                    myColloids[i].bar_n*myColloids[i].chi*
                    pow(myColloids[i].bar_n*myColloids[i].chi+1,ns-1);
    double result2=result1/myColloids[i].bar_n
                +result1*(ns-1)*myColloids[i].chi/(myColloids[i].bar_n*myColloids[i].chi+1);//satisfied only by ns>2
    return {result1,result2};
}

vector<double> System::sum_kq_od(int& i, int& order_j){
    double result1,result2;
    double part1,part2;
    double Ji,Jj;
    double factor;
    int j=myColloids[i].neighbors[order_j];
    factor=myColloids[i].V_b_ij[order_j]*exp(mu_L/T)*ns*myColloids[i].bar_n*myColloids[i].chi;
    Ji=myColloids[i].bar_n*myColloids[i].chi+1;
    Jj=myColloids[j].bar_n*myColloids[j].chi+1;
    if(myColloids[i].type == myColloids[j].type){
        part1=factor*pow(Ji+Jj-1.,ns-1);
        part2=factor*pow(Ji,ns-1);
        result1=part1-part2;
        result2=result1/myColloids[i].bar_n+(ns-1)*myColloids[i].chi*(
            part1/(Ji+Jj-1)-part2/(Ji)
        );
    }else{
        result1=factor*pow(Ji,ns-1)*(pow(Jj,ns)-1);
        result2=result1/myColloids[i].bar_n+result1*(ns-1)*myColloids[i].chi/(Ji);
    }
    return {result1,result2};
}

void System::cal_sum_m(int& i){
    myColloids[i].sum_m=myColloids[i].V_free_r*exp(mu_L/T)*
                    (pow(myColloids[i].bar_n*myColloids[i].chi+1,ns)-1);   
}

void System::cal_sum_q(int& i,int& j){
    vector<int>::iterator pos_j=find(myColloids[i].neighbors.begin(),myColloids[i].neighbors.end(),j);
    int position_j=pos_j-myColloids[i].neighbors.begin();
    vector<int>::iterator pos_i=find(myColloids[j].neighbors.begin(),myColloids[j].neighbors.end(),i);
    int position_i=pos_i-myColloids[j].neighbors.begin();
    double Ji,Jj;
    double factor=myColloids[i].V_b_ij[position_j]*exp(mu_L/T);
    Ji=myColloids[i].bar_n*myColloids[i].chi+1;
    Jj=myColloids[j].bar_n*myColloids[j].chi+1;

    if(myColloids[i].type == myColloids[j].type){
        myColloids[i].sum_q[position_j]=factor*
                 (pow(Ji+Jj-1,ns)-pow(Ji,ns)-pow(Jj,ns)+1);
    }else{
        myColloids[i].sum_q[position_j]=factor*
                 (pow(Ji,ns)*pow(Jj,ns)-pow(Ji,ns)-pow(Jj,ns)+1);
    }

    myColloids[j].sum_q[position_i]=myColloids[i].sum_q[position_j];
}


double System::cal_sum_Q_iL(int& i){
    double Q_iL=0;
    for(int order_j=0;order_j<myColloids[i].neighbors.size();order_j++){
        if(i>myColloids[i].neighbors[order_j]){
            cal_sum_q(i,myColloids[i].neighbors[order_j]);
            Q_iL+=myColloids[i].sum_q[order_j];
        }
    }
    return Q_iL;
}

void System::remove_from_cell(int& i){ // remove the colloid from the cell
    int cell_x,cell_y,cell_z;
    cell_x=floor(myColloids[i].boxs[0]/cellLength);
    cell_y=floor(myColloids[i].boxs[1]/cellLength);
    cell_z=floor(myColloids[i].boxs[2]/cellLength);
    cells[cell_x][cell_y][cell_z].remove(&myColloids[i]);
}

void System::extract_colloid(int& i){// extract the colloid from the cell and update its neighbors' information
    copy(temp_boxs,myColloids[i].boxs);
    temp_V_free_r=myColloids[i].V_free_r;
    temp_chi=myColloids[i].chi;
    temp_bar_n=myColloids[i].bar_n;

    temp_V_b_ij=myColloids[i].V_b_ij;
    temp_distances=myColloids[i].distances;
    temp_neighbors=myColloids[i].neighbors;
    temp_sum_q=myColloids[i].sum_q;
    for(int order_j=0;order_j<myColloids[i].neighbors.size();order_j++){
        myColloids[myColloids[i].neighbors[order_j]].backup_and_erase(i,G,rho_0);
    }
    remove_from_cell(i);
    myColloids[i].prepare_update();
    
    
}

void System::put_in_cell(int& i){// put the colloid in the cell
    int cell_x,cell_y,cell_z;
    cell_x=floor(myColloids[i].boxs[0]/cellLength);
    cell_y=floor(myColloids[i].boxs[1]/cellLength);
    cell_z=floor(myColloids[i].boxs[2]/cellLength);
    cells[cell_x][cell_y][cell_z].insert(&myColloids[i]);
}

void System::update_mutual_info(int& i,int& j,double& distance){// update the information between two colloids
    // colloid i is the principal colloid
    double core_exclusion;
    double mutual_volume;
    
    core_exclusion=cal_ex_volume(i,j,distance);
    mutual_volume=cal_mutual_volume(i,j,distance);
    myColloids[i].V_free_r=myColloids[i].V_free_r-core_exclusion;
    myColloids[j].V_free_r=myColloids[j].V_free_r-core_exclusion;
    myColloids[i].neighbors.push_back(j);
    myColloids[j].neighbors.push_back(i);
    myColloids[i].distances.push_back(distance);
    myColloids[j].distances.push_back(distance);
    myColloids[i].V_b_ij.push_back(mutual_volume);
    myColloids[j].V_b_ij.push_back(mutual_volume);
    myColloids[i].sum_q.push_back(1);
    myColloids[j].sum_q.push_back(1);
    myColloids[j].cal_chi(G,rho_0);
}

int System::update_neighbors(int& i){// update info of all the neighbors of the colloid i
    int neighborX,neighborY,neighborZ;
    int cell_x,cell_y,cell_z;
    double distance;
    Colloid *theCol=&myColloids[i];
    double dr[3],iter_box_adj[3];
    int periodic_adj[3];
    cell_x=floor(theCol->boxs[0]/cellLength);
    cell_y=floor(theCol->boxs[1]/cellLength);
    cell_z=floor(theCol->boxs[2]/cellLength);
    for(Colloid* iter= theCol->nextPtr; iter !=nullptr; iter=iter->nextPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }else if(distance<theCol->sigma+2.*theCol->rc){
            
            update_mutual_info(i,iter->index,distance);
            
        }
    }
    for(Colloid* iter= theCol->previousPtr; iter !=nullptr; iter=iter->previousPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }else if(distance<theCol->sigma+2.*theCol->rc){
        
            update_mutual_info(i,iter->index,distance);
     
        }
    }
    for(int q=-1; q<2;q++)// check its neighboring cells
    for(int n=-1; n<2;n++)
    for(int p=-1; p<2;p++){
        if(q!=0 || n!=0 || p !=0){
            
            neighborX=cell_x+q;neighborY=cell_y+n;neighborZ=cell_z+p;
            periodic_adj[0]=0;periodic_adj[1]=0;periodic_adj[2]=0;

            if(neighborX<0) periodic_adj[0]=-1;
            else if(neighborX>=Ncell) periodic_adj[0]=1;

            if(neighborY<0) periodic_adj[1]=-1;
            else if(neighborY>=Ncell) periodic_adj[1]=1;
                
            if(neighborZ<0) periodic_adj[2]=-1;
            else if(neighborZ>=Ncell) periodic_adj[2]=1;

            neighborX-=Ncell*periodic_adj[0]; neighborY-=Ncell*periodic_adj[1]; neighborZ-=Ncell*periodic_adj[2];
            
            Colloid* ngbPtr= cells[neighborX][neighborY][neighborZ].headPtr;
            for(Colloid* iter=ngbPtr;iter!=nullptr;iter=iter->nextPtr){
                periodic_boundary(iter_box_adj,iter->boxs,periodic_adj);
                Substract(dr,iter_box_adj,theCol->boxs);
                distance=Norm(dr)*lunit;
                
                if(distance<theCol->sigma){
                    return 1;
                }else if(distance<theCol->sigma+2.*theCol->rc){
                    
                    update_mutual_info(i,iter->index,distance);
    
                }
            }  
        }
    }
    myColloids[i].cal_chi(G,rho_0);
    return 0; // no collision
}

void System::extract_colloid_nobackup(int& i){// retract from the trial cell
    int neighbor;
    for(int order_j=0;order_j<myColloids[i].neighbors.size();order_j++){
        neighbor=myColloids[i].neighbors[order_j];
        myColloids[neighbor].reject_addition(i,G,rho_0);
    }
    myColloids[i].prepare_update();
    remove_from_cell(i);
}

void System::move_colloid_back(int& i){// put the colloid back to its original position
    int neighbor;
    copy(myColloids[i].boxs,temp_boxs);
    myColloids[i].V_free_r=temp_V_free_r;
    myColloids[i].chi=temp_chi;
    myColloids[i].bar_n=temp_bar_n;

    myColloids[i].V_b_ij=temp_V_b_ij;
    myColloids[i].distances=temp_distances;
    myColloids[i].neighbors=temp_neighbors;
    myColloids[i].sum_q=temp_sum_q;
    for(int order_j=0;order_j<myColloids[i].neighbors.size();order_j++){
        neighbor=myColloids[i].neighbors[order_j];
        myColloids[neighbor].restore_to_backup();
    }
    put_in_cell(i);
}

void System::reject_MC(int& i){// reject the MC move
    extract_colloid_nobackup(i);
    move_colloid_back(i);
}

void System::find_solution(void){ // solve the Nc self-consistent equations
    
    vector<double> medium(Nc);
    vector<double> error(Nc);
    double max_error=1;

    while (max_error>0.0000001){

        #pragma omp parallel for
        for(int i=0;i<Nc;i++){
            vector<double> fun_der=function_and_derivative(i);
            medium[i]=myColloids[i].bar_n-fun_der[0]/fun_der[1];
            error[i]=fabs(medium[i]-myColloids[i].bar_n);

        }
        auto max_it=max_element(error.begin(),error.end());
        max_error=*max_it;
        #pragma omp parallel for
        for(int i=0;i<Nc;i++){
            myColloids[i].bar_n=medium[i];
        }

    }
 
}

void System::energy(void){// give the energy of the system
    double Fatt,Frep;
    Frep=0; Fatt=0;
    find_solution();
    #pragma omp parallel for reduction(+:Fatt,Frep)
    for(int i=0;i<Nc;i++){
        Frep+=myColloids[i].repulsion();
        cal_sum_m(i);
        Fatt+=(myColloids[i].ni*log(myColloids[i].bar_n)-myColloids[i].bar_n-myColloids[i].sum_m-cal_sum_Q_iL(i));
    }
    F=Fatt+Frep;
    // cout<<"absorbed: "<<myColloids[0].sum_m<<endl;
    // int i=1;
    // cout<<"bridged: "<<cal_sum_Q_iL(i)<<endl;
    
}

void System::MCmove(void){// do a MC
    int chosen_index=randNum()%Nc;
    int collision;
    Ntot++;
    extract_colloid(chosen_index);
    myColloids[chosen_index].random_move(ds);
    put_in_cell(chosen_index);
    collision=update_neighbors(chosen_index);
    if(collision ==1){
        reject_MC(chosen_index);
    }else{
        backup_F=F;
        energy();
        if(randCV()<exp(-1*(F-backup_F)/T)){
            Nacc++;
            
        }else{
            reject_MC(chosen_index);
            F=backup_F;
        }

    }


}

void System::adapt_step(void){
    if(Ntot>200){
        double acceptance_rate;
        acceptance_rate=Nacc/Ntot;
        if(acceptance_rate>0.6 && ds<0.03){
            ds=ds*1.01;
        }else if(acceptance_rate<0.3 && ds>0.0003){
            ds=ds/1.01;
        }
        Nacc=0;Ntot=0;
    }
}

vector<double> System::function_and_derivative(int& i){
    vector<double> result={0,0};
    vector<double> km,kq,each_kq;
    km={0,0}; kq={0,0};
    km=sum_km_od(i);
    for(int q=0;q<myColloids[i].neighbors.size();q++){
        each_kq=sum_kq_od(i,q);
        kq[0]+=each_kq[0];
        kq[1]+=each_kq[1];
    }
    result[0]=km[0]+kq[0]+myColloids[i].bar_n-myColloids[i].ni;
    result[1]=1+km[1]+kq[1];
    return result;
}

int System::update_own_info(int& i){
    int neighborX,neighborY,neighborZ;
    int cell_x,cell_y,cell_z;
    double distance;
    Colloid *theCol=&myColloids[i];
    double dr[3],iter_box_adj[3];
    int periodic_adj[3];
    double core_exclusion,mutual_volume;
    cell_x=floor(theCol->boxs[0]/cellLength);
    cell_y=floor(theCol->boxs[1]/cellLength);
    cell_z=floor(theCol->boxs[2]/cellLength);
    for(Colloid* iter= theCol->nextPtr; iter !=nullptr; iter=iter->nextPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }else if(distance<theCol->sigma+2.*theCol->rc){
            
            core_exclusion=cal_ex_volume(i,iter->index,distance);
            mutual_volume=cal_mutual_volume(i,iter->index,distance);
            myColloids[i].V_free_r=myColloids[i].V_free_r-core_exclusion;

            myColloids[i].distances.push_back(distance);
            myColloids[i].V_b_ij.push_back(mutual_volume);
            myColloids[i].neighbors.push_back(iter->index);
            myColloids[i].sum_q.push_back(1);
           
        }
    }
    for(Colloid* iter= theCol->previousPtr; iter !=nullptr; iter=iter->previousPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }else if(distance<theCol->sigma+2.*theCol->rc){

            core_exclusion=cal_ex_volume(i,iter->index,distance);
            mutual_volume=cal_mutual_volume(i,iter->index,distance);
            myColloids[i].V_free_r=myColloids[i].V_free_r-core_exclusion;
            
            myColloids[i].neighbors.push_back(iter->index);
            myColloids[i].distances.push_back(distance);
            myColloids[i].V_b_ij.push_back(mutual_volume);
            myColloids[i].sum_q.push_back(1);
            
        }
    }
    for(int q=-1; q<2;q++)// check its neighboring cells
    for(int n=-1; n<2;n++)
    for(int p=-1; p<2;p++){
        if(q!=0 || n!=0 || p !=0){
            
            neighborX=cell_x+q;neighborY=cell_y+n;neighborZ=cell_z+p;
            periodic_adj[0]=0;periodic_adj[1]=0;periodic_adj[2]=0;

            if(neighborX<0) periodic_adj[0]=-1;
            else if(neighborX>=Ncell) periodic_adj[0]=1;

            if(neighborY<0) periodic_adj[1]=-1;
            else if(neighborY>=Ncell) periodic_adj[1]=1;
                
            if(neighborZ<0) periodic_adj[2]=-1;
            else if(neighborZ>=Ncell) periodic_adj[2]=1;

            neighborX-=Ncell*periodic_adj[0]; neighborY-=Ncell*periodic_adj[1]; neighborZ-=Ncell*periodic_adj[2];
            
            Colloid* ngbPtr= cells[neighborX][neighborY][neighborZ].headPtr;
            for(Colloid* iter=ngbPtr;iter!=nullptr;iter=iter->nextPtr){
                periodic_boundary(iter_box_adj,iter->boxs,periodic_adj);
                Substract(dr,iter_box_adj,theCol->boxs);
                distance=Norm(dr)*lunit;
                if(distance<theCol->sigma){
                    return 1;
                }else if(distance<theCol->sigma+2.*theCol->rc){
                    
                    core_exclusion=cal_ex_volume(i,iter->index,distance);
                    mutual_volume=cal_mutual_volume(i,iter->index,distance);
                    myColloids[i].V_free_r=myColloids[i].V_free_r-core_exclusion;
                      
                    myColloids[i].neighbors.push_back(iter->index);  
                    myColloids[i].distances.push_back(distance);
                    myColloids[i].V_b_ij.push_back(mutual_volume);
                    myColloids[i].sum_q.push_back(1);
            
                }
            }  
        }
    }
    myColloids[i].cal_chi(G,rho_0);
    return 0; // no collision

}

int System::check_overlap(int& i){// find the overlap volume of colloid i
    int neighborX,neighborY,neighborZ;
    int cell_x,cell_y,cell_z;
    double distance;
    Colloid *theCol=&myColloids[i];
    double dr[3], iter_box_adj[3];
    int periodic_adj[3];
    cell_x=floor(theCol->boxs[0]/cellLength);
    cell_y=floor(theCol->boxs[1]/cellLength);
    cell_z=floor(theCol->boxs[2]/cellLength);

    for(Colloid* iter= theCol->nextPtr; iter !=nullptr; iter=iter->nextPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }
    }
    for(Colloid* iter= theCol->previousPtr; iter !=nullptr; iter=iter->previousPtr){
        Substract(dr,iter->boxs,theCol->boxs);
        distance=Norm(dr)*lunit;
        if(distance<theCol->sigma){
            return 1;
        }
    }
    for(int q=-1; q<2;q++)// check its neighboring cells
    for(int n=-1; n<2;n++)
    for(int p=-1; p<2;p++){
        if(q!=0 || n!=0 || p !=0){
            
            neighborX=cell_x+q;neighborY=cell_y+n;neighborZ=cell_z+p;
            periodic_adj[0]=0;periodic_adj[1]=0;periodic_adj[2]=0;

            if(neighborX<0) periodic_adj[0]=-1;
            else if(neighborX>=Ncell) periodic_adj[0]=1;

            if(neighborY<0) periodic_adj[1]=-1;
            else if(neighborY>=Ncell) periodic_adj[1]=1;
                
            if(neighborZ<0) periodic_adj[2]=-1;
            else if(neighborZ>=Ncell) periodic_adj[2]=1;

            neighborX-=Ncell*periodic_adj[0]; neighborY-=Ncell*periodic_adj[1]; neighborZ-=Ncell*periodic_adj[2];
            
            Colloid* ngbPtr= cells[neighborX][neighborY][neighborZ].headPtr;
            for(Colloid* iter=ngbPtr;iter!=nullptr;iter=iter->nextPtr){
                periodic_boundary(iter_box_adj,iter->boxs,periodic_adj);
                Substract(dr,iter_box_adj,theCol->boxs);
                distance=Norm(dr)*lunit;
                
                if(distance<theCol->sigma){
                    return 1;
                }
            }  
        }
    }
    return 0; // no collision
}

vector<double> System::statistics(void){// do statistics on singlets, interaction between different types of colloids
    double Na_bound,Nb_bound;
    Na_bound=0;Nb_bound=0;
    int Nsing,Naa,Nbb,Nab,Nba;
    Nsing=0;Naa=0;Nbb=0;Nab=0;Nba=0;
    for(int i=0;i<Nc;i++){
        if(myColloids[i].neighbors.size()==0){
            Nsing++;
        }else{
            if(myColloids[i].type==0){
                Na_bound++;
            }else{
                Nb_bound++;
            }
            for(int order_j=0;order_j<myColloids[i].neighbors.size();order_j++){
                int neighbor=myColloids[i].neighbors[order_j];
                if(myColloids[i].type==0){
                    if(myColloids[neighbor].type==0){
                        Naa++;
                    }else{
                        Nab++;
                    }
                }else{
                    if(myColloids[neighbor].type==0){
                        Nba++;
                    }else{
                        Nbb++;
                    }
                }

            }
        }
    }
    double Aaa,Abb,Aab,Aba;
    Aaa=(double)Naa/(Na_bound);
    Abb=(double)Nbb/(Nb_bound);
    Aab=(double)Nab/(Na_bound);
    Aba=(double)Nba/(Nb_bound);
    return {(double)Nsing,Aaa,Abb,Aab,Aba};
}

void System::takePicture(ofstream& position){
    position<<Nc<<endl;
    position<<" "<<endl;
    for(int i=0;i<Nc;i++){
        if(myColloids[i].type==0){
            position<<"Cs"<<" "<<myColloids[i].boxs[0]*lunit<<" "<<myColloids[i].boxs[1]*lunit<<" "<<myColloids[i].boxs[2]*lunit<<endl;
        }else if(myColloids[i].type == 1){
            position<<"Cl"<<" "<<myColloids[i].boxs[0]*lunit<<" "<<myColloids[i].boxs[1]*lunit<<" "<<myColloids[i].boxs[2]*lunit<<endl;
        }
        
    }
}

