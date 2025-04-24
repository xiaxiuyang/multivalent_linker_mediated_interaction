#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>


using namespace std;


double overlap(double R1,double R2,double d){ // calculate the overlap volume of 2 spheres
    double result=0;
    if(d > R1+R2){
        result=0.;
    }else if(d<fabs(R1-R2)){
        result=4*M_PI*pow(min(R1,R2),3)/(double)3;
    }else{
        result= M_PI*pow(R1+R2-d,2)*(pow(d,2)+2.*d*(R1+R2)-3.*(pow(R1,2)+pow(R2,2))+6.*R1*R2)/(12.*d);
    }
    return result;
}

double shell_overlap(double R1,double R2,double thickness1,double thickness2,double d){// overlap of 2 shells
    double result=0.;
    result=overlap(R1+thickness1,R2+thickness2,d)-overlap(R1+thickness1,R2,d)-overlap(R1,R2+thickness2,d)+overlap(R1,R2,d);
    return result;
}

double cal_shell(double R,double delta){// calculate the volume of a shell
    return (4.*M_PI/3.)*(pow(R+delta,3)-pow(R,3));
}

double cal_ex_shell( double R1,double R2,double thickness1,double d){
    return overlap(R1+thickness1,R2,d)-overlap(R1,R2,d);
}

