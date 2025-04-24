// #include <iostream>
#include <cmath>
#include <algorithm>


using namespace std;

void copy(double x1[],double x2[])
{
	for(int i = 0 ; i < 3; ++i)
	{
		x1[i] = x2[i];
	}
}
void Plus(double a[],double b[], double c[]){
    for( int i=0;i<3;i++){
        a[i]=b[i]+c[i];
    }
}
void Substract(double a[],double b[],double c[]){
	for(int i=0;i<3;i++){
		a[i]=b[i]-c[i];
	}
}
double Norm(double a[]){
	double result=0;
	
	for(int i=0;i<3;i++){
        result+= pow(a[i],2);
	}
	return sqrt(result);
}

void initialize_1d_array(double a[],double initial_value,int number){
    for(int i=0;i<number;i++){
        a[i]=initial_value;
    }
}

void periodic_boundary(double mir_pos[],double org_pos[], int adj[]){
    for(int i=0;i<3;i++){
        mir_pos[i]=org_pos[i]+adj[i];
	}
}