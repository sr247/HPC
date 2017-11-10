#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <ctime>

using namespace std;

bool isprime(int n);
void Exercice1(int, int);
void Exercice2(int, int);
void Mat2vecProduct(float*, float*, int, int, float*, int);

int main(int argc, char ** argv){
	
	//~ cout << argc << endl;
	if(argc < 4){
		if(atoi(argv[1]) == 1)
			Exercice1(atoi(argv[2]),atoi(argv[3]));
		if(atoi(argv[1]) == 2)
			Exercice2(atoi(argv[2]), atoi(argv[3]));
	}
	else
	{
		cout << "Too much arguments provided !" << endl;
	}
	return 0;
	
}




void Exercice1(int t, int n){	
	omp_set_num_threads(t);
	float* A = (float*) malloc (n*n*sizeof(float));
	float* X = (float*) malloc (n*sizeof(float));
	float* B = (float*) malloc (n*sizeof(float));
	
	#pragma omp parallel for
	for (int i=0; i<n; i++){
		X[i] = 1.0;
		B[i] = 0.0;
		for(int j=0; j<n; j++)
			A[i*n+j] = 1.0;
	}	
	//~ for (int i = 0; i < n; i++){
		//~ for(int j = 0; j < n; j++)
			//~ cout << " " << A[i] << " ";
		//~ cout << endl;
	//~ }
	clock_t tmp;
	Mat2vecProduct(A, X, n, n, B, t);
	tmp=clock();
	cout << " " << (double)tmp/CLOCKS_PER_SEC << "s " << endl;
}



void Mat2vecProduct(float *A, float *X, int n, int m, float *B, int t){
	
	//~ Profile time
	#pragma omp parallel for num_threads(t) //schelude(runtime);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			B[i] += A[i*n+j] * X[j];
		}
	}
	//~ Profile end 

	for(int i = 0; i < n; i++){
		//~ printf(" %lf ", B[i]);
		cout << " " << B[i] << " ";
	}
	cout << endl;
}



void Exercice2(int t, int n){
	
	double s;
	double fxi = 0.0;
	double res = 0.0;
	
    s = (1.0 / n);
    
    clock_t tmp;
    #pragma omp reduce(+:fxi) for num_threads(t)
    for(int i = 0; i < n; i++){				
		fxi += ((4.0 / (1.0 + pow(i*s, 2))) + (4.0 / (1.0 + pow((i+1)*s, 2))))/2.0;
	}
	tmp=clock();
	
	fxi = fxi * s;
	cout << "Pi = " << (double)fxi << std::endl;
	cout << " " << (double)tmp/CLOCKS_PER_SEC << "s " << endl;
}


void Exercice3(int t, int n){
	

	cout << isprime(n) << endl;
}




bool isprime(int n){
	
	int sqn = ceil(sqrt(n));
	for(int i=2; i<=sqn; i++){
		if(i % n == 0)
			return false;
	}	
	return true;
}
