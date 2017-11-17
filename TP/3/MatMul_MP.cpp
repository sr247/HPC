#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

using namespace std;

void Mat4MatProduct(float**, float**, float**, int, int);
void sequentialMatMult4(const int&);
void printMat(float **, int);


int main(int argc, char** argv) {

	sequentialMatMult4(0);
	
    return 0;
}

void Mat2MatProduct(float **A, float **B, float **C, int n, int m){
	
	#pragma opm parallel for
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			for(int k = 0; k < n; k++){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}



void sequentialMatMult4(const int &v){
	
	int n = 16;
	float** A = (float**) malloc (n*sizeof(float*));
    float** B = (float**) malloc (n*sizeof(float*));
    float** C = (float**) malloc (n*sizeof(float*));
    
	for (int i = 0; i < n; i++){
		A[i] = (float*) malloc (n*sizeof(float));
		B[i] = (float*) malloc (n*sizeof(float));
		C[i] = (float*) malloc (n*sizeof(float));
	}
    
    for (int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			A[i][j] = 1.0;
			B[i][j] = 1.0;
			C[i][j] = 0.0;
		}
	}
	
    clock_t start, end;
	start = clock();
	Mat2MatProduct(A, B, C, n, n);
	end = clock();
	cout << "Time elasped:" << float(end - start)/float(CLOCKS_PER_SEC) << endl;
	
	if(n < 32){
		printMat(C, n);
	}
	free(A);
	free(B);
	free(C);
	
}


void printMat(float ** A, int n){
	for (int i = 0; i < n ; i++){
		for(int j = 0; j < n; j++)
			cout << " " << A[i][j] << " ";
		cout << endl;
	}
}
