#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>
using namespace std;

int clusteredMatMul(int, char**, int);
void printMat(float *, int);
void Mat2vecProduct(float *A, float *X, int n, int m, float *Y);

int main(int argc, char** argv) {
	
	clusteredMatMul(argc, argv, 0);
    return 0;
}

int clusteredMatMul(int argc, char **argv, int v){
	
	int n = 1024;
	int id, nbprocs, err;
	MPI_Status status;
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD , &id);
    err = MPI_Comm_size(MPI_COMM_WORLD , &nbprocs);		
	
	int size = n/nbprocs;
	float* A = (float*) malloc (n*n*sizeof(float));
    float* A_local = (float*) malloc (n*size*sizeof(float));
    float* B = (float*) malloc (n*n*sizeof(float));
    float* C = (float*) malloc (n*n*sizeof(float));
	float* C_local = (float*) malloc (n*size*sizeof(float));
    
    if(id == 0){
		for (unsigned int i = 0; i < n; i++){		
			for(unsigned int j=0; j < n;j++){
				A[i*n+j] = 1.f;
				B[i*n+j] = 1.f;
				C[i*n+j] = 0.f;
				if(i < size)
					C_local[i*size+j] = 0.f;
			}
		}
	}
	else{
		
		for (unsigned int i = 0; i < size; i++){		
			for(unsigned int j = 0; j < n; j++){
				C_local[i*size+j] = 0.f;
			}
		}
	}
	MPI_Bcast(B, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	MPI_Scatter(A, n*size, MPI_FLOAT, A_local, n*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	//~ printf("Haha %d %d\n", nbprocs, id);
	
	for(int i = 0; i < size; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
				//~ printf("%d %d %d\n", i, j, k);
				C_local[i*n+j] += A_local[i * size + k] * B[k*size+j];
			}
		}
	}
	
	MPI_Gather(C_local, n*size, MPI_FLOAT, C, n*size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	if(id == 0){
		if( n < 32)
			printMat(C, n);
	}
    MPI_Finalize();
	return err;
}


void Mat2MatProduct(float **A, float **B, float **C, int n, int m){
	
	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = 0; j < m; j++){
			for(unsigned int k = 0; k < n; k++){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void printMat(float * A, int n){
	for (int i = 0; i < n ; i++){
		for(int j = 0; j < n; j++)
			cout << " " << A[i*n+j] << " ";
		cout << endl;
	}
}
