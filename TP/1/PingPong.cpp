#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>
using namespace std;

int Exercice1(int, char **);
int Exercice2(int, char **);
int Exercice3(int, char **, int);
int Exercice4(int, char **, int);
void Mat2vecProduct(float *A, float *X, int n, int m, float *Y);

int main(int argc, char** argv) {

	//~ Exercice1(argc, argv);
	
	//~ Exercice2(argc, argv);
	
	//~ Exercice3(argc, argv, 1000000);
	
	Exercice4(argc, argv, 10);

	

    return 0;
}

int Exercice1(int argc, char **argv){	
	int id, nbprocs;
	int msg;
	int err;
	MPI_Status status;
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD , &id);
    err = MPI_Comm_size(MPI_COMM_WORLD , &nbprocs);    
    if(id % 2 == 0){
		err = MPI_Send(&id, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD);		
		err = MPI_Recv(&msg, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD, &status);
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Send " << id << " to " << id+1 << std::endl 
		<< "Receive " << msg << " from " << id+1 << std::endl;
		std::cout << std::endl;
	}
	else
	{
		err = MPI_Recv(&msg, 1, MPI_INT, id-1, 0, MPI_COMM_WORLD, &status);
		err = MPI_Send(&id, 1, MPI_INT, id-1, 0, MPI_COMM_WORLD);
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Send " << id << " to " << id-1 << std::endl 
		<< "Receive " << msg << " from " << id-1 << std::endl;
		std::cout << std::endl;
	}    
    
    MPI_Finalize();
	return err;
}

int Exercice2(int argc, char **argv){
	int id, nbprocs;
	int msg;
	int err;
	MPI_Status status;
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD , &id);
    err = MPI_Comm_size(MPI_COMM_WORLD , &nbprocs);    
    if(id == 0){
		err = MPI_Send(&id, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD);		
		err = MPI_Recv(&msg, 1, MPI_INT, nbprocs-1, 0, MPI_COMM_WORLD, &status);
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Send " << id << " to " << id+1 << std::endl 
		<< "Receive " << msg << " from " << nbprocs-1 << std::endl;
		std::cout << std::endl;
	}
	else if(id < nbprocs-1)
	{
		err = MPI_Recv(&msg, 1, MPI_INT, id-1, 0, MPI_COMM_WORLD, &status);
		err = MPI_Send(&id, 1, MPI_INT, id+1, 0, MPI_COMM_WORLD);
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Send " << id << " to " << id+1 << std::endl 
		<< "Receive " << msg << " from " << id-1 << std::endl;
		std::cout << std::endl;
	}
	else
	{
		err = MPI_Recv(&msg, 1, MPI_INT, id-1, 0, MPI_COMM_WORLD, &status);
		err = MPI_Send(&id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Send " << id << " to " << 0 << std::endl 
		<< "Receive " << msg << " from " << id-1 << std::endl;
		std::cout << std::endl;
	}
    
    MPI_Finalize();
	return err;
}

int Exercice3(int argc, char **argv, int n){
	
	int id, nbprocs;
	double s;
	double fxi = 0.0;
	double res = 0.0;
	double offset = 0.01;
	int err;
	MPI_Status status;
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD , &id);
    err = MPI_Comm_size(MPI_COMM_WORLD , &nbprocs);
    
    double i0, imax;
    s = (1.0 / n);
    i0 = id * n / nbprocs;
    imax = i0 + n / nbprocs;
    
    for(int i = i0; i < imax; i++){				
		fxi += ((4.0 / (1.0 + pow(i*s, 2))) + (4.0 / (1.0 + pow((i+1)*s, 2))))/2.0;
	}
	
	MPI_Reduce(&fxi, &res, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(id == 0){
		res = res * s;
		std::cout << "Hello ! " << "I'm the thread " << id << " of " << nbprocs << std::endl
		<< "Pi = " << res << std::endl;
		std::cout << std::endl;
	}
    MPI_Finalize();
	return err;
}


/**  Exercice 4 : Produit Matriciel
 * [       chaque procs a n/p ligne de la matrice 
 * 			et calcule n éléments de y 
 * 			On Scatter La matrice a tous le monde et on broadcast X à tous le monde
 * 			Ensuite on gather les Yp et 0 affiches
 * 			A -> tableau de taille n*n  (chaque Ai a une taille n/p * n)
 * 			X -> vecteur de taille n
 * ]
 * */

int Exercice4(int argc, char **argv, int n){
	
	int id, nbprocs, err;
	MPI_Status status;
    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_rank(MPI_COMM_WORLD , &id);
    err = MPI_Comm_size(MPI_COMM_WORLD , &nbprocs);
    //~ Version c++
    //~ std::auto_ptr <float> A(new float[n*n]);
    //~ std::auto_ptr <float> A_local(new float[n*n/nbprocs]);
    //~ std::auto_ptr <float> X (new float[n]);
    //~ std::auto_ptr <float> Y (new float[n]);
    //~ std::auto_ptr <float> Y_local (new float[n/nbprocs]);
		
	
	float* A = (float*) malloc (n*n*sizeof(float));
    float* A_local = (float*) malloc (n*n/nbprocs*sizeof(float));
    float* X = (float*) malloc (n*sizeof(float));
    float* Y = (float*) malloc (n*sizeof(float));
    float* Y_local = (float*) malloc (n/nbprocs*sizeof(float));
    if(id == 0){
		for (int i=0; i<n;i++){
			X[i] = 1.0;
			Y[i] = 0.0;
			for(int j=0;j<n;j++)
				A[i*n+j] = 1.0;
		}
		
		for (int i = 0; i < n ; i++){
			for(int j = 0; j < n; j++)
				cout << " " << A[i] << " ";
			cout << endl;
		}
		
		Mat2vecProduct(A, X, n, n, Y);
		cout << nbprocs << endl;
	}
    MPI_Finalize();
	return err;
}


void Mat2vecProduct(float *A, float *X, int n, int m, float *Y){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			Y[j] += A[i*n+j] * X[j];
		}
	}
	for(int i = 0; i < n; i++){
		cout << " " << Y[i] << " ";
	}
	cout << endl;
}
