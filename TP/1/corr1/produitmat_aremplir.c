/*****************************************************************/
/*                                                               */
/* Produit de matrice par double diffusion                       */
/*                                                               */
/*****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>

/******* Fonctions d'affichage ***********/
#define VOIRD(expr) do {printf("P%d (%d,%d) <%.3d> : \t{ " #expr " = %d }\n", \
my_id,i_row,i_col,__LINE__,expr);} while(0)
#define PRINT_MESSAGE(expr) do {printf("P%d (%d,%d) <%.3d> : \t" #expr "\n", \
my_id,i_row,i_col,__LINE__);} while(0)

typedef struct { 
	int row;   /* le nombre de lignes   */
	int col;   /* le nombre de colonnes */
	int* data; /* les valeurs           */
} matrix_t;

int my_id, nb_proc;
int nb_col,nb_row;       /* Largeur et hauteur de la grille */
MPI_Comm MPI_HORIZONTAL; /* Communicateur pour diffusions horizontales */
MPI_Comm MPI_VERTICAL;   /* Communicateur pour diffusions verticales */
int i_col,i_row;         /* Position dans la grille */


/*******************************************/
/* initialisation aléatoire de la matrice  */
/*******************************************/
void mat_init_alea(matrix_t *mat,int width, int height)
{
	int i;
	
	mat->row = height;
	mat->col = width;
	
	mat->data=(int*)malloc(sizeof(int)*height*width);
	for(i=0 ; i<height*width ; i++)
		mat->data[i]=rand()%2;
}

/*******************************************/
/* initialisation à 0 de la matrice        */
/*******************************************/
void mat_init_empty(matrix_t *mat,int width, int height)
{
	mat->row = height;
	mat->col = width;
	mat->data=(int*)calloc((mat->row)*(mat->col),sizeof(int));
}

/*******************************************/
/* affichage de la matrice                 */
/*******************************************/
void mat_display(matrix_t A)
{
	int i,j,t=0;
	
	printf("      ");
	for(j=0;j<A.col;j++)
		printf("%.3d ",j);
	printf("\n");
	printf("    __");
	for(j=0;j<A.col;j++)
		printf("____");
	printf("_\n");
	
	for(i=0;i<A.row;i++)
    {
		printf("%.3d | ",i);
		for(j=0;j<A.col;j++)
		/*printf("%d ",A.data[i*A.row+j]);*/
			printf("%.3d ",A.data[t++]);
		printf("|\n");
    }
	printf("    --");
	for(j=0;j<A.col;j++)
		printf("----");
	printf("-\n");
}

/*******************************************/
/* C+=A*B                                  */
/*******************************************/
void mat_mult(matrix_t A, matrix_t B, matrix_t C)
{
	int i,j,k,M,N,K,*_A,*_B,*_C;
	
	_A=A.data;
	_B=B.data;
	_C=C.data;
	
	M = C.row;
	N = C.col;
	K = A.col;
	
	if((M!=A.row) || (N!=B.col) || (K!=B.row)) {
		PRINT_MESSAGE("Attention, tailles incompatibles");
		VOIRD(A.row);VOIRD(A.col);VOIRD(B.row);
		VOIRD(C.col);VOIRD(C.row);VOIRD(C.col);
		exit(1);
	}
	
	for(i=0 ; i<M ; i++)
		for(j=0 ; j<N ; j++)
			for(k=0 ; k<K ; k++) 
				_C[i*N+j]+=_A[i*K+k]*_B[k*N+j];
}

/*******************************************/
/* Initialisations la grille               */
/*******************************************/
void init_communicateurs()
{
	nb_col=nb_proc; nb_row=1;
	/* Améliorez la disposition */
	/*A REMPLIR */
	
	/* Partitionnez MPI_COMM_WORLD en communicateurs pour les lignes 
     et en communicateurs pour les colonnes */
     /*A REMPLIR */
	
	/* Vérifiez que votre rang dans le nouveau communicateur 
     correspond bien à celui que vous attendez. */
     /*A REMPLIR */
}

/*******************************************/
/* Produit de matrice par double diffusion */
/*******************************************/
void parallel_mat_mult(matrix_t A, matrix_t B, matrix_t C)
{
	int k;  
	
	for(k=0; k<(nb_col*C.col) ; k++) {
		/* diffusion horizontale */
		/*A REMPLIR */
		/* diffusion verticale */
		/*A REMPLIR */
		/* produit de matrices en local */
		/*A REMPLIR */
	}
}


int main(int argc, char **argv) 
{
	int taille;
	
	/* initialisation mpi */
	
	/*A REMPLIR */
	
	if(argc==0)
		fprintf(stderr,"Usage : %s <taille de matrice>",argv[0]);
	else 
		taille = atoi(argv[1]);
	
	init_communicateurs();
	
	{
		matrix_t A,B,C;
		
		mat_init_alea(&A,taille/nb_col,taille/nb_row);
		mat_init_alea(&B,taille/nb_col,taille/nb_row);
		mat_init_empty(&C,taille/nb_col,taille/nb_row);
		parallel_mat_mult(A,B,C);
		if(my_id==0) mat_display(C);
	}
	MPI_Finalize();
	
	return (0);
}
