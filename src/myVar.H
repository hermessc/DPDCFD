//H
/*
alcune variabili utilizzate nel codice, tra cui la lunghezza delle stringhe dello script di lammps, un char contentente la linea, una matrice di storage e un flag per verificare se eseguire o meno le simulazioni
*/
#include "string.h"
int n;
char line[1024];

//int flagRunMe = 0;
/*double **storageMatrix;
//storageMatrix = (double **)malloc(20*sizeof(double));

storageMatrix = (double **)malloc(20*sizeof(double*));
MPI_Bcast(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
*/
int **storageMatrix;
//storageMatrix = (double **)malloc(20*sizeof(double));

storageMatrix = alloc_2d_int(10,2);
//cout<<"STO QUA";
//MPI_Bcast(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
for(int b=0;b<10;b++)
for(int c= 0; c<2;c++)
	storageMatrix[b][c]=0;
