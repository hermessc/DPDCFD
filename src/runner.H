






 	double  shear_current = 0.0;	
for (int i=0; i < non_zeros ; i++) {	

	MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
	shear_current = bigMM[i][0];
	#include "shearConverter.H"
	#include "lammpsRun.H"	
		
		//MPI_Barrier(MPI_COMM_WORLD);
		


