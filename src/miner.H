MPI_Barrier(MPI_COMM_WORLD);
double *array_val = NULL;
array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);

double viscosity = *array_val;


for (int i = 0; i< 100; i++){
	if (shear_current == bigMM[i][0]){
	bigMM[i][1] = viscosity ;
	}
}



delete lmp;
}
