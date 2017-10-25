MPI_Barrier(MPI_COMM_WORLD);

	double *array_val = NULL;
	array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	cout<< *array_val;
	double viscosity = *array_val;
 	cout<<viscosity;
		
flagRunMe = 0;
delete lmp;

