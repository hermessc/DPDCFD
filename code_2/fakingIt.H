MPI_Barrier(MPI_COMM_WORLD);

if (me == 0)
	{
	double *array_val = NULL;
	array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	cout<< *array_val;
	double viscosity = *array_val;
 	cout<<viscosity;
	}
	
	if (std::round(shearRate[celli]) == 1)
		{
		viscoDPD[celli] = 0.5;
		storageMatrix[0][0] = 1.0;
		storageMatrix[0][1] = 0.5;		
		}
else
{
	double *array_val = NULL;
	array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	cout<< *array_val;
	double viscosity = *array_val;
 	cout<<viscosity;
	}
	
	if (std::round(shearRate[celli]) == 1)
		{
		viscoDPD[celli] = 0.5;
		storageMatrix[0][0] = 1.0;
		storageMatrix[0][1] = 0.5;		
		}
delete lmp;

