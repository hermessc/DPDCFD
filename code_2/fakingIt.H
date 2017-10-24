MPI_Barrier(MPI_COMM_WORLD);

if (me == 0)
	{
	double *array_val = NULL;
	array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	cout<< *array_val;
	double viscosity = *array_val;
 	cout<<viscosity;

		
	
		
	}

forAll(mesh.C(), celli) 
	{	//dimensionedScalar dummyNu ("dummyNu",dimensionSet(0,2,-1,0,0,0,0),1); 
		//scaling factor here
	if((shearRate[celli] < 1) && (shearRate[celli] > 0.8))
	
	{	cout<<"here";
		viscoDPD[celli] = 0.025;
	}//*dummyNu;
	if((shearRate[celli] < 0.8) && (shearRate[celli] > 0.6))
	{cout<<"hereQOQOQOQ";
		viscoDPD[celli] = 0.03;
	}
	if((shearRate[celli] < 0.4) && (shearRate[celli] > 0.6))
	{cout<<"hereAAAAA";
		viscoDPD[celli] = 0.055;
	}

	}
delete lmp;

