shearRate = mag(fvc::grad(U)); 
MPI_Barrier(MPI_COMM_WORLD);

if (TimeCounter == DPD_Sim_Every_X_Timestep)
	{		
		//while(g!=0)
	forAll(mesh.C(), celli)	{
		#include "initialSetup.H"
		#include "runOnlyIfMissing.H"
		#include "lammpsRun.H"
		#include "fakingIt.H"
		//	g = 0;		
		//dimensionedScalar dummyNu ("dummyNu",dimensionSet(0,2,-1,0,0,0,0),1); 
		//scaling factor here
	
	//	}			
	}}
