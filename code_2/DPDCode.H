shearRate = mag(fvc::grad(U)); 
MPI_Barrier(MPI_COMM_WORLD);

if (TimeCounter == DPD_Sim_Every_X_Timestep)
	{		
		//while(g!=0)
	forAll(mesh.C(), celli)	{
		#include "runOnlyIfMissing.H"
		#include "initialSetup.H"
		#include "lammpsRun.H"
		#include "fakingIt.H"
		//	g = 0;		
		//dimensionedScalar dummyNu ("dummyNu",dimensionSet(0,2,-1,0,0,0,0),1); 
		//scaling factor here
		if((shearRate[celli] < 1) && (shearRate[celli] > 0.8)) {
			viscoDPD[celli] = 0.025;
			}
		if((shearRate[celli] < 0.8) && (shearRate[celli] > 0.6)) {
			viscoDPD[celli] = 0.03;
			}
		if((shearRate[celli] < 0.4) && (shearRate[celli] > 0.6)) {
			viscoDPD[celli] = 0.055;
			}
	//	}			
	}}
