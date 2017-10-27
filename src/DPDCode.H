//H
//File che contiene il ciclo di operazioni per far funzionare le simulazioni di LAMMPS
// 1) calcolo dello shear cella per cella
// 2) MPI barrier per sincronizzare tutti i singoli processori.
// 3) TimeCounter è una variabile definita fuori da questo file che controlla quanti timesteps CFD sono trascorsi.
// 4) Verificare che timecounter sia uguale all'intervallo che noi poniamo nel dizionario DPDictionary. Significa che le simulazioni Lammps sono effettuate ogni $(DPD_Sim_Every_X_Timestep).
// 5) ciclo su ogni cella del dominio CFD, viene lanciato lammps
// 6) vedi i singoli headers.
// 7) Quando la simulazione CFD è terminata, TimeCounter torna a 0 (viene incrementato alla fine di ogni timestep CFD)
// 8) Sincro dei procs 


shearRate = mag(fvc::grad(U)); 
MPI_Barrier(MPI_COMM_WORLD);



if (TimeCounter == DPD_Sim_Every_X_Timestep) {		
	//while(g!=0)
	forAll(mesh.C(), celli)	{
		#include "initialSetup.H"
		//#include "runOnlyIfMissing.H"
		#include "lammpsRun.H"
		//MPI_Bcast(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#include "fakingIt.H"


		//	g = 0;		
		//dimensionedScalar dummyNu ("dummyNu",dimensionSet(0,2,-1,0,0,0,0),1); 
		//scaling factor here
	
	//	}			
	}

TimeCounter = 0;
MPI_Barrier(MPI_COMM_WORLD);
}

	//MPI_Barrier(MPI_COMM_WORLD);
