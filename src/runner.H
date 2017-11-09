for (int i=0; i< cutter_length; i++) { 
		//cout<<"HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERMES";
		MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
		#include "lammpsRun.H"	
		
		//MPI_Barrier(MPI_COMM_WORLD);
		}
