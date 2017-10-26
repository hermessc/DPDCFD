shearRate = mag(fvc::grad(U)); 
//cout<<"Prima di lammps";
#include "shearConverter.H"

//cout<<"Dopo shearConverter";
      //  std::sprintf(shearStr,"fix		shear all deform 1 xy erate ${srate} remap v flip yes units box",controlVariable);
MPI_Barrier(MPI_COMM_WORLD);
if (me == 0){
//cout<<"proc num"<<me<<"\n";
	if (flagRunMe == 0){
//		cout<<"Dentro if\n\n";
		while (1) {
			MPI_Barrier(MPI_COMM_WORLD);
			if (me == 0) {
				if (fgets(line,1024,sp) == NULL) n = 0;
				else n = strlen(line) + 1;
				if (n == 0) fclose(sp);
				}
			MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
			if (n == 0) break;
			MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
			std::string controlstr(line);
			MPI_Barrier(MPI_COMM_WORLD);
			if (controlstr == "#shear\n") {
				if (lammps == 1) lmp->input->one(controlShear);
				}
			else lmp->input->one(line);
			MPI_Barrier(MPI_COMM_WORLD);
		//if (lammps == 1) lmp->input->one(line);
		//extract viscosity
		}
				
	}
	else
	continue;	
	MPI_Barrier(MPI_COMM_WORLD);
}

else {
//cout<<"proc num"<<me<<"\n";
	if (flagRunMe == 0){
//		cout<<"Dentro if\n\n";
		while (1) {
			MPI_Barrier(MPI_COMM_WORLD);
			if (me == 0) {
				if (fgets(line,1024,sp) == NULL) n = 0;
				else n = strlen(line) + 1;
				if (n == 0) fclose(sp);
				}
			MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
			if (n == 0) break;
			MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
			std::string controlstr(line);
			MPI_Barrier(MPI_COMM_WORLD);
			if (controlstr == "#shear\n") {
				if (lammps == 1) lmp->input->one(controlShear);
				}
			else lmp->input->one(line);
			MPI_Barrier(MPI_COMM_WORLD);
		//if (lammps == 1) lmp->input->one(line);
		//extract viscosity
			}
				
	}
	else
	continue;
	MPI_Barrier(MPI_COMM_WORLD);	
}
