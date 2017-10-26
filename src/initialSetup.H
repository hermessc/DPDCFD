int myid;
int numProcs;
//MPI_Init(NULL,NULL);
MPI_Comm_rank(MPI_COMM_WORLD, &myid);   // rank del processore in uso, tipo 1°, 2° ecc proc.
MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // numero totale di processori usati

label nprocs = Pstream::nProcs(); // label nprocs = 4;

int me = Pstream::myProcNo();

//Info<<"numProcs "<<numProcs<<" nprocs "<<nprocs<<endl;
if (nprocs != numProcs) {
        Info<<"ERROR: problem with MPI parallel initialization of LAMMPS"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
}
else Info<<"Using all "<<nprocs<<" for LAMMPS simulation"<<endl;

MPI_Comm comm_lammps;  //set up new comm protocol for lammps
int lammps  = 1;   //check to confirm lammps is running on me
MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
LAMMPS_NS::LAMMPS *lmp;
char **lmparg = new char*[3];
  lmparg[0] = NULL;        
  lmparg[1] =   (char *) "-screen";
  lmparg[2] =   (char *) "none";
if (lammps == 1) {lmp = new LAMMPS_NS::LAMMPS(3,lmparg,comm_lammps);}


 sp = fopen("/home/hermes/Desktop/gitDPD/DPDCFD/poiseuilleflow/in.supertest","r");
    		if (sp == NULL)
		{
			printf("ERROR: Could not open LAMMPS input script\n");
      			MPI_Abort(MPI_COMM_WORLD,1);
    		}
