//H
//Falsa conversione della viscosità da Lammps a OF
// 1) Sincro dei procs
// 2) dichiarazione di un pointer che servirà a trovare il valore ultimo di viscosità ottenuto in Lammps dopo la sua run
// 3) estrazione tramite funzione di Lammps di un fix definito nel file di input di lammps. Viene castato su un double, usato extract fix per trovare dal nostro indicatore lmp, il valore chiamato "vave" (viscosity average, che ho definito io così sul file di input), gli altri numeri servono ad indicare ad esempio se si estrare un valore singolo o un vettore, da verificare gli altri flags
// 4) manipolazione fake della matrice di storage per vedere se fosse funzionante
// 5) flag per effettuare simulazioni di lammps resettata
// 6) pulizia della memoria e del pointer che si riferisce a lammps.

//cout<<"MI TROVO IN FAKINGIT.H";
MPI_Barrier(MPI_COMM_WORLD);
int a;
double *array_val = NULL;
array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	//cout<< *array_val;
double viscosity = *array_val;
 	//cout<<viscosity;
MPI_Status status;

if (me == 0) {
	for(a=0;a<5;a++) {
		if (storageMatrix[a][0] == 0) {
		storageMatrix[a][0] = 1.0;
		flagRunMe =0;		
		}
	}
}
else {
	for(a=0;a<5;a++) {
		if (storageMatrix[a][0] == 0) {
		storageMatrix[a][0] = 2.0;
		flagRunMe = 1;
		}
	}
}

/*
if (TimeCounter == 10){

for(int b=0;b<10;b++)
for(int c= 0; c<2;c++)
	storageMatrix[b][c]=0;


if(me == 0) {
	//cout<<"Dentro il processore 00000000000000";
	for( a = 0; a < 5; a++) {
		if (storageMatrix[a][0] == 0) {
		//cout<<"Dentro if per la "<<a<<" volta";
		storageMatrix[a][0] = 1.0;
       			//MPI_Bcast(storageMatrix,10*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
		//break;		
		}
	
	//	MPI_Send(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	//   	MPI_Recv(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
	
	
	}
//MPI_Bcast(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

if (me == 1) {
 	 flagRunMe = 1;
	 MPI_Bcast(&flagRunMe,1,MPI_INT,1,MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);
cout<<"\nflagRunMe "<<flagRunMe<<"from proc "<<me<<"\n";


else {
for (int c = 0; c<10;c++)
		for(int d = 0; d<2; d++) {
		MPI_Recv(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
   		MPI_Send(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	//for(int d = 5; d < 10; d++) {

	//	if (storageMatrix[d][0] == 0) {
	//	storageMatrix[d][0] = 2.0; 
		
//MPI_Recv(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);		
		//break;		
	//	}
	//}
}

//MPI_Barrier(MPI_COMM_WORLD);
//cout<<"HERE";

if(me==1)
{
	MPI_Recv(&(storageMatrix[1][0]), 10*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Send(&(storageMatrix[1][0]), 10*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	cout<<"Proc 0, Storage val: "<<storageMatrix[0][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[1][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[2][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[3][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[4][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[5][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[6][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[7][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[8][0];
	cout<<"\nProc 0, Storage val: "<<storageMatrix[9][0]<<" \n";
       // MPI_Bcast(storageMatrix,10*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

//MPI_Barrier(MPI_COMM_WORLD);

} */

//storageMatrix[0][0] = 1.0;
//cout<<storageMatrix[0][0];
//cout<<flagRunMe;


//flagRunMe = 0;
delete lmp;

