//H
//Falsa conversione della viscosità da Lammps a OF
// 1) Sincro dei procs
// 2) dichiarazione di un pointer che servirà a trovare il valore ultimo di viscosità ottenuto in Lammps dopo la sua run
// 3) estrazione tramite funzione di Lammps di un fix definito nel file di input di lammps. Viene castato su un double, usato extract fix per trovare dal nostro indicatore lmp, il valore chiamato "vave" (viscosity average, che ho definito io così sul file di input), gli altri numeri servono ad indicare ad esempio se si estrare un valore singolo o un vettore, da verificare gli altri flags
// 4) manipolazione fake della matrice di storage per vedere se fosse funzionante
// 5) flag per effettuare simulazioni di lammps resettata
// 6) pulizia della memoria e del pointer che si riferisce a lammps.


MPI_Barrier(MPI_COMM_WORLD);

	double *array_val = NULL;
	array_val =  (double *) lammps_extract_fix(lmp,(char *) "vave",0,0,1,1);
	//cout<< *array_val;
	double viscosity = *array_val;
 	//cout<<viscosity;
			viscoDPD[celli] =viscosity ;
storageMatrix[0][0] = 1.0;
//cout<<storageMatrix[0][0];
//cout<<flagRunMe;
flagRunMe = 0;
delete lmp;

