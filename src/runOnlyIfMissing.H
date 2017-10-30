//H
/*
creazione di una matrice di storage che conserva i valori di shear già trovati e la viscosità calcolata con lammps. In questo modo si vuole ridurre il calcolo di viscosità per valori di shear già presenti nel sistema.

1) inizializzazione del flag
2) la matrice è grande 10, ciclo sulle righe per vedere se il valore di shear è stato già trovato in precedenza
3) se positivo, il flag cambia valore e la viscosità della cella viene semplicemente letta
4) proc sincro
*/
flagRunMe = 0;
for (int i=0;i<10;i++)
	{
	if(std::round(shearRate[celli]) == storageMatrix[i][0])
		{
		 viscoDPD[celli] = storageMatrix[i][1];
		 flagRunMe = 1;	
		//MPI_Bcast(&flagRunMe,1,MPI_INT,0,MPI_COMM_WORLD);
		//cout << "sto sabotando tutto "<<storageMatrix[0][0];	
		} 		
		//cout << "sto sabotando tutto "<<storageMatrix[0][0];	
	MPI_Barrier(MPI_COMM_WORLD);
	}

for(int b=0;b<10;b++){
	for(int c= 0; c<2;c++) {
		cout<<storageMatrix[b][c]<<" ";
	}
cout<<"\n";	
}

//flagRunMe = 1;
//cout<<"FlagRunMe on proc "<<me<<"equal to "<<flagRunMe;
//MPI_Bcast(&flagRunMe,1,MPI_INT,0,MPI_COMM_WORLD);
//flagRunMe = 1;
//MPI_Barrier(MPI_COMM_WORLD);


