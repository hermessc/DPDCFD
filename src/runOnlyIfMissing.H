

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
//flagRunMe = 1;
//cout<<"FlagRunMe on proc "<<me<<"equal to "<<flagRunMe;
//MPI_Bcast(&flagRunMe,1,MPI_INT,0,MPI_COMM_WORLD);
//flagRunMe = 1;
MPI_Barrier(MPI_COMM_WORLD);


