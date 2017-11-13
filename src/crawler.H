
		//shearRate = mag(fvc::grad(U)); 
		forAll(mesh.C(),celli) {
			for(int b = 0; b < int(rows/nSplits); b++) {
				if (bigM[b][0] == 0) {
					bigM[b][0] = std::round(shearRate[celli]);
					MPI_Bcast(&bigM[b][0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					break;			
				}
			} 
		}
	

//Step 2 rimozione valori doppi e sorting delle matrici:
		int cutter_length=0;
		for(int i = 0; i < (rows-1); i++) {
	 		cutter[cutter_length] = bigM[i][0];
			if (bigM[i][0] != bigM[i+1][0]) {
				cutter_length++;		
			}
		}

		if(me ==0){
			cout<<"cutter: ";
			for (int j = 0; j < rows; j++) {
			cout<<cutter[j];
			}
		}
	

//	cout<<"The value of l is:     "<< cutter_length<<"\n\n\n\n\n\n\n";


	int total_length = 0;
	MPI_Reduce(&cutter_length, &total_length, 1, MPI_INT, MPI_SUM, 0,
        	   MPI_COMM_WORLD);
	
	MPI_Bcast(&total_length,1,MPI_INT,0,MPI_COMM_WORLD);
//float *big_m_clone = NULL;
  //big_m_clone = malloc(sizeof(float) * cutter[rows]);

int	big_m_clone[100] = {0};
//da sostituire 25 con il numero di celle per processore
	MPI_Allgather(&cutter,cutter_length,MPI_INT,&big_m_clone,25,MPI_INT,MPI_COMM_WORLD);				
	//MPI_Barrier(MPI_COMM_WORLD);
	

 		//   printf("Process 1 received number %d from process 0\n",
           //number);
		
//	cout<<"The value of total length is   "<<total_length<<"\n\n\n\n\n";

	

//step 3 creare una matrice unica che raccoglie i valori degli shear

/*
float big_m_clone[rows][2] ={{0.0}};



	for ( int j = (me*cells_per_proc); j< ((me+1)*cells_per_proc); j++) {
		if (me == 0) {
 MPI_Recv(&(bigM[j][0]), 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD,
 	            MPI_STATUS_IGNORE);
			big_m_clone[j][0] = bigM[j][0];
			MPI_Bcast(&(big_m_clone[j][0]),1,MPI_DOUBLE,0, MPI_COMM_WORLD);  
		}
 		else if (me != 0) {
 		     			MPI_Send(&(bigM[j][0]), 1, MPI_FLOAT, me, 0, MPI_COMM_WORLD);
 		//   printf("Process 1 received number %d from process 0\n",
           //number);
		}

*/
		/*
			int cutter_adjusted_position = (me * nCells / nSplits) + j;
			big_m_clone[j][0] = cutter[cutter_adjusted_position];
			MPI_Bcast(&big_m_clone[j][0], 1, MPI_DOUBLE, me, MPI_COMM_WORLD);
		//int cutter_adjusted_position = (me * nCells / nSplits) + j;			
		//MPI_Recv(&cutter[cutter_adjusted_position],1,MPI_FLOAT,me,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
		//big_m_clone[j][0] = cutter[cutter_adjusted_position];
		//MPI_Bcast(&big_m_clone[j][0], 1, MPI_FLOAT,0,MPI_COMM_WORLD);
		//sleep(1);

		cout<<"\nme equals "<<me;
			cout<<"\n j equals "<< j;
			cout<<"\n matrix equals "<< bigM[j][0];
			cout<<"\n new matrix equals "<< big_m_clone[j][0];		
*/
		


	if(me == 0){
			cout<<"collector: ";
			for (int j = 0; j < 50; j++) {
			cout<<big_m_clone[j]<<" ";
			}
		}



