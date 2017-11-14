
forAll(mesh.C(),celli) {
	for(int b = 0; b < int(rows/nSplits); b++) {
		if (bigM[b][0] == 0) {
			bigM[b][0] = std::round(shearRate[celli]);
			MPI_Bcast(&bigM[b][0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			break;			
		}
	} 
}
	
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

int total_length = 0;
MPI_Reduce(&cutter_length, &total_length, 1, MPI_INT, MPI_SUM, 0,
        	   MPI_COMM_WORLD);
	
MPI_Bcast(&total_length,1,MPI_INT,0,MPI_COMM_WORLD);

MPI_Allgather(&cutter,cutter_length,MPI_INT,&big_m_clone,25,MPI_INT,MPI_COMM_WORLD);				

if(me == 0){
	cout<<"collector: ";
	for (int j = 0; j < 50; j++) {
		cout<<big_m_clone[j]<<" ";
	}
}



