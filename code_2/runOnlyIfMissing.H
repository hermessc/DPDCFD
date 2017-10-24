
int storageMatrix[10][2] = {{0}};
int flagRunMe = 0;
flagRunMe = 0;
for (int i=0;i<10;i++)
	{
	if(std::round(shearRate[celli]) == storageMatrix[i][0])
		{
		 viscoDPD[celli] = storageMatrix[i][1];
		 flagRunMe = 1;		
		} 		
cout<<"round equal to: "<<std::round(shearRate[celli]);
	}

