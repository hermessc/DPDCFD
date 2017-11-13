//aggiungere ciclo for per ridurre le ripetizioni

//nome matrice big_m_clone

int bigMM[100][2] = {0};
int new_counter = 0;
for(int i = 0; i < 100 ; i ++ )
{
	if(big_m_clone[i] != 0) {
		bigMM[new_counter][0] = big_m_clone[i];
		for (int j= 0; j< 100; j++) {
			if(bigMM[new_counter][0] == big_m_clone[j]) {
			big_m_clone[j] = 0;			
			}

		}		new_counter++;	
	}
}

if(me == 0){
			cout<<"collector: ";
			for (int j = 0; j < 100; j++) {
			cout<<bigMM[j][0]<<" ";
			}
		}
