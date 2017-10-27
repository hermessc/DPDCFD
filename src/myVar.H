//H
/*
alcune variabili utilizzate nel codice, tra cui la lunghezza delle stringhe dello script di lammps, un char contentente la linea, una matrice di storage e un flag per verificare se eseguire o meno le simulazioni
*/
#include "string.h"
int n;
char line[1024];
int storageMatrix[10][2] = {{0}};
int flagRunMe = 0;
