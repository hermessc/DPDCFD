//H
/*
alcune variabili utilizzate nel codice, tra cui la lunghezza delle stringhe dello script di lammps, un char contentente la linea, una matrice di storage e un flag per verificare se eseguire o meno le simulazioni
*/
#include "string.h"
int n;
char line[1024];
int storageMatrix[10][2] = {{0}};
int flagRunMe = 0;
char checkStr[1024] ="#fix		shear all deform 1 xy erate ${srate} remap v flip yes units box\n";
int rows = nCells;

cout<<"rows:: "<< rows;
double bigM[rows][2] = {{0.0}};
double cutter[rows] = {0.0};

