//convert shear cell into shear DPD
// max val = 0.1, min val 0.01 DPD
// poiseuille flow max shear rate = 20

//H 
/* convertire shear ottenuto dalle celle cfd in valori da passare a lammps
1) creare una variabile shearDPD double e in questo caso è assegnata una legge casuale che converta questo numero in un valore fittizio dpd
2) definizione di stringhe di controllo
3) costruzione della stringa che viene usata come comando di input al codice di lammps per indicare lo shear da applicare sul sistema
*/

double shearDPD = shear_current *  0.1 / 20;

char shearStr[1024];
char controlShear[1024];
std::sprintf(controlShear,"fix		shear all deform 1 xy erate %f remap v flip yes units box\n",shearDPD);

//cout<<controlShear;
