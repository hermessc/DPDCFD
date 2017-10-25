//convert shear cell into shear DPD
// max val = 0.1, min val 0.01 DPD
// poiseuille flow max shear rate = 20


double shearDPD = shearRate[celli] *  0.1 / 20;

char shearStr[1024];
char controlShear[1024];
std::sprintf(controlShear,"fix		shear all deform 1 xy erate %f remap v flip yes units box",shearDPD);

cout<<controlShear;
