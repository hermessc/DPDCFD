//H
/*
1) ricalcolo dello shear, per sicurezza.
2) header che converte il valore di shear CFD in un valore fittizio di shear DPD
3) Sincro dei procs
4) ci sono due cicli che sono stati inseriti per verificare che tutto succedesse correttamente sui diversi procs, ma in realtà il primo if è sufficiente al funzionamento del codice.
5) cout di verifica processore
6) flag derivante dal controllo dello shear in matrici di storage per evitare simulazioni ripetute
7) run di lammps finchè ci sono righe nel file di testo
8) il primo processore controlla se ci sono line di lunghezza diverse da 0 e Bcasta il valore della lunghezza di n agli altr procs (lammps utilizza entrambi i processori per una singola simulazione invece openFoam divide il dominio e i calcolo sui due procs. Quindi significa che metà celle sono risolte da proc 0 di OF e metà da proc 1 di OF se nprocs = 2. Invece lammps utilizza entrambi i procs per una singola simulazione e poi passa alla cella succesiva).
9) se la linea è vuota esci dal ciclo e Bcasta questa info agli altri procs.
10) è stata definita una stringa fittizia che controlla se siamo arrivati al punto dello script di lammps in cui troviamo il termine #shear. Se questo è verificato si invia come linea di comando una linea costruita artificialmente in shearConverter.H invece che la linea successiva del file di input
11) altrimenti si continua a leggere il file di input riga per riga. 
12) se lo shear era già presente nella matrice di storage non serve effettuare una simulazione di lammps ma si può continuare con il calcolo cfd
*/
shearRate = mag(fvc::grad(U)); 
//cout<<"Prima di lammps";
#include "shearConverter.H"

//cout<<"Dopo shearConverter";
      //  std::sprintf(shearStr,"fix		shear all deform 1 xy erate ${srate} remap v flip yes units box",controlVariable);
MPI_Barrier(MPI_COMM_WORLD);
if (me == 0){
	cout<<"Runnin Lammps from proc: "<<me<<"\n\n";
//cout<<"proc num"<<me<<"\n";
	if (flagRunMe == 0){
//		cout<<"Dentro if\n\n";
		while (1) {
			if (me == 0) {
				if (fgets(line,1024,sp) == NULL) n = 0;
				else n = strlen(line) + 1;
				if (n == 0) fclose(sp);
				}
			MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
			if (n == 0) break;
			MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
			std::string controlstr(line);
			if (controlstr == "#shear\n") {
				if (lammps == 1) lmp->input->one(controlShear);
				}
			else lmp->input->one(line);
		//if (lammps == 1) lmp->input->one(line);
		//extract viscosity
		}
				
	}
	else
	continue;	

}

else {
//cout<<"proc num"<<me<<"\n";
	if (flagRunMe == 0){
		cout<<"Runnin Lammps from proc: "<<me<<"\n\n";
		while (1) {
			if (me == 0) {
				if (fgets(line,1024,sp) == NULL) n = 0;
				else n = strlen(line) + 1;
				if (n == 0) fclose(sp);
				}
			MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
			if (n == 0) break;
			MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
			std::string controlstr(line);
			if (controlstr == "#shear\n") {
				if (lammps == 1) lmp->input->one(controlShear);
				}
			else lmp->input->one(line);
		//if (lammps == 1) lmp->input->one(line);
		//extract viscosity
			}
				
	}
	else
	continue;
}

