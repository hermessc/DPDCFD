//H
//procedura MPI per l'accoppiamento dei due codici in parallelo.
// 1) definizione di variabili su cui salvare l'identità dei processori, in particolare il loro rank e il size generale del world communicator.
// 2) labelling dei procs e assegnazione di un numero a tutti. Es. me = 0, è il processore master
// 3) condizione che verificare se tutti i processori in parallelo vengono condivisi o meno dai due codici
// 4) definizione di un sottoprotocollo per assegnazione degli stessi processori a Lammps.
// 5) split in questo caso non è utilizzato, ma può servire per far partire istanze multiple di lammps e dividere i processori su diverse istanze.
// 6) creazione del pointer con cui utilizzare i comandi e le funzioni di lammps.
// 7) lmparg è un array contenete i flag con cui far partire lammps, nel mio caso voglio che non stampi nulla a schermo (-screen none).
// 8) se tutto è andato bene, il pointer punterà alle librerie di lammps, sarà pronto ad eseguire quello che chiediamo su lammps e leggerà in input i 3 argomenti che gli abbiamo passato, run sul world di processori che abbiamo deciso
// 9) apertura del file (versione 1), legge il file di input dalla cartella di test. Si potrebbe aggiungere come valore nel dizionario. Apre e controlla che ci sia il file di input di lammps per sapere i comandi da eseguire. (Versione 2: si possono aggiungere due ulteriori flag a lammps "-in" e "in.supertest" che direttamente il file di input dalla cartella"

int myid;
int numProcs;
//MPI_Init(NULL,NULL);
MPI_Comm_rank(MPI_COMM_WORLD, &myid);   // rank del processore in uso, tipo 1°, 2° ecc proc.
MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // numero totale di processori usati

label nprocs = Pstream::nProcs(); // label nprocs = 4;

int me = Pstream::myProcNo();

//Info<<"numProcs "<<numProcs<<" nprocs "<<nprocs<<endl;
if (nprocs != numProcs) {
        Info<<"ERROR: problem with MPI parallel initialization of LAMMPS"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
}
else Info<<"Using all "<<nprocs<<" for LAMMPS simulation"<<endl;

MPI_Comm comm_lammps;  //set up new comm protocol for lammps
int lammps  = 1;   //check to confirm lammps is running on me
MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
LAMMPS_NS::LAMMPS *lmp;
int flagRunMe = 0;


//MPI_Bcast(&(storageMatrix[0][0]), 10*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);

