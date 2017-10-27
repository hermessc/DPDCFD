//H
//Dizionario aggiunto ai parametri di simulazione
//Per il momento contiene soltanto il numero di CFD timestep da aspettare prima di effettuare una sim Lammps ma può essere aggiunta roba se servisse
IOdictionary DPDParameters
    (
        IOobject
        (
            "DPDParameters",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
//Some examples

scalar DPD_Sim_Every_X_Timestep      = readScalar(DPDParameters.subDict("Simulation").lookup("DPD_Sim_Every_X_Timestep"));





