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





