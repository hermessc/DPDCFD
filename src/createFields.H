Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//H
//Prima variabile di testing usata su icofoam (fluidi newtoniani) per provare a inserire una viscosità fittizia
dimensionedScalar DPDnu (
	"DPDnu",
	dimViscosity,
	0.1
);
//variabile per cambio dimensionale
dimensionedScalar dummyNu (
	"dummyNu",
	dimViscosity,
	1
);
//field contenente le viscosità ricavate da DPD per OGNI CELLA del dominio computazionale
Info<< "Hey! I am faking the viscosity \n" << endl;
    volScalarField viscoDPD
    (
        IOobject
        (
	    "viscoDPD",
	     runTime.timeName(),
	     mesh,
	     IOobject::MUST_READ,
	     IOobject::AUTO_WRITE
        ),
         mesh     
    );

volScalarField shearRate
    (
	IOobject
	(
	"shearRate",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
	mesh	
    ); 
#include "createPhi.H"


singlePhaseTransportModel fluid(U, phi);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());





//Matrici


