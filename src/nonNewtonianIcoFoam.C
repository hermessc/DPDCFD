/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    nonNewtonianIcoFoam

Description
    Transient solver for incompressible, laminar flow of non-Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//H
// includere mpi.h per parallelo
#include "mpi.h"
#include <string>
#include "singlePhaseTransportModel.H"
#include "pisoControl.H"
//H
//includere gli headers di lammps per funzioni
#include "lammpsHeaders.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshNoClear.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "DPDictionary.H"
	//H
	//definire le mie variabili e dichiarare il dizionario DPD 
    #include "myVar.H"
   double bigM[rows][2] = {{0.0}};
    #include "initialSetup.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    int TimeCounter = 0;

    while (runTime.loop())
    {

        
	Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
	fluid.correct();
	shearRate = mag(fvc::grad(U)); 

	forAll(mesh.C(),celli) {	
	#include "shearConverter.H"	
	cout<<"I am the shear "<<shearRate[celli]<<" from proc "<< me<<"\n";
	}


/*
for(int i = 0; i < rows; i++) {
		for(int j = 0; j < 2 ; j++ ) {
		cout<<" "<<bigM[i][j]<<" ";
		}
		cout<<"\n";
	}*/
//	MPI_Barrier(MPI_COMM_WORLD);

//Step 1: Creazione prima matrice di accumulo shear //
	if (TimeCounter == DPD_Sim_Every_X_Timestep) {		
		//shearRate = mag(fvc::grad(U)); 
		forAll(mesh.C(),celli) {
			for(int b = 0; b < int(rows/nSplits); b++) {
				if (bigM[b][0] == 0) {
				//MPI_Bcast(&bigM[j][0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				bigM[b][0] = std::round(shearRate[celli]);
				MPI_Bcast(&bigM[b][0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				break;			
			}
		} 
		
		}
	}

//Step 2 rimozione valori doppi e sorting delle matrici:


/* 

Insert Here

*/ 


//
	//MPI_Barrier(MPI_COMM_WORLD);
	
	
//CheckPoint: stampa valorI
	for(int i = 0; i < int(rows/nSplits); i++) {
		for(int j = 0; j < 2 ; j++ ) {
		cout<<" "<<bigM[i][j]<<" ";
		}
		cout<<"\n";
	}


//Step 4: Running lammps
	/* REMOVE COMMENTS TO ENABLE LAMMPS	
	for (int i=0; i< 5; i++) { 
	MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
	#include "lammpsRun.H"	
		
		//MPI_Barrier(MPI_COMM_WORLD);
		}
	*/









//H
// parte che contiene tutto quello che serve per far partire le simulazioni di lammps prima del calcolo del campo di velocità
	//FILE *sp; // open LAMMPS input script
   //if (me == 0) 
  	//{
    
	//}
	//DPD simulation
	//if (me == 0)
        //{
	//#include "DPDCode.H"
	//}
	//else
	//{
	//#include "DPDCode.H"
	//}
	//#include "viscosityMod.H"
        // Momentum predictor
// fine

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(viscoDPD, U)
          - (fvc::grad(U) & fvc::grad(viscoDPD))
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
	       
	}
	//H
	//variabile che conta il numero di timestep che sono passati nel calcolo cfd
	TimeCounter++; 	
	cout<<TimeCounter<<"\n\n\n\n";        
	runTime.write();
	//delete[] bigM;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
	
    return 0;
}


// ************************************************************************* //
