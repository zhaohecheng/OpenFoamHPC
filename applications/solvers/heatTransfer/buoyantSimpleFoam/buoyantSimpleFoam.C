/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    buoyantSimpleFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent flow of compressible fluids,
    including radiation, for ventilation and heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Steady-state solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation."
    );

#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
    Foam::Info << "Create time\n" << Foam::endl;

    Foam::Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
#include "createControl.H"
    Info << "Reading thermophysical properties\n" << endl;

    autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
    rhoThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    Info << "Reading field U\n" << endl;
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

    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );

    Info << "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );


    Info<< "\nReading g" << endl;
    const meshObjects::gravity& g = meshObjects::gravity::New(runTime);
    Info<< "\nReading hRef" << endl;
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(word::null, dimLength, Zero)
    );
    Info<< "Calculating field g.h\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );
    const int oldLocal = volVectorField::Boundary::localConsistency;
    volVectorField::Boundary::localConsistency = 0;

    volScalarField gh("gh", (g & mesh.C()) - ghRef);
    surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

    volVectorField::Boundary::localConsistency = oldLocal;


    Info << "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Force p_rgh to be consistent with p
    p_rgh = p - rho * gh;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        simple.dict(),
        pRefCell,
        pRefValue
    );

    mesh.setFluxRequired(p_rgh.name());

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
    dimensionedScalar totalVolume = sum(mesh.V());

    IOMRFZoneList MRF(mesh);
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(thermo.T())
    );

    const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, simple.dict());
    const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, simple.dict());

    fv::options& fvOptions(fv::options::New(mesh));

    if (!fvOptions.optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }

    const volScalarField& psi = thermo.psi();
#include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
            // Solve the Momentum equation

            MRF.correctBoundaryVelocity(U);

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + MRF.DDt(rho, U)
              + turbulence->divDevRhoReff(U)
             ==
                fvOptions(rho, U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            if (simple.momentumPredictor())
            {
                solve
                (
                    UEqn
                 ==
                    fvc::reconstruct
                    (
                        (
                          - ghf*fvc::snGrad(rho)
                          - fvc::snGrad(p_rgh)
                        )*mesh.magSf()
                    )
                );

                fvOptions.correct(U);
            }

#include "EEqn.H"
#include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
