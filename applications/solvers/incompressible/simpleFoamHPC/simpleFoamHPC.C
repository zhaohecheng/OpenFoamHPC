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
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
    Foam::Info<< "Create time\n" << Foam::endl;

    Foam::Time runTime(Foam::Time::controlDictName, args);
    Info << "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(args, runTime));

    dynamicFvMesh& mesh = meshPtr();

    simpleControl simple(mesh);
    Info << "Reading field p\n" << endl;
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
        fvc::flux(U)
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    IOMRFZoneList MRF(mesh);

    fv::options& fvOptions(fv::options::New(mesh));

    if (!fvOptions.optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }
#include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.controlledUpdate();

        if (mesh.changing())
        {
            MRF.update();
        }

        // --- Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            MRF.correctBoundaryVelocity(U);

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
                + MRF.DDt(U)
                + turbulence->divDevReff(U)
                ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            if (simple.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));

                fvOptions.correct(U);
            }
            {
                volScalarField rAU(1.0 / UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                MRF.makeRelative(phiHbyA);
                adjustPhi(phiHbyA, U, p);

                tmp<volScalarField> rAtU(rAU);

                if (simple.consistent())
                {
                    rAtU = 1.0 / (1.0 / rAU - UEqn.H1());
                    phiHbyA +=
                        fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf();
                    HbyA -= (rAU - rAtU()) * fvc::grad(p);
                }

                tUEqn.clear();

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAtU(), MRF);

                // Non-orthogonal pressure corrector loop
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve();

                    if (simple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

#include "continuityErrs.H"

                // Explicitly relax pressure for momentum corrector
                p.relax();

                // Momentum corrector
                U = HbyA - rAtU() * fvc::grad(p);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //