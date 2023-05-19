/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "solidificationSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidificationSystem, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationSystem::solidificationSystem
(
    const volVectorField& U
)
:
    Foam::twoPhaseMixture::IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh()),

    rho1_
    (
        "rho",
        dimensionSet(1, -3, 0, 0, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    rho2_
    (
        "rho",
        dimensionSet(1, -3, 0, 0, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),

    Cp1_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    Cp2_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),

    kappa1_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    kappa2_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),

    mu1_
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    mu2_
    (
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),

    D1_
    (
        "D",
        dimensionSet(0, 2, -1, 0, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    D2_
    (
        "D",
        dimensionSet(0, 2, -1, 0, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),
    DAS_
    (
        "DAS",
        dimensionSet(0, 1, 0, 0, 0),
        twoPhaseMixture::subDict(this->phase1Name())
    ),
    betaT_
    (
        "betaT",
        dimensionSet(0, 0, 0, -1, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),
    betaC_
    (
        "betaC",
        dimensionSet(0, 0, 0, 0, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),
    TRef_
    (
        "TRef",
        dimensionSet(0, 0, 0, 1, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),
    CRef_
    (
        "CRef",
        dimensionSet(0, 0, 0, 0, 0),
        twoPhaseMixture::subDict(this->phase2Name())
    ),

    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidificationSystem::read()
{
    const dictionary& phaseDict1(twoPhaseMixture::subDict(this->phase1Name()));
    const dictionary& phaseDict2(twoPhaseMixture::subDict(this->phase2Name()));

    rho1_ = phaseDict1.lookup("rho");
    rho2_ = phaseDict2.lookup("rho");

    Cp1_ = phaseDict1.lookup("Cp");
    Cp2_ = phaseDict2.lookup("Cp");

    kappa1_ = phaseDict1.lookup("kappa");
    kappa2_ = phaseDict2.lookup("kappa");

    mu1_ = phaseDict1.lookup("mu");
    mu2_ = phaseDict2.lookup("mu");

    D1_ = phaseDict1.lookup("D");
    D2_ = phaseDict2.lookup("D");

    DAS_ = phaseDict1.lookup("DAS");

    betaT_ = phaseDict2.lookup("betaT");
    betaC_ = phaseDict2.lookup("betaC");

    TRef_ = phaseDict2.lookup("TRef");
    CRef_ = phaseDict2.lookup("CRef");

    return true;
}


// ************************************************************************* //
