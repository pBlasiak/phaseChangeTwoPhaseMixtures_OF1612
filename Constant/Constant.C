/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Constant.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Constant, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Constant::Constant
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    mCondFlux_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("condMassFlux")),
    mEvapFlux_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("evapMassFlux")),
    vmCondFlux_
    (
        IOobject
        (
            "vmCondFlux",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("vmCondFlux", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    vmEvapFlux_
    (
        IOobject
        (
            "vmEvapFlux",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("vmEvapFlux", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    )
{
	Info<< "Constant model settings:  " << endl;
	Info<< "Condensation mass flow rate per unit area: " << mCondFlux_ << endl;
	Info<< "Evaporation mass flow rate per unit area: "  << mEvapFlux_ << endl;
	vmCondFlux_ = mCondFlux_*calcGradAlphal();
	vmEvapFlux_ = mEvapFlux_*calcGradAlphal();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Constant::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotAlphal() const
{

	return Pair<tmp<volScalarField> >
	(
		vmCondFlux_*scalar(1),
	   -vmEvapFlux_*scalar(1) 
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotP() const
{
	return Pair<tmp<volScalarField> >
	(
	//	vmCondFlux_*pos(p_ - pSat_)/max(p_-pSat_,1E-6*pSat_),
	//   -vmEvapFlux_*neg(p_ - pSat_)/max(pSat_-p_,1E-6*pSat_)
		vmCondFlux_*scalar(1),
	   -vmEvapFlux_*scalar(1)
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotT() const
{
	Info<< "vmEvapFlux_ = " << vmEvapFlux_ << endl;
	return Pair<tmp<volScalarField> >
	(
		-vmCondFlux_*neg(T_ - TSat_)/max(TSat_ - T_,1E-6*TSat_),
	     vmEvapFlux_*pos(T_ - TSat_)/max(T_ - TSat_,1E-6*TSat_)
		//-vmCondFlux_*scalar(1),
	    // vmEvapFlux_*scalar(1)
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::vmDot() const
{
}

bool Foam::phaseChangeTwoPhaseMixtures::Constant::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condMassFlux") >> mCondFlux_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evapMassFlux") >> mEvapFlux_;

		// to raczej nie jest dobrze
        vmCondFlux_ = mCondFlux_*calcGradAlphal();
	    vmEvapFlux_ = mEvapFlux_*calcGradAlphal();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
