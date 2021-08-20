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

    mCondFlux_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("condMassFlux")),
    mEvapFlux_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("evapMassFlux"))
{
	Info<< "Constant model settings:  " << endl;
	Info<< "Condensation mass flow rate per unit area: " << mCondFlux_ << endl;
	Info<< "Evaporation mass flow rate per unit area: "  << mEvapFlux_ << endl;
	
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	volScalarField gradAlphal = mag(fvc::grad(limitedAlpha1));

	// minus sign "-" to provide mc > 0  and mv < 0
    mCondNoAlphal_ =  mCondFlux_*gradAlphal*pos(1-limitedAlpha1)/max(1-limitedAlpha1, 1.0);
    mEvapNoAlphal_ = -mEvapFlux_*gradAlphal*pos(limitedAlpha1)/max(limitedAlpha1, 1.0);

	mCondAlphal_   = mCondNoAlphal_;//*(1-limitedAlpha1);
	mEvapAlphal_   = mEvapNoAlphal_;//*limitedAlpha1;

	// minus sign to provide mc > 0  and mv > 0
	mCondNoTmTSat_ =  mCondAlphal_/T1;
	mEvapNoTmTSat_ = -mEvapAlphal_/T1;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Constant::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotAlphal() 
{
	return Pair<tmp<volScalarField> >
	(
		mCondNoAlphal_*scalar(1),
	    mEvapNoAlphal_*scalar(1) 
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotP() const
{
	return Pair<tmp<volScalarField> >
	(
		mCondAlphal_*scalar(1)*pos(p_ - pSat_)/max(p_-pSat_,1E-6*pSat_),
	    mEvapAlphal_*scalar(1)//*neg(p_ - pSat_)/max(pSat_-p_,1E-6*pSat_)
//    	vmCondFlux_*scalar(1),
//       -vmEvapFlux_*scalar(1)
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Constant::mDotT() const
{
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
	return Pair<tmp<volScalarField> >
	(
		 mCondNoTmTSat_*neg(T_ - TSat_)/max(TSat_ - T_,1E-6*TSat_)*T1,
	     mEvapNoTmTSat_*pos(T_ - TSat_)/max(T_ - TSat_,1E-6*TSat_)*T1
		//-vmCondFlux_*scalar(1),
	    // vmEvapFlux_*scalar(1)
	);
}

bool Foam::phaseChangeTwoPhaseMixtures::Constant::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condMassFlux") >> mCondFlux_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evapMassFlux") >> mEvapFlux_;

		// to raczej nie jest dobrze
        //vmCondFlux_ = mCondFlux_*calcGradAlphal();
	    //vmEvapFlux_ = mEvapFlux_*calcGradAlphal();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
