/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::thermalConductivityModels::Linear

Description
TODO poprawic
    Linear phase change model. It is based on saturation temperature.
	For T < TSat => condensation occurs
	and T > TSat => evaporation occurs

    Reference:
    \verbatim
        Linear, W.H.
        "A Pressure Iteration Scheme for Two-Phase Modeling. Technical Report LA-UR 79-975",
        Los Alamos Scientific Laboratory, Los Alamos, New Mexico, 1979.
    \endverbatim

SourceFiles
    Linear.C

\*--------------------------------------------------------------------*/

#ifndef Linear_H
#define Linear_H

#include "thermalConductivity.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalConductivityModels
{

/*--------------------------------------------------------------------*\
                              Class Linear
\*--------------------------------------------------------------------*/

class Linear
:
    public thermalConductivity
{
    // Private data

        //dimensionedScalar Cc_;
        //dimensionedScalar Cv_;

        //dimensionedScalar mcCoeff_;
        //dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("Linear");


    // Constructors

        //- construct from components
        Linear
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Linear()
    {}


    // Member Functions
	
        //- Return the effective laminar thermal conductivity
        virtual tmp<volScalarField> k(const thermalIncompressibleTwoPhaseMixture* titpm) const;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalConductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
