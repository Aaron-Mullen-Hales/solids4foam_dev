/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    nonLinGeomTotalLagTotalDispGultekinSolid

Description
    Isolated finite-volume mixed displacement-pressure model for the original
    Gultekin-Dal-Holzapfel Q1P0+WAS continuum equations. This is not a finite
    element implementation. It retains the complete passive Cauchy stress and
    uses the finite-K constraint associated with

        U(J) = K*(J - ln(J) - 1).

    With the solids4foam pressure sign convention:

        sigma = sigmaPassive - p*I

        R_p = -p/K - (J - 1)/J + S_p.

    Pressure and momentum stabilisation remain numerical finite-volume terms.

SourceFiles
    nonLinGeomTotalLagTotalDispGultekinSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomTotalLagTotalDispGultekinSolid_H
#define nonLinGeomTotalLagTotalDispGultekinSolid_H

#include "nonLinGeomTotalLagTotalDispSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

class nonLinGeomTotalLagTotalDispGultekinSolid
:
    public nonLinGeomTotalLagTotalDispSolid
{
    // Private data

        //- Write optional formulation-decomposition fields
        const Switch writeGultekinMixedDiagnostics_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        nonLinGeomTotalLagTotalDispGultekinSolid
        (
            const nonLinGeomTotalLagTotalDispGultekinSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomTotalLagTotalDispGultekinSolid&);


protected:

    // Protected Member Functions

        //- Retain the full matrix-plus-unsplit-fibre passive stress
        virtual bool retainFullPassiveStressInMixedSplit() const
        {
            return true;
        }

        //- g(J) = (J - 1)/J = 1 - 1/J
        virtual tmp<volScalarField> mixedVolumetricConstraint
        (
            const volScalarField& J
        ) const;

        //- Scalar g(J)
        virtual scalar mixedVolumetricConstraintValue(const scalar J) const;

        //- dg/dJ = 1/J^2
        virtual scalar mixedVolumetricConstraintDerivative
        (
            const scalar J
        ) const;

        //- Human-readable g(J) expression
        virtual string mixedVolumetricConstraintDescription() const;


public:

    //- Runtime type information
    TypeName("nonLinGeomTotalLagTotalDispGultekinSolid");

    // Constructors

        nonLinGeomTotalLagTotalDispGultekinSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    //- Destructor
    virtual ~nonLinGeomTotalLagTotalDispGultekinSolid()
    {}

    //- Write standard fields and optional Gultekin mixed diagnostics
    virtual void writeFields(const Time& runTime);
};

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //

#endif

