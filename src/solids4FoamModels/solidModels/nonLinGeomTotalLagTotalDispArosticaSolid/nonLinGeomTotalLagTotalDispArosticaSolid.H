/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    nonLinGeomTotalLagTotalDispArosticaSolid

Description
    Total-Lagrangian mixed displacement-pressure solid model for the
    Aróstica benchmark. It reuses the complete base implementation and
    retains the full non-volumetric constitutive Cauchy stress before adding
    the independent mixed pressure contribution.

SourceFiles
    nonLinGeomTotalLagTotalDispArosticaSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomTotalLagTotalDispArosticaSolid_H
#define nonLinGeomTotalLagTotalDispArosticaSolid_H

#include "nonLinGeomTotalLagTotalDispSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
              Class nonLinGeomTotalLagTotalDispArosticaSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomTotalLagTotalDispArosticaSolid
:
    public nonLinGeomTotalLagTotalDispSolid
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        nonLinGeomTotalLagTotalDispArosticaSolid
        (
            const nonLinGeomTotalLagTotalDispArosticaSolid&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const nonLinGeomTotalLagTotalDispArosticaSolid&
        );


protected:

    // Protected Member Functions

        //- Retain the complete non-volumetric constitutive stress.
        virtual bool retainFullPassiveStressInMixedSplit() const
        {
            return true;
        }

        //- Insert the physical mixed pressure in direct face stress.
        virtual void completeDirectFaceCauchyStress
        (
            surfaceSymmTensorField& sigmaFace
        ) const;


public:

    //- Runtime type information
    TypeName("arosticaNonLinearGeometryTotalLagrangianTotalDisplacement");


    // Constructors

        //- Construct from components
        nonLinGeomTotalLagTotalDispArosticaSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~nonLinGeomTotalLagTotalDispArosticaSolid()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
