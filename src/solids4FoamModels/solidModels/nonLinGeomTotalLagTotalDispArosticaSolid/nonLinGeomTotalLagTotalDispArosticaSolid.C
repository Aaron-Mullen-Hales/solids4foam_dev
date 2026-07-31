/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nonLinGeomTotalLagTotalDispArosticaSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispArosticaSolid, 0);
addToRunTimeSelectionTable
(
    solidModel,
    nonLinGeomTotalLagTotalDispArosticaSolid,
    dictionary
);


void nonLinGeomTotalLagTotalDispArosticaSolid::
completeDirectFaceCauchyStress
(
    surfaceSymmTensorField& sigmaFace
) const
{
    sigmaFace -=
        fvc::interpolate(p())
       *dimensionedSymmTensor("I", dimless, symmTensor::I);
}


nonLinGeomTotalLagTotalDispArosticaSolid::
nonLinGeomTotalLagTotalDispArosticaSolid
(
    Time& runTime,
    const word& region
)
:
    nonLinGeomTotalLagTotalDispSolid(typeName, runTime, region)
{
    if (!solvePressure())
    {
        FatalErrorInFunction
            << typeName << " requires solvePressure true" << nl
            << exit(FatalError);
    }

    Info<< nl
        << "Aróstica mixed stress split:" << nl
        << "    sigma = sigmaNonVolumetricFull - p*I" << nl
        << "    pressure is owned by the mixed solid model" << endl;
}


} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
