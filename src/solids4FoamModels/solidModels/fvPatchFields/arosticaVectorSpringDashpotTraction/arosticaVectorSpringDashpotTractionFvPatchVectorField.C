/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "arosticaVectorSpringDashpotTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF)
{}


arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF, dict)
{}


arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, p, iF, mapper)
{}


#ifndef OPENFOAM_ORG
arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf)
{}
#endif


arosticaVectorSpringDashpotTractionFvPatchVectorField::
arosticaVectorSpringDashpotTractionFvPatchVectorField
(
    const arosticaVectorSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, iF)
{}


void arosticaVectorSpringDashpotTractionFvPatchVectorField::
calculateContributions
(
    vectorField& spring,
    vectorField& dashpot,
    const vectorField& displacement,
    const vectorField& velocity
) const
{
    spring = -springCoefficient_.value()*displacement;
    dashpot = -dashpotCoefficient_.value()*velocity;
}


makePatchTypeField
(
    fvPatchVectorField,
    arosticaVectorSpringDashpotTractionFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
