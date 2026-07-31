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

#include "arosticaNormalSpringDashpotTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam
{

arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF)
{}


arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(p, iF, dict)
{
    validateTangentialTraction(dict);
}


arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, p, iF, mapper)
{}


#ifndef OPENFOAM_ORG
arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf)
{}
#endif


arosticaNormalSpringDashpotTractionFvPatchVectorField::
arosticaNormalSpringDashpotTractionFvPatchVectorField
(
    const arosticaNormalSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    arosticaSpringDashpotTractionFvPatchVectorField(pvf, iF)
{}


void arosticaNormalSpringDashpotTractionFvPatchVectorField::
validateTangentialTraction(const dictionary& dict) const
{
    if (!dict.found("tangentialTraction"))
    {
        return;
    }

    const vectorField tangential
    (
        "tangentialTraction",
        dict,
        patch().size()
    );

    if (max(mag(tangential)) > SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "The Aróstica epicardial condition has zero tangential "
            << "traction; tangentialTraction must be zero"
            << exit(FatalIOError);
    }
}


void arosticaNormalSpringDashpotTractionFvPatchVectorField::
calculateContributions
(
    vectorField& spring,
    vectorField& dashpot,
    const vectorField& displacement,
    const vectorField& velocity
) const
{
    const vectorField N(patch().nf());
    const scalarField normalDisplacement(displacement & N);
    const scalarField normalVelocity(velocity & N);

    spring =
        -springCoefficient_.value()*normalDisplacement*N;

    dashpot =
        -dashpotCoefficient_.value()*normalVelocity*N;
}


makePatchTypeField
(
    fvPatchVectorField,
    arosticaNormalSpringDashpotTractionFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
