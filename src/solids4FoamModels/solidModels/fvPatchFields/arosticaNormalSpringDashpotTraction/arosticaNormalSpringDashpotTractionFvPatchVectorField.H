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

Class
    arosticaNormalSpringDashpotTractionFvPatchVectorField

Description
    Aróstica epicardial normal spring-dashpot traction.

\*---------------------------------------------------------------------------*/

#ifndef arosticaNormalSpringDashpotTractionFvPatchVectorField_H
#define arosticaNormalSpringDashpotTractionFvPatchVectorField_H

#include "arosticaSpringDashpotTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class arosticaNormalSpringDashpotTractionFvPatchVectorField
\*---------------------------------------------------------------------------*/

class arosticaNormalSpringDashpotTractionFvPatchVectorField
:
    public arosticaSpringDashpotTractionFvPatchVectorField
{

    // Private Member Functions

        //- Calculate normal-only spring-dashpot traction
        virtual void calculateContributions
        (
            vectorField& spring,
            vectorField& dashpot,
            const vectorField&,
            const vectorField&
        ) const;

        //- Validate the optional legacy zero tangential entry
        void validateTangentialTraction(const dictionary&) const;


public:

    //- Runtime type information
    TypeName("arosticaNormalSpringDashpotTraction");


    // Constructors

        arosticaNormalSpringDashpotTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        arosticaNormalSpringDashpotTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        arosticaNormalSpringDashpotTractionFvPatchVectorField
        (
            const arosticaNormalSpringDashpotTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        arosticaNormalSpringDashpotTractionFvPatchVectorField
        (
            const arosticaNormalSpringDashpotTractionFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );
#endif

        arosticaNormalSpringDashpotTractionFvPatchVectorField
        (
            const arosticaNormalSpringDashpotTractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new arosticaNormalSpringDashpotTractionFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting the internal field
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new arosticaNormalSpringDashpotTractionFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
