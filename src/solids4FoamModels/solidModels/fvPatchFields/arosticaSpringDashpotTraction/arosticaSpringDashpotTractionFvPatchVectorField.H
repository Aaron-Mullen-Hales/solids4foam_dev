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
    arosticaSpringDashpotTractionFvPatchVectorField

Description
    Shared stateless implementation for the Aróstica spring-dashpot
    traction boundary conditions.

\*---------------------------------------------------------------------------*/

#ifndef arosticaSpringDashpotTractionFvPatchVectorField_H
#define arosticaSpringDashpotTractionFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class arosticaSpringDashpotTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class arosticaSpringDashpotTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

protected:

    // Protected data

        //- Spring coefficient with dimensions pressure/length
        const dimensionedScalar springCoefficient_;

        //- Dashpot coefficient with dimensions pressure*time/length
        const dimensionedScalar dashpotCoefficient_;

        //- Force reference-area traction, required by the benchmark
        const Switch forceReferenceArea_;

        //- Optional diagnostic reporting
        const Switch writeDiagnostics_;


    // Protected Member Functions

        //- Construct a base dictionary with zero static traction/pressure
        static dictionary withTractionDefaults(const dictionary&);

        //- Read a dimensioned coefficient or legacy scalar alpha/beta entry
        static dimensionedScalar readCoefficient
        (
            const dictionary&,
            const word& dimensionedName,
            const word& legacyName,
            const dimensionSet&
        );

        //- Require reference-area integration for these conditions
        static Switch requireReferenceArea(const dictionary&);

        //- Calculate spring and dashpot contributions from current trial D
        //  and Ddot
        virtual void calculateContributions
        (
            vectorField& spring,
            vectorField& dashpot,
            const vectorField&,
            const vectorField&
        ) const = 0;

        //- Report optional trial contributions without changing state
        void reportDiagnostics
        (
            const vectorField&,
            const vectorField&,
            const vectorField&,
            const vectorField&,
            const vectorField&
        ) const;


public:

    // Constructors

        //- Construct from patch and internal field
        arosticaSpringDashpotTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        arosticaSpringDashpotTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        arosticaSpringDashpotTractionFvPatchVectorField
        (
            const arosticaSpringDashpotTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        arosticaSpringDashpotTractionFvPatchVectorField
        (
            const arosticaSpringDashpotTractionFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );
#endif

        //- Construct as copy setting internal field reference
        arosticaSpringDashpotTractionFvPatchVectorField
        (
            const arosticaSpringDashpotTractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );


    // Member functions

        // Access

            //- The Aróstica conditions always use reference-area traction
            virtual const Switch& useUndeformedArea() const
            {
                return forceReferenceArea_;
            }

            //- Spring coefficient used by the compact preconditioner
            const dimensionedScalar& springCoefficient() const
            {
                return springCoefficient_;
            }

            //- Dashpot coefficient used by the compact preconditioner
            const dimensionedScalar& dashpotCoefficient() const
            {
                return dashpotCoefficient_;
            }

            //- Re-evaluate the stateless trial contributions for diagnostics
            //  without changing the stored patch traction.
            void diagnosticContributions
            (
                vectorField& spring,
                vectorField& dashpot,
                const vectorField& displacement,
                const vectorField& velocity
            ) const
            {
                calculateContributions
                (
                    spring,
                    dashpot,
                    displacement,
                    velocity
                );
            }

        //- Update the stateless trial traction
        virtual void updateCoeffs();

        //- Write the common dictionary entries
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
