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
    RhieChowZeroForceStab

Description
    Rhie-Chow stabilisation with a projection that removes the global
    translational resultant from the extensive cell-vector stabilisation force.

SourceFiles
    RhieChowZeroForceStab.C

Author
    Aaron Mullen-Hales, UCD.

\*---------------------------------------------------------------------------*/

#ifndef RhieChowZeroForceStab_H
#define RhieChowZeroForceStab_H

#include "RhieChowStab.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RhieChowZeroForceStab Declaration
\*---------------------------------------------------------------------------*/

class RhieChowZeroForceStab
:
    public RhieChowStab
{
    // Private data

        //- Report zero-force projection diagnostics
        const Switch reportZeroForceDiagnostics_;

        //- Absolute tolerance for the post-projection resultant force
        const scalar zeroForceAbsTolerance_;

        //- Relative tolerance for the post-projection resultant force
        const scalar zeroForceRelTolerance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RhieChowZeroForceStab(const RhieChowZeroForceStab&);

        //- Disallow default bitwise assignment
        void operator=(const RhieChowZeroForceStab&);


public:

    //- Runtime type information
    TypeName("rhiechowZeroForce");

    // Constructors

        //- Construct from dictionary
        RhieChowZeroForceStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Member Functions

        //- Project an extensive cell-vector stabilisation force
        virtual void projectExtensiveVectorForce
        (
            vectorField& stabilisationForce
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
