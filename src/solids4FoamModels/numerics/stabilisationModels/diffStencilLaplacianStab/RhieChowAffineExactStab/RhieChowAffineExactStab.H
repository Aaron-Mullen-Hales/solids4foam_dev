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
    RhieChowAffineExactStab

Description
    Affine-exact high-pass vector stabilisation.

    The explicit residual contribution is face based.  For an internal face
    between owner P and neighbour N, the model forms the displacement jump
    defect that cannot be represented by the supplied smooth gradient:

        D_N - D_P - (C_N - C_P) & gradD_f

    and scales it by mesh.deltaCoeffs() and scaleFactor.  This has the same
    dimensions as fvc::snGrad(D).  The inherited Laplacian vectorJacobian() is
    retained as an approximate preconditioner; it is not the exact derivative
    of this high-pass residual.

    Coupled patches use neighbour displacement, neighbour cell-centre and
    owner/neighbour gradient data from patchNeighbourField().  Non-coupled
    physical patches use the boundary displacement value minus the owner-cell
    affine extrapolation to the face centre.  The nonlinear solid model later
    applies its normal traction-boundary enforcement to the scaled face force.

SourceFiles
    RhieChowAffineExactStab.C

Author
    Aaron Mullen-Hales, UCD.

\*---------------------------------------------------------------------------*/

#ifndef RhieChowAffineExactStab_H
#define RhieChowAffineExactStab_H

#include "diffStencilLaplacianStab.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class RhieChowAffineExactStab Declaration
\*---------------------------------------------------------------------------*/

class RhieChowAffineExactStab
:
    public diffStencilLaplacianStab
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

        //- Write the affine-exact face field at normal write times
        const Switch writeAffineExactFields_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RhieChowAffineExactStab(const RhieChowAffineExactStab&);

        //- Disallow default bitwise assignment
        void operator=(const RhieChowAffineExactStab&);


public:

    //- Runtime type information
    TypeName("rhiechowAffineExact");

    // Constructors

        //- Construct from dictionary
        RhieChowAffineExactStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Member Functions

        // Access

            //- Re-calculate the vector stabilisation
            virtual void updateVector
            (
                const volVectorField& field,
                const volTensorField* gradPtr = nullptr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
