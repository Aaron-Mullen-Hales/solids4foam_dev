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
    RhieChowHighPassStab

Description
    Experimental high-pass Rhie-Chow vector stabilisation.

    The model first applies a cached, geometry-based moving least-squares
    polynomial smoothing operator to the cell displacement field,

        Dhigh = D - P_smooth(D)

    and then applies the standard Rhie-Chow/difference-stencil face operator to
    Dhigh only. This differs from ordinary Rhie-Chow, which acts on the
    complete displacement field, and from the affine-defect experiment, which
    directly rewrites the neighbour jump. For any polynomial field reproduced
    by the cached local reconstruction, Dhigh is zero and no momentum
    stabilisation is produced.

    The inherited Laplacian vectorJacobian() is retained as an approximate
    preconditioner; it is not the exact derivative of RC((I - P_smooth)D).

SourceFiles
    RhieChowHighPassStab.C

Author
    Aaron Mullen-Hales, UCD.

\*---------------------------------------------------------------------------*/

#ifndef RhieChowHighPassStab_H
#define RhieChowHighPassStab_H

#include "diffStencilLaplacianStab.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToCellStencil;

/*---------------------------------------------------------------------------*\
                   Class RhieChowHighPassStab Declaration
\*---------------------------------------------------------------------------*/

class RhieChowHighPassStab
:
    public diffStencilLaplacianStab
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

        //- Requested polynomial order for the smooth reconstruction
        const label highPassPolynomialOrder_;

        //- Requested topological stencil layers
        const label highPassStencilLayers_;

        //- Exponent in the distance-based least-squares weight
        const scalar highPassWeightExponent_;

        //- Report high-pass diagnostics
        const Switch reportHighPassDiagnostics_;

        //- Report interval in updateVector calls
        const label highPassReportInterval_;

        //- Write high-pass diagnostic fields
        const Switch writeHighPassFields_;

        //- Relative singular-value cutoff for local SVD fits
        const scalar svdMinCondition_;

        //- Active geometric directions used in the polynomial basis
        labelList activeDirections_;

        //- Parallel distribution map for the compacted cell-only stencil
        autoPtr<mapDistribute> mapPtr_;

        //- Compacted, cell-only smoothing stencil
        labelListList stencil_;

        //- Cached smoothing weights corresponding to stencil_
        List<scalarList> filterWeights_;

        //- Number of cells that used a lower polynomial order
        label nOrderFallbacks_;

        //- Minimum polynomial order used by any local filter
        label minOrderUsed_;

        //- Maximum polynomial reproduction error in the cached filters
        scalar maxReproductionError_;

        //- Cached smooth displacement from the latest update
        mutable autoPtr<volVectorField> smoothFieldPtr_;

        //- Cached high-pass displacement from the latest update
        mutable autoPtr<volVectorField> highFieldPtr_;

        //- Count updateVector calls for throttled diagnostics
        mutable label updateCalls_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RhieChowHighPassStab(const RhieChowHighPassStab&);

        //- Disallow default bitwise assignment
        void operator=(const RhieChowHighPassStab&);

        //- Return the number of polynomial basis entries for the order
        label nBasis(const label order) const;

        //- Evaluate the dimensionless polynomial basis
        scalarList basis(const vector& xi, const label order) const;

        //- Build the smoothing stencil and cached weights
        void buildFilter();

        //- Build from OpenFOAM cell-to-cell stencils
        void buildFilter(const cellToCellStencil& rawStencil);

        //- Build from OpenFOAM cell-to-cell stencils
        void buildFilter
        (
            const cellToCellStencil& rawStencil,
            const cellToCellStencil* supplementalStencil
        );

        //- Apply the cached smoothing filter to a vector field
        void applyFilter
        (
            const volVectorField& field,
            volVectorField& smoothField
        ) const;

        //- Ensure the diagnostic/working fields exist
        void ensureWorkingFields(const volVectorField& field) const;

        //- Zero all non-coupled physical boundary patch values
        void zeroPhysicalBoundaries(volVectorField& field) const;

        //- Zero non-coupled physical boundary face stabilisation
        void zeroPhysicalBoundaryFaces(surfaceVectorField& field) const;

        //- Optionally report filter/update diagnostics
        void reportDiagnostics
        (
            const volVectorField& field,
            const volVectorField& smoothField,
            const volVectorField& highField,
            const surfaceVectorField& faceStab
        ) const;


public:

    //- Runtime type information
    TypeName("rhiechowHighPass");

    // Constructors

        //- Construct from dictionary
        RhieChowHighPassStab
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

            //- Smooth displacement from the latest update
            const volVectorField& smoothField() const;

            //- High-pass displacement from the latest update
            const volVectorField& highField() const;

            //- Number of rank/order fallbacks in the cached filter
            label nOrderFallbacks() const
            {
                return nOrderFallbacks_;
            }

            //- Maximum polynomial reproduction error in the cached filter
            scalar maxReproductionError() const
            {
                return maxReproductionError_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
