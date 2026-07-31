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

\*---------------------------------------------------------------------------*/

#include "RhieChowAffineExactStab.H"
#include "addToRunTimeSelectionTable.H"
#include "compatibilityFunctions.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RhieChowAffineExactStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, RhieChowAffineExactStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RhieChowAffineExactStab::RhieChowAffineExactStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    diffStencilLaplacianStab(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
    writeAffineExactFields_
    (
        dict.lookupOrDefault<Switch>("writeAffineExactFields", false)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RhieChowAffineExactStab::updateVector
(
    const volVectorField& field,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();

    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << field.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    const volTensorField& gradField = *gradPtr;

    if (faceVectorPtr().empty())
    {
        faceVectorPtr().set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "faceStabilisation(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    writeAffineExactFields_
                  ? IOobject::AUTO_WRITE
                  : IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("0", dims(), vector::zero)
            )
        );
    }

    surfaceVectorField& faceStab = autoPtrRef(faceVectorPtr());

    // The affine defect is built from owner-to-neighbour differences and
    // therefore changes sign if the face orientation is reversed.
    faceStab.setOriented();

    const surfaceTensorField gradF(fvc::interpolate(gradField));

    const vectorField& fieldI = primitiveField(field);
    const vectorField& C = mesh().C();
    const tensorField& gradFI = primitiveField(gradF);
    const scalarField& deltaCoeffs = mesh().deltaCoeffs();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    vectorField& faceStabI = primitiveFieldRef(faceStab);

    forAll(faceStabI, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector affineDefect =
            fieldI[nei]
          - fieldI[own]
          - ((C[nei] - C[own]) & gradFI[faceI]);

        faceStabI[faceI] = scaleFactor_*deltaCoeffs[faceI]*affineDefect;
    }

    surfaceVectorField::Boundary& faceStabB = faceStab.boundaryFieldRef();

    forAll(faceStabB, patchI)
    {
        const fvPatch& patch = mesh().boundary()[patchI];
        const labelUList& faceCells = patch.faceCells();
        const scalarField& patchDeltaCoeffs =
            mesh().deltaCoeffs().boundaryField()[patchI];

        vectorField& patchFaceStab = faceStabB[patchI];
        const tmp<vectorField> tFieldOwn =
            field.boundaryField()[patchI].patchInternalField();
        const vectorField& fieldOwn = tFieldOwn();

        if (field.boundaryField()[patchI].coupled())
        {
            const tmp<vectorField> tFieldNei =
                field.boundaryField()[patchI].patchNeighbourField();
            const vectorField& fieldNei = tFieldNei();

            const tmp<tensorField> tGradOwn =
                gradField.boundaryField()[patchI].patchInternalField();
            const tensorField& gradOwn = tGradOwn();

            const tmp<tensorField> tGradNei =
                gradField.boundaryField()[patchI].patchNeighbourField();
            const tensorField& gradNei = tGradNei();

            const tmp<vectorField> tCNei =
                mesh().C().boundaryField()[patchI].patchNeighbourField();
            const vectorField& CNei = tCNei();

            const tmp<vectorField> tCOwn =
                mesh().C().boundaryField()[patchI].patchInternalField();
            const vectorField& COwn = tCOwn();

            const scalarField& patchWeights = patch.weights();

            forAll(patchFaceStab, faceI)
            {
                const tensor gradFace =
                    patchWeights[faceI]*gradOwn[faceI]
                  + (1.0 - patchWeights[faceI])*gradNei[faceI];

                const vector affineDefect =
                    fieldNei[faceI]
                  - fieldOwn[faceI]
                  - ((CNei[faceI] - COwn[faceI]) & gradFace);

                patchFaceStab[faceI] =
                    scaleFactor_*patchDeltaCoeffs[faceI]*affineDefect;
            }
        }
        else
        {
            const vectorField& fieldB = field.boundaryField()[patchI];
            const vectorField& Cf = mesh().Cf().boundaryField()[patchI];
            const tensorField& gradFieldI = primitiveField(gradField);

            forAll(patchFaceStab, faceI)
            {
                const label own = faceCells[faceI];

                const vector affineDefect =
                    fieldB[faceI]
                  - fieldOwn[faceI]
                  - ((Cf[faceI] - C[own]) & gradFieldI[own]);

                patchFaceStab[faceI] =
                    scaleFactor_*patchDeltaCoeffs[faceI]*affineDefect;
            }
        }
    }

}


// ************************************************************************* //
