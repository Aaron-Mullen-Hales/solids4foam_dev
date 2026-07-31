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

#include "RhieChowHighPassStab.H"
#include "addToRunTimeSelectionTable.H"
#include "CFCCellToCellStencil.H"
#include "CPCCellToCellStencil.H"
#include "SVD.H"
#include "compatibilityFunctions.H"
#include "extendedCellToFaceStencil.H"
#include "fvc.H"
#include "globalIndex.H"
#include "volFields.H"

#include <cmath>

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RhieChowHighPassStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, RhieChowHighPassStab, stabModel
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::label Foam::RhieChowHighPassStab::nBasis(const label order) const
{
    const label nDims = activeDirections_.size();

    label n = 1;

    if (order >= 1)
    {
        n += nDims;
    }

    if (order >= 2)
    {
        n += nDims;
        n += (nDims*(nDims - 1))/2;
    }

    return n;
}


Foam::scalarList Foam::RhieChowHighPassStab::basis
(
    const vector& xi,
    const label order
) const
{
    scalarList values(nBasis(order), 0.0);

    label n = 0;
    values[n++] = 1.0;

    if (order >= 1)
    {
        forAll(activeDirections_, dirI)
        {
            values[n++] = xi[activeDirections_[dirI]];
        }
    }

    if (order >= 2)
    {
        forAll(activeDirections_, dirI)
        {
            const scalar x = xi[activeDirections_[dirI]];
            values[n++] = sqr(x);
        }

        forAll(activeDirections_, dirI)
        {
            const scalar x = xi[activeDirections_[dirI]];

            for (label dirJ = dirI + 1; dirJ < activeDirections_.size(); dirJ++)
            {
                const scalar y = xi[activeDirections_[dirJ]];
                values[n++] = x*y;
            }
        }
    }

    return values;
}


void Foam::RhieChowHighPassStab::buildFilter()
{
    const vectorField& C = mesh().C();

    vector minC(GREAT, GREAT, GREAT);
    vector maxC(-GREAT, -GREAT, -GREAT);

    forAll(C, cellI)
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            minC[cmpt] = min(minC[cmpt], C[cellI][cmpt]);
            maxC[cmpt] = max(maxC[cmpt], C[cellI][cmpt]);
        }
    }

    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        reduce(minC[cmpt], minOp<scalar>());
        reduce(maxC[cmpt], maxOp<scalar>());
    }

    const vector span(maxC - minC);
    const scalar maxSpan = max(max(span.x(), span.y()), span.z());
    DynamicList<label> activeDirs(vector::nComponents);

    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        if (span[cmpt] > max(SMALL, 1e-8*maxSpan))
        {
            activeDirs.append(cmpt);
        }
    }

    activeDirections_.transfer(activeDirs);

    if (activeDirections_.empty())
    {
        activeDirections_.setSize(1);
        activeDirections_[0] = 0;
    }

    CFCCellToCellStencil rawStencil(mesh());
    CPCCellToCellStencil supplementalStencil(mesh());

    buildFilter(rawStencil, &supplementalStencil);
}


void Foam::RhieChowHighPassStab::buildFilter
(
    const cellToCellStencil& rawStencil
)
{
    buildFilter(rawStencil, nullptr);
}


void Foam::RhieChowHighPassStab::buildFilter
(
    const cellToCellStencil& rawStencil,
    const cellToCellStencil* supplementalStencil
)
{
    const globalIndex& numbering = rawStencil.globalNumbering();
    const labelList nCellsPerProc(UPstream::allGatherValues(mesh().nCells()));

    labelListList oneLayerStencil(rawStencil.size());

    forAll(rawStencil, cellI)
    {
        const labelList& raw = rawStencil[cellI];
        DynamicList<label> cells(raw.size());

        forAll(raw, i)
        {
            const label globalI = raw[i];
            const label proci = numbering.whichProcID(globalI);
            const label localI = numbering.toLocal(proci, globalI);

            if (localI < nCellsPerProc[proci] && !cells.found(globalI))
            {
                cells.append(globalI);
            }
        }

        if (cells.empty())
        {
            cells.append(numbering.toGlobal(cellI));
        }

        oneLayerStencil[cellI].transfer(cells);
    }

    labelListList supplementalOneLayer;

    if (supplementalStencil != nullptr)
    {
        supplementalOneLayer.setSize(supplementalStencil->size());

        forAll((*supplementalStencil), cellI)
        {
            const labelList& raw = (*supplementalStencil)[cellI];
            DynamicList<label> cells(raw.size());

            forAll(raw, i)
            {
                const label globalI = raw[i];
                const label proci = numbering.whichProcID(globalI);
                const label localI = numbering.toLocal(proci, globalI);

                if (localI < nCellsPerProc[proci] && !cells.found(globalI))
                {
                    cells.append(globalI);
                }
            }

            supplementalOneLayer[cellI].transfer(cells);
        }
    }

    labelListList globalStencil(oneLayerStencil.size());
    const label effectiveStencilLayers =
        max(highPassStencilLayers_, highPassPolynomialOrder_ >= 2 ? 2 : 1);

    List<labelHashSet> visited(oneLayerStencil.size());
    labelListList frontier(oneLayerStencil.size());

    forAll(oneLayerStencil, cellI)
    {
        const label ownGlobal = numbering.toGlobal(cellI);

        visited[cellI].insert(ownGlobal);
        frontier[cellI].setSize(1, ownGlobal);
    }

    for (label layerI = 0; layerI < effectiveStencilLayers; layerI++)
    {
        labelListList compactFrontier(frontier);
        List<Map<label>> compactMap(Pstream::nProcs());
        mapDistribute frontierMap(numbering, compactFrontier, compactMap);

        List<labelList> nbrData(frontierMap.constructSize());

        forAll(oneLayerStencil, localCellI)
        {
            nbrData[localCellI] = oneLayerStencil[localCellI];
        }

        // Exchange one-layer neighbour lists so layer expansion does not stop
        // at processor boundaries.
        frontierMap.distribute(nbrData);

        forAll(frontier, cellI)
        {
            DynamicList<label> nextFrontier;
            const labelList& compactCells = compactFrontier[cellI];

            forAll(compactCells, compactI)
            {
                const label compactCellI = compactCells[compactI];
                const labelList& nbrs = nbrData[compactCellI];

                forAll(nbrs, nbrI)
                {
                    const label nbrGlobal = nbrs[nbrI];

                    if (!visited[cellI].found(nbrGlobal))
                    {
                        visited[cellI].insert(nbrGlobal);
                        nextFrontier.append(nbrGlobal);
                    }
                }
            }

            frontier[cellI].transfer(nextFrontier);
        }
    }

    forAll(oneLayerStencil, cellI)
    {
        if (supplementalStencil != nullptr)
        {
            const labelList& supplemental = supplementalOneLayer[cellI];

            forAll(supplemental, suppI)
            {
                visited[cellI].insert(supplemental[suppI]);
            }
        }

        globalStencil[cellI] = visited[cellI].sortedToc();
    }

    stencil_.transfer(globalStencil);

    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset(new mapDistribute(numbering, stencil_, compactMap));

    List<List<vector>> stencilC;
    extendedCellToFaceStencil::collectData
    (
        mapPtr_(),
        stencil_,
        mesh().C(),
        stencilC
    );

    filterWeights_.setSize(stencil_.size());
    nOrderFallbacks_ = 0;
    minOrderUsed_ = highPassPolynomialOrder_;
    maxReproductionError_ = 0.0;

    const vectorField& C = mesh().C();
    const scalarField& V = mesh().V();

    forAll(stencil_, cellI)
    {
        const List<vector>& localC = stencilC[cellI];

        scalar h = 0.0;

        forAll(localC, i)
        {
            h = max(h, mag(localC[i] - C[cellI]));
        }

        h = max(h, Foam::pow(max(V[cellI], VSMALL), 1.0/3.0));

        label orderUsed = min(max(highPassPolynomialOrder_, label(0)), label(2));
        scalarList weights(localC.size(), 0.0);
        bool accepted = false;

        while (!accepted && orderUsed >= 0)
        {
            const label nCoeffs = nBasis(orderUsed);

            if (localC.size() >= nCoeffs)
            {
                scalarRectangularMatrix M(localC.size(), nCoeffs, 0.0);
                scalarList sqrtW(localC.size(), 0.0);

                forAll(localC, i)
                {
                    const vector xi((localC[i] - C[cellI])/h);
                    const scalar r = mag(xi);
                    const scalar w =
                        1.0/(1.0 + Foam::pow(r, highPassWeightExponent_));

                    sqrtW[i] = Foam::sqrt(max(w, VSMALL));

                    const scalarList phi(basis(xi, orderUsed));

                    forAll(phi, coeffI)
                    {
                        M(i, coeffI) = sqrtW[i]*phi[coeffI];
                    }
                }

                SVD svd(M, svdMinCondition_);

                if (svd.converged() && svd.nZeros() == 0)
                {
                    const scalarRectangularMatrix pinv(svd.VSinvUt());

                    forAll(weights, i)
                    {
                        weights[i] = pinv(0, i)*sqrtW[i];
                    }

                    scalar maxLocalError = 0.0;

                    for (label coeffI = 0; coeffI < nCoeffs; coeffI++)
                    {
                        scalar reproduced = 0.0;

                        forAll(localC, i)
                        {
                            const vector xi((localC[i] - C[cellI])/h);
                            reproduced += weights[i]*basis(xi, orderUsed)[coeffI];
                        }

                        const scalar exact = coeffI == 0 ? 1.0 : 0.0;
                        maxLocalError =
                            max(maxLocalError, mag(reproduced - exact));
                    }

                    maxReproductionError_ =
                        max(maxReproductionError_, maxLocalError);

                    accepted = true;
                }
            }

            if (!accepted)
            {
                orderUsed--;
            }
        }

        if (!accepted)
        {
            FatalErrorInFunction
                << "Unable to construct even a constant high-pass filter for "
                << "cell " << cellI << abort(FatalError);
        }

        if (orderUsed < highPassPolynomialOrder_)
        {
            nOrderFallbacks_++;
        }

        minOrderUsed_ = min(minOrderUsed_, orderUsed);
        filterWeights_[cellI].transfer(weights);
    }

    reduce(nOrderFallbacks_, sumOp<label>());
    reduce(minOrderUsed_, minOp<label>());
    reduce(maxReproductionError_, maxOp<scalar>());

    if (reportHighPassDiagnostics_)
    {
        Info<< "rhiechowHighPass filter setup" << nl
            << "    requested polynomial order = "
            << highPassPolynomialOrder_ << nl
            << "    minimum order used = " << minOrderUsed_ << nl
            << "    requested stencil layers = "
            << highPassStencilLayers_ << nl
            << "    active geometric directions = " << activeDirections_ << nl
            << "    rank/order fallbacks = " << nOrderFallbacks_ << nl
            << "    max polynomial reproduction error = "
            << maxReproductionError_ << endl;
    }
}


void Foam::RhieChowHighPassStab::applyFilter
(
    const volVectorField& field,
    volVectorField& smoothField
) const
{
    List<List<vector>> stencilValues;
    extendedCellToFaceStencil::collectData
    (
        mapPtr_(),
        stencil_,
        field,
        stencilValues
    );

    vectorField& smoothI = primitiveFieldRef(smoothField);

    forAll(smoothI, cellI)
    {
        const List<vector>& values = stencilValues[cellI];
        const scalarList& weights = filterWeights_[cellI];

        smoothI[cellI] = vector::zero;

        forAll(values, i)
        {
            smoothI[cellI] += weights[i]*values[i];
        }
    }
}


void Foam::RhieChowHighPassStab::ensureWorkingFields
(
    const volVectorField& field
) const
{
    wordList scratchPatchTypes
    (
        mesh().boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    forAll(scratchPatchTypes, patchI)
    {
        const fvPatch& patch = mesh().boundary()[patchI];

        if (patch.coupled() || fvPatch::constraintType(patch.type()))
        {
            scratchPatchTypes[patchI] = field.boundaryField()[patchI].type();
        }
    }

    if (smoothFieldPtr_.empty())
    {
        smoothFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "Dsmooth(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    writeHighPassFields_
                  ? IOobject::AUTO_WRITE
                  : IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("0", field.dimensions(), vector::zero),
                scratchPatchTypes
            )
        );
    }

    if (highFieldPtr_.empty())
    {
        highFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "Dhigh(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    writeHighPassFields_
                  ? IOobject::AUTO_WRITE
                  : IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("0", field.dimensions(), vector::zero),
                scratchPatchTypes
            )
        );
    }
}


void Foam::RhieChowHighPassStab::zeroPhysicalBoundaries
(
    volVectorField& field
) const
{
    volVectorField::Boundary& bf = field.boundaryFieldRef();

    forAll(bf, patchI)
    {
        if (!bf[patchI].coupled())
        {
            bf[patchI] == vector::zero;
        }
    }
}


void Foam::RhieChowHighPassStab::zeroPhysicalBoundaryFaces
(
    surfaceVectorField& field
) const
{
    surfaceVectorField::Boundary& bf = field.boundaryFieldRef();

    forAll(bf, patchI)
    {
        if (!bf[patchI].coupled())
        {
            bf[patchI] == vector::zero;
        }
    }
}


void Foam::RhieChowHighPassStab::reportDiagnostics
(
    const volVectorField& field,
    const volVectorField& smoothField,
    const volVectorField& highField,
    const surfaceVectorField& faceStab
) const
{
    if
    (
        !reportHighPassDiagnostics_
     || highPassReportInterval_ <= 0
     || (updateCalls_ % highPassReportInterval_) != 0
    )
    {
        return;
    }

    const vectorField& fieldI = primitiveField(field);
    const vectorField& smoothI = primitiveField(smoothField);
    const vectorField& highI = primitiveField(highField);

    scalar sumD = sum(mag(fieldI));
    scalar sumSmooth = sum(mag(smoothI));
    scalar sumHigh = sum(mag(highI));
    scalar sumD2 = sum(magSqr(fieldI));
    scalar sumHigh2 = sum(magSqr(highI));
    scalar maxHigh = highI.size() ? max(mag(highI)) : 0.0;
    label nCells = highI.size();

    reduce(sumD, sumOp<scalar>());
    reduce(sumSmooth, sumOp<scalar>());
    reduce(sumHigh, sumOp<scalar>());
    reduce(sumD2, sumOp<scalar>());
    reduce(sumHigh2, sumOp<scalar>());
    reduce(maxHigh, maxOp<scalar>());
    reduce(nCells, sumOp<label>());

    const scalar rmsD = nCells ? Foam::sqrt(sumD2/scalar(nCells)) : 0.0;
    const scalar rmsHigh = nCells ? Foam::sqrt(sumHigh2/scalar(nCells)) : 0.0;

    tmp<volVectorField> tStabDensity =
        fvc::div(mesh().magSf()*faceStab);

    vectorField stabCellForce(primitiveField(tStabDensity()));

    forAll(stabCellForce, cellI)
    {
        stabCellForce[cellI] *= mesh().V()[cellI];
    }

    const vector netForce = gSum(stabCellForce);
    const scalar sumMagForce = gSum(mag(stabCellForce));

    Info<< "rhiechowHighPass diagnostics" << nl
        << "    sumMag(D) = " << sumD << nl
        << "    sumMag(Dsmooth) = " << sumSmooth << nl
        << "    sumMag(Dhigh) = " << sumHigh << nl
        << "    sumMag(Dhigh)/sumMag(D) = "
        << sumHigh/(sumD + VSMALL) << nl
        << "    max Dhigh = " << maxHigh << nl
        << "    RMS D = " << rmsD << nl
        << "    RMS Dhigh = " << rmsHigh << nl
        << "    net extensive stabilisation force = " << netForce << nl
        << "    sum of stabilisation-force magnitudes = " << sumMagForce << nl
        << "    rank/order fallbacks = " << nOrderFallbacks_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RhieChowHighPassStab::RhieChowHighPassStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    diffStencilLaplacianStab(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
    highPassPolynomialOrder_
    (
        min
        (
            max(dict.lookupOrDefault<label>("highPassPolynomialOrder", 2), 0),
            2
        )
    ),
    highPassStencilLayers_
    (
        max(dict.lookupOrDefault<label>("highPassStencilLayers", 2), 1)
    ),
    highPassWeightExponent_
    (
        max(dict.lookupOrDefault<scalar>("highPassWeightExponent", 2.0), SMALL)
    ),
    reportHighPassDiagnostics_
    (
        dict.lookupOrDefault<Switch>("reportHighPassDiagnostics", false)
    ),
    highPassReportInterval_
    (
        max(dict.lookupOrDefault<label>("highPassReportInterval", 100), 1)
    ),
    writeHighPassFields_
    (
        dict.lookupOrDefault<Switch>("writeHighPassFields", false)
    ),
    svdMinCondition_
    (
        dict.lookupOrDefault<scalar>("highPassSvdMinCondition", 1e-12)
    ),
    activeDirections_(),
    mapPtr_(),
    stencil_(),
    filterWeights_(),
    nOrderFallbacks_(0),
    minOrderUsed_(highPassPolynomialOrder_),
    maxReproductionError_(0.0),
    smoothFieldPtr_(),
    highFieldPtr_(),
    updateCalls_(0)
{
    buildFilter();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RhieChowHighPassStab::updateVector
(
    const volVectorField& field,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();
    updateCalls_++;

    ensureWorkingFields(field);

    volVectorField& smoothField = smoothFieldPtr_();
    volVectorField& highField = highFieldPtr_();

    applyFilter(field, smoothField);

    primitiveFieldRef(highField) =
        primitiveField(field) - primitiveField(smoothField);

    highField.correctBoundaryConditions();
    zeroPhysicalBoundaries(highField);

    tmp<volTensorField> tGradHigh =
        fvc::grad(highField, "grad(" + field.name() + ")");
    volTensorField& gradHigh = tmpRef(tGradHigh);
    gradHigh.correctBoundaryConditions();

    if (faceVectorPtr().empty())
    {
        faceVectorPtr().set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "highPassFaceStabilisation(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    writeHighPassFields_
                  ? IOobject::AUTO_WRITE
                  : IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("0", dims(), vector::zero)
            )
        );
    }

    surfaceVectorField& faceStab = autoPtrRef(faceVectorPtr());

    // The difference stencil is an owner-to-neighbour face quantity and
    // changes sign if the face orientation is reversed.
    faceStab.setOriented();

    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    faceStab =
        scaleFactor_
       *(
            fvc::snGrad(highField, "snGrad(" + field.name() + ")")
          - (
                n
              & fvc::interpolate
                (
                    gradHigh,
                    "interpolate(grad(" + field.name() + "))"
                )
            )
        );

    zeroPhysicalBoundaryFaces(faceStab);

    reportDiagnostics(field, smoothField, highField, faceStab);

    if (writeHighPassFields_ && mesh().time().outputTime())
    {
        smoothField.write();
        highField.write();
        faceStab.write();

        tmp<volVectorField> tStabDensity =
            fvc::div(mesh().magSf()*faceStab);
        tmpRef(tStabDensity).rename("highPassStabForceDensity");
        tStabDensity().write();
    }
}


const Foam::volVectorField& Foam::RhieChowHighPassStab::smoothField() const
{
    if (smoothFieldPtr_.empty())
    {
        FatalErrorInFunction
            << "updateVector(...) must be called before smoothField()"
            << abort(FatalError);
    }

    return smoothFieldPtr_();
}


const Foam::volVectorField& Foam::RhieChowHighPassStab::highField() const
{
    if (highFieldPtr_.empty())
    {
        FatalErrorInFunction
            << "updateVector(...) must be called before highField()"
            << abort(FatalError);
    }

    return highFieldPtr_();
}


// ************************************************************************* //
