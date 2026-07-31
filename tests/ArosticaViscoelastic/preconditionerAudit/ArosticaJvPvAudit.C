/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#ifdef USE_PETSC
    #include "petscErrorHandling.H"
#endif
#include "fvCFD.H"
#include "IOdictionary.H"
#include "solidModel.H"
#include "nonLinGeomTotalLagTotalDispSolid.H"

using namespace Foam;

#ifdef USE_PETSC
namespace
{

scalar norm(Vec x)
{
    PetscReal value = 0.0;
    AssertPETSc(VecNorm(x, NORM_2, &value));
    return value;
}


word modeName(const label mode)
{
    switch (mode)
    {
        case 0: return "smooth displacement";
        case 1: return "pressure only";
        case 2: return "coupled displacement-pressure";
        case 3: return "rigid translation";
        case 4: return "approximate rotation";
        case 5: return "localised displacement";
        case 6: return "smooth pressure gradient";
        case 7: return "localised pressure";
        case 8: return "checkerboard pressure";
        case 9: return "one-cell pressure basis";
        case 10: return "volumetric displacement";
        case 11: return "fibre-direction deformation";
        case 12: return "sheet-direction deformation";
        case 13: return "shear deformation";
        case 14: return "constant pressure";
        case 15: return "deterministic random pressure";
        case 16: return "boundary-adjacent pressure basis";
        case 17: return "deterministic random displacement";
        case 18: return "boundary-supported displacement";
        case 19: return "smooth global displacement";
        default: return "unknown";
    }
}


void setMode
(
    Vec v,
    const label mode,
    const fvMesh& mesh
)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(v, &size));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    AssertPETSc(VecGetArray(v, &values));

    vector meshCentre(vector::zero);
    scalar meshVolume = 0.0;
    if (mode == 17 || mode == 18)
    {
        forAll(mesh.V(), cellI)
        {
            meshCentre += mesh.V()[cellI]*mesh.C()[cellI];
            meshVolume += mesh.V()[cellI];
        }
        reduce(meshCentre, sumOp<vector>());
        reduce(meshVolume, sumOp<scalar>());
        meshCentre /= meshVolume + VSMALL;
    }
    const scalar meshLength = Foam::cbrt(meshVolume + VSMALL);

    for (PetscInt i = 0; i < size; i += blockSize)
    {
        const PetscInt cellI = i/blockSize;
        values[i] = 0.0;
        values[i + 1] = 0.0;
        values[i + 2] = 0.0;
        values[i + 3] = 0.0;

        if (mode == 0) // smooth displacement
        {
            values[i] = 1.0;
            values[i + 1] = -0.7;
            values[i + 2] = 0.4;
        }
        else if (mode == 1) // pressure only
        {
            values[i + 3] = 1.0;
        }
        else if (mode == 2) // coupled displacement-pressure
        {
            values[i] = 1.0;
            values[i + 1] = 0.3;
            values[i + 2] = -0.5;
            values[i + 3] = 0.8;
        }
        else if (mode == 3) // rigid translation
        {
            values[i] = 1.0;
        }
        else if (mode == 4) // approximate rotation about the mesh origin
        {
            const vector omega(0, 0, 1);
            const vector displacement(omega ^ (mesh.C()[cellI]));
            values[i] = displacement.x();
            values[i + 1] = displacement.y();
            values[i + 2] = displacement.z();
        }
        else if (mode == 5) // localised displacement
        {
            if (cellI == 0)
            {
                values[i] = 1.0;
                values[i + 1] = -0.4;
                values[i + 2] = 0.2;
            }
        }
        else if (mode == 6) // smooth pressure gradient
        {
            values[i + 3] = 0.2 + 0.8*scalar(cellI % 101)/100.0;
        }
        else if (mode == 7) // localised pressure
        {
            if (cellI < 16)
            {
                values[i + 3] = 1.0;
            }
        }
        else if (mode == 8) // checkerboard pressure
        {
            values[i + 3] = (cellI % 2) ? -1.0 : 1.0;
        }
        else if (mode == 9) // one-cell pressure basis
        {
            if (cellI == 0)
            {
                values[i + 3] = 1.0;
            }
        }
        else if (mode == 10) // volumetric displacement
        {
            const vector displacement(mesh.C()[cellI]);
            values[i] = displacement.x();
            values[i + 1] = displacement.y();
            values[i + 2] = displacement.z();
        }
        else if (mode == 11) // fibre-direction deformation
        {
            const volVectorField& f0 =
                mesh.lookupObject<volVectorField>("f0");
            const vector fibreAxis = f0[0];
            const vector displacement =
                fibreAxis*(fibreAxis & mesh.C()[cellI]);
            values[i] = displacement.x();
            values[i + 1] = displacement.y();
            values[i + 2] = displacement.z();
        }
        else if (mode == 12) // sheet-direction deformation
        {
            const volVectorField& s0 =
                mesh.lookupObject<volVectorField>("s0");
            const vector sheetAxis = s0[0];
            const vector displacement =
                sheetAxis*(sheetAxis & mesh.C()[cellI]);
            values[i] = displacement.x();
            values[i + 1] = displacement.y();
            values[i + 2] = displacement.z();
        }
        else if (mode == 13) // fibre-sheet shear
        {
            const volVectorField& f0 =
                mesh.lookupObject<volVectorField>("f0");
            const volVectorField& s0 =
                mesh.lookupObject<volVectorField>("s0");
            const vector fibreAxis = f0[0];
            const vector sheetAxis = s0[0];
            const vector displacement =
                fibreAxis*(sheetAxis & mesh.C()[cellI]);
            values[i] = displacement.x();
            values[i + 1] = displacement.y();
            values[i + 2] = displacement.z();
        }
        else if (mode == 14) // constant pressure
        {
            values[i + 3] = 1.0;
        }
        else if (mode == 15) // deterministic random pressure
        {
            values[i + 3] =
                scalar((37*cellI + 17*(cellI % 11) + 5) % 101 - 50)
               /50.0;
        }
        else if (mode == 16) // boundary-adjacent pressure basis
        {
            label boundaryCell = -1;
            forAll(mesh.boundary(), patchI)
            {
                const labelUList& patchCells =
                    mesh.boundary()[patchI].faceCells();
                if (!patchCells.empty())
                {
                    boundaryCell = patchCells[0];
                    break;
                }
            }
            if (cellI == boundaryCell)
            {
                values[i + 3] = 1.0;
            }
        }
        else if (mode == 17) // deterministic random displacement
        {
            // A fixed-seed, smooth pseudo-random combination of spatial
            // modes.  Unlike cellwise white noise, this is an admissible
            // displacement perturbation on successively refined meshes.
            const vector q = (mesh.C()[cellI] - meshCentre)/meshLength;
            values[i] = Foam::sin(1.37*q.x() - 0.83*q.y() + 0.41*q.z());
            values[i + 1] =
                Foam::sin(-0.59*q.x() + 1.61*q.y() + 0.73*q.z());
            values[i + 2] =
                Foam::sin(0.47*q.x() + 0.31*q.y() - 1.43*q.z());
        }
        else if (mode == 18) // boundary-supported displacement
        {
            // Smoothly weight the field towards the exterior while keeping
            // a continuous cell-to-cell direction for the boundary tangent.
            const vector q = (mesh.C()[cellI] - meshCentre)/meshLength;
            const scalar weight = sqr(magSqr(q))/(1.0 + sqr(magSqr(q)));
            values[i] = weight*(0.4 + 0.2*q.x());
            values[i + 1] = weight*(-0.3 + 0.1*q.y());
            values[i + 2] = weight*(0.2 - 0.15*q.z());
        }
        else if (mode == 19) // smooth global displacement
        {
            const vector& c = mesh.C()[cellI];
            values[i] = 0.5*c.x() + 0.1*c.y() - 0.2*c.z();
            values[i + 1] = -0.15*c.x() + 0.4*c.y() + 0.05*c.z();
            values[i + 2] = 0.1*c.x() - 0.2*c.y() + 0.3*c.z();
        }
    }

    AssertPETSc(VecRestoreArray(v, &values));
}


scalar checksum(const volSymmTensorField& field)
{
    scalar result = 0.0;
    forAll(field, cellI)
    {
        result += (cellI + 1)*magSqr(field[cellI]);
    }

    forAll(field.boundaryField(), patchI)
    {
        const symmTensorField& patchField = field.boundaryField()[patchI];
        forAll(patchField, faceI)
        {
            result +=
                (field.size() + faceI + patchI + 1)
                *magSqr(patchField[faceI]);
        }
    }

    return result;
}


scalar checksum(const surfaceSymmTensorField& field)
{
    scalar result = 0.0;
    forAll(field, faceI)
    {
        result += (faceI + 1)*magSqr(field[faceI]);
    }

    forAll(field.boundaryField(), patchI)
    {
        const symmTensorField& patchField = field.boundaryField()[patchI];
        forAll(patchField, faceI)
        {
            result +=
                (field.size() + faceI + patchI + 1)
                *magSqr(patchField[faceI]);
        }
    }

    return result;
}


void historyNames
(
    const fvMesh& mesh,
    wordList& names
)
{
    if (!mesh.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return;
    }

    const IOdictionary& mechanicalProperties =
        mesh.lookupObject<IOdictionary>("mechanicalProperties");

    if (!mechanicalProperties.found("mechanical"))
    {
        return;
    }

    const PtrList<entry> lawEntries
    (
        mechanicalProperties.lookup("mechanical")
    );
    const wordList prefixes
    (
        makeList<word>({"EOld", "EOldOld", "EfOld", "EfOldOld"})
    );

    forAll(lawEntries, lawI)
    {
        forAll(prefixes, prefixI)
        {
            const word fieldName
            (
                "Arostica" + prefixes[prefixI] + '_'
              + lawEntries[lawI].keyword()
            );

            if
            (
                mesh.foundObject<volSymmTensorField>(fieldName)
             || mesh.foundObject<surfaceSymmTensorField>(fieldName)
            )
            {
                names.append(fieldName);
            }
        }
    }
}


scalar historyChecksum
(
    const fvMesh& mesh,
    const word& name
)
{
    if (mesh.foundObject<volSymmTensorField>(name))
    {
        return checksum(mesh.lookupObject<volSymmTensorField>(name));
    }

    if (mesh.foundObject<surfaceSymmTensorField>(name))
    {
        return checksum(mesh.lookupObject<surfaceSymmTensorField>(name));
    }

    return 0.0;
}


scalar maximumHistoryDifference
(
    const fvMesh& mesh,
    const wordList& names,
    const scalarField& reference
)
{
    scalar maximum = 0.0;
    forAll(names, nameI)
    {
        maximum = max
        (
            maximum,
            mag(historyChecksum(mesh, names[nameI]) - reference[nameI])
        );
    }
    return maximum;
}


scalar dotCosine(Vec a, Vec b)
{
    PetscScalar dot = 0.0;
    const scalar na = norm(a);
    const scalar nb = norm(b);
    AssertPETSc(VecDot(a, b, &dot));
    return PetscRealPart(dot)/(na*nb + VSMALL);
}


scalar firstComponentsNorm(const scalarField& components)
{
    scalar result = 0.0;
    const label nComponents = components.size() < 3 ? components.size() : 3;
    for (label cmptI = 0; cmptI < nComponents; ++cmptI)
    {
        result += sqr(components[cmptI]);
    }
    return Foam::sqrt(result);
}


void componentNorms(Vec v, scalarField& result)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(v, &size));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    result.setSize(blockSize, 0.0);
    AssertPETSc(VecGetArrayRead(v, &values));
    for (PetscInt i = 0; i < size; i += blockSize)
    {
        for (PetscInt cmptI = 0; cmptI < blockSize; ++cmptI)
        {
            result[cmptI] += sqr(PetscRealPart(values[i + cmptI]));
        }
    }
    AssertPETSc(VecRestoreArrayRead(v, &values));
    forAll(result, cmptI)
    {
        result[cmptI] = Foam::sqrt(result[cmptI]);
    }
}


void pressureComponents(Vec v, scalarField& result)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(v, &size));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    result.setSize(size/blockSize, 0.0);
    AssertPETSc(VecGetArrayRead(v, &values));
    forAll(result, cellI)
    {
        result[cellI] =
            PetscRealPart(values[cellI*blockSize + blockSize - 1]);
    }
    AssertPETSc(VecRestoreArrayRead(v, &values));
}


scalar fieldNorm(const scalarField& values)
{
    scalar sumSqr = 0.0;
    forAll(values, valueI)
    {
        sumSqr += sqr(values[valueI]);
    }
    reduce(sumSqr, sumOp<scalar>());
    return Foam::sqrt(sumSqr);
}


void reportComponentStage
(
    const word& stage,
    const scalarField& trueLocal,
    const scalarField& trueStabilisation,
    const scalarField& trueTotal,
    const scalarField& approximateLocal,
    const scalarField& approximateStabilisation,
    const scalarField& approximateTotal
)
{
    Info<< "        " << stage
        << ": true(local,stab,total)=("
        << fieldNorm(trueLocal) << ','
        << fieldNorm(trueStabilisation) << ','
        << fieldNorm(trueTotal) << "), P(local,stab,total)=("
        << fieldNorm(approximateLocal) << ','
        << fieldNorm(approximateStabilisation) << ','
        << fieldNorm(approximateTotal) << ')' << endl;
}


void setDisplacementPart(Vec destination, Vec source)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    const PetscScalar* sourceValues = nullptr;
    AssertPETSc(VecGetLocalSize(destination, &size));
    AssertPETSc(VecGetBlockSize(destination, &blockSize));
    AssertPETSc(VecGetArray(destination, &values));
    AssertPETSc(VecGetArrayRead(source, &sourceValues));
    for (PetscInt i = 0; i < size; i += blockSize)
    {
        for (PetscInt cmptI = 0; cmptI < blockSize; ++cmptI)
        {
            values[i + cmptI] =
                cmptI < 3 ? sourceValues[i + cmptI] : 0.0;
        }
    }
    AssertPETSc(VecRestoreArrayRead(source, &sourceValues));
    AssertPETSc(VecRestoreArray(destination, &values));
}


void setPrestressedState
(
    Vec x,
    const fvMesh& mesh,
    const scalar scale
)
{
    const scalar coordinateScale = scale;

    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(x, &size));
    AssertPETSc(VecGetBlockSize(x, &blockSize));
    AssertPETSc(VecGetArray(x, &values));

    for (PetscInt cellI = 0; cellI < mesh.nCells(); ++cellI)
    {
        const vector& c = mesh.C()[cellI];
        const vector displacement
        (
            coordinateScale*(0.5*c.x() + 0.1*c.y()),
            coordinateScale*(-0.05*c.x() + 0.4*c.y()),
            coordinateScale*(0.1*c.x() + 0.15*c.z())
        );
        const PetscInt offset = blockSize*cellI;
        values[offset] = displacement.x();
        values[offset + 1] = displacement.y();
        values[offset + 2] = displacement.z();
        values[offset + 3] = 0.0;
    }

    AssertPETSc(VecRestoreArray(x, &values));
}


bool auditPrestressedState(const IOdictionary& solidProperties)
{
    const word coeffsName
    (
        "arosticaNonLinearGeometryTotalLagrangianTotalDisplacementCoeffs"
    );

    return
        solidProperties.found(coeffsName)
     && solidProperties.subDict(coeffsName).lookupOrDefault<Switch>
        (
            "auditPrestressedState", false
        );
}


void setCentralTerm
(
    Vec destination,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& plus,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& minus,
    const label termI,
    const scalar epsilon,
    const fvMesh& mesh
)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(destination, &size));
    AssertPETSc(VecGetBlockSize(destination, &blockSize));
    AssertPETSc(VecGetArray(destination, &values));
    for (PetscInt cellI = 0; cellI < mesh.nCells(); ++cellI)
    {
        const vector action
        (
            (
                plus.cellActions[termI][cellI]
              - minus.cellActions[termI][cellI]
            )/(2.0*epsilon)*mesh.V()[cellI]
        );
        const PetscInt offset = blockSize*cellI;
        values[offset] = action.x();
        values[offset + 1] = action.y();
        values[offset + 2] = action.z();
        values[offset + 3] = 0.0;
    }
    AssertPETSc(VecRestoreArray(destination, &values));
}


scalar centralFaceActionNorm
(
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& plus,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& minus,
    const label termI,
    const scalar epsilon
)
{
    scalar sumSqr = 0.0;
    forAll(plus.faceForces[termI], faceI)
    {
        const vector action
        (
            (
                plus.faceForces[termI][faceI]
              - minus.faceForces[termI][faceI]
            )/(2.0*epsilon)
        );
        sumSqr += magSqr(action);
    }
    reduce(sumSqr, sumOp<scalar>());
    return Foam::sqrt(sumSqr);
}


void centralFaceCategoryNorms
(
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& plus,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& minus,
    const label termI,
    const scalar epsilon,
    const fvMesh& mesh,
    scalar& internalFaceNorm,
    scalar& boundaryFaceNorm
)
{
    scalar internalSqr = 0.0;
    scalar boundarySqr = 0.0;
    forAll(plus.faceForces[termI], faceI)
    {
        const vector action
        (
            (
                plus.faceForces[termI][faceI]
              - minus.faceForces[termI][faceI]
            )/(2.0*epsilon)
        );
        if (faceI < mesh.nInternalFaces())
        {
            internalSqr += magSqr(action);
        }
        else
        {
            boundarySqr += magSqr(action);
        }
    }
    reduce(internalSqr, sumOp<scalar>());
    reduce(boundarySqr, sumOp<scalar>());
    internalFaceNorm = Foam::sqrt(internalSqr);
    boundaryFaceNorm = Foam::sqrt(boundarySqr);
}


void reportLargestCellsAndFaces
(
    Vec jvD,
    Vec pv,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& plus,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& minus,
    const scalar epsilon,
    const fvMesh& mesh
)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    const PetscScalar* jvValues = nullptr;
    const PetscScalar* pvValues = nullptr;
    AssertPETSc(VecGetLocalSize(jvD, &size));
    AssertPETSc(VecGetBlockSize(jvD, &blockSize));
    AssertPETSc(VecGetArrayRead(jvD, &jvValues));
    AssertPETSc(VecGetArrayRead(pv, &pvValues));

    scalar maxCellError = 0.0;
    label maxCell = -1;
    for (label cellI = 0; cellI < mesh.nCells(); ++cellI)
    {
        const label offset = blockSize*cellI;
        const vector error
        (
            jvValues[offset] - pvValues[offset],
            jvValues[offset + 1] - pvValues[offset + 1],
            jvValues[offset + 2] - pvValues[offset + 2]
        );
        if (mag(error) > maxCellError)
        {
            maxCellError = mag(error);
            maxCell = cellI;
        }
    }
    AssertPETSc(VecRestoreArrayRead(pv, &pvValues));
    AssertPETSc(VecRestoreArrayRead(jvD, &jvValues));

    reduce(maxCellError, maxOp<scalar>());

    scalar maxFaceAction = 0.0;
    label maxFace = -1;
    forAll(plus.faceForces[0], faceI)
    {
        vector action = vector::zero;
        for (label termI = 0; termI < plus.names.size(); ++termI)
        {
            action +=
                (
                    plus.faceForces[termI][faceI]
                  - minus.faceForces[termI][faceI]
                )/(2.0*epsilon);
        }
        if (mag(action) > maxFaceAction)
        {
            maxFaceAction = mag(action);
            maxFace = faceI;
        }
    }

    Info<< "    largest localised Jv-Pv cell action = " << maxCellError
        << " at local cell " << maxCell << nl
        << "    largest decomposed face action = " << maxFaceAction
        << " at global face " << maxFace;
    if (maxFace >= 0 && maxFace < mesh.nInternalFaces())
    {
        Info<< " (internal, owner=" << mesh.owner()[maxFace]
            << ", neighbour=" << mesh.neighbour()[maxFace] << ')';
    }
    Info<< endl;
}


scalar bestFitScalar(Vec target, Vec trial)
{
    PetscScalar dot = 0.0;
    AssertPETSc(VecDot(target, trial, &dot));
    const scalar denominator = sqr(norm(trial));
    return denominator > VSMALL ? PetscRealPart(dot)/denominator : 0.0;
}


scalar bestFitError(Vec target, Vec trial, const scalar scale)
{
    Vec difference = nullptr;
    AssertPETSc(VecDuplicate(target, &difference));
    AssertPETSc(VecCopy(target, difference));
    AssertPETSc(VecAXPY(difference, -scale, trial));
    const scalar error = norm(difference)/(norm(target) + VSMALL);
    AssertPETSc(VecDestroy(&difference));
    return error;
}


void reportGeometrySplit
(
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& base,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& plus,
    const solidModels::nonLinGeomTotalLagTotalDispSolid::
        momentumResidualDecompositionData& minus,
    const scalar epsilon,
    const fvMesh& mesh
)
{
    scalar materialSqr = 0.0;
    scalar geometrySqr = 0.0;
    scalar crossSqr = 0.0;
    const label nInternalFaces = mesh.nInternalFaces();

    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        const symmTensor dStress
        (
            (
                plus.directPassiveFaceStress[faceI]
              - minus.directPassiveFaceStress[faceI]
            )/(2.0*epsilon)
        );
        const vector dSf
        (
            (
                plus.faceAreaVectors[faceI]
              - minus.faceAreaVectors[faceI]
            )/(2.0*epsilon)
        );
        const vector material = base.faceAreaVectors[faceI] & dStress;
        const vector geometry = dSf & base.directPassiveFaceStress[faceI];
        const vector exact
        (
            (
                plus.faceForces[0][faceI]
              - minus.faceForces[0][faceI]
            )/(2.0*epsilon)
        );
        materialSqr += magSqr(material);
        geometrySqr += magSqr(geometry);
        crossSqr += magSqr(exact - material - geometry);
    }

    reduce(materialSqr, sumOp<scalar>());
    reduce(geometrySqr, sumOp<scalar>());
    reduce(crossSqr, sumOp<scalar>());

    Info<< "    passive direct-face material/geometric norms = "
        << Foam::sqrt(materialSqr) << " / " << Foam::sqrt(geometrySqr)
        << ", cross remainder = " << Foam::sqrt(crossSqr) << endl;
}


void reportPressureToDisplacementStats(const Mat P)
{
    PetscInt rows = 0;
    PetscInt cols = 0;
    AssertPETSc(MatGetSize(P, &rows, &cols));

    PetscInt nonZero = 0;
    scalar sumSqr = 0.0;
    scalar maxAbs = 0.0;
    scalar minNonZero = GREAT;
    scalar minRowSum = GREAT;
    scalar maxRowSum = -GREAT;
    scalar minColSum = GREAT;
    scalar maxColSum = -GREAT;
    List<scalar> columnSums(cols, 0.0);

    for (PetscInt row = 0; row < rows; ++row)
    {
        if (row % 4 == 3)
        {
            continue;
        }

        PetscInt nCols = 0;
        const PetscInt* colIndices = nullptr;
        const PetscScalar* values = nullptr;
        AssertPETSc(MatGetRow(P, row, &nCols, &colIndices, &values));
        scalar rowSum = 0.0;
        for (PetscInt entryI = 0; entryI < nCols; ++entryI)
        {
            const PetscInt col = colIndices[entryI];
            if (col % 4 == 3)
            {
                const scalar value = PetscRealPart(values[entryI]);
                rowSum += value;
                columnSums[col] += value;
                sumSqr += sqr(value);
                maxAbs = max(maxAbs, mag(value));
                if (mag(value) > SMALL)
                {
                    ++nonZero;
                    minNonZero = min(minNonZero, mag(value));
                }
            }
        }
        minRowSum = min(minRowSum, rowSum);
        maxRowSum = max(maxRowSum, rowSum);
        AssertPETSc(MatRestoreRow(P, row, &nCols, &colIndices, &values));
    }

    forAll(columnSums, colI)
    {
        minColSum = min(minColSum, columnSums[colI]);
        maxColSum = max(maxColSum, columnSums[colI]);
    }

    Info<< "P_Dp matrix statistics:" << nl
        << "    nonzero count = " << nonZero << nl
        << "    Frobenius norm = " << Foam::sqrt(sumSqr) << nl
        << "    max abs coefficient = " << maxAbs << nl
        << "    min nonzero coefficient = " << minNonZero << nl
        << "    row-sum range = [" << minRowSum << ", " << maxRowSum
        << "]" << nl
        << "    column-sum range = [" << minColSum << ", " << maxColSum
        << "]" << endl;
}


void setBlockPart
(
    Vec destination,
    Vec source,
    const bool pressurePart
)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* destinationValues = nullptr;
    const PetscScalar* sourceValues = nullptr;
    AssertPETSc(VecGetLocalSize(destination, &size));
    AssertPETSc(VecGetBlockSize(destination, &blockSize));
    AssertPETSc(VecGetArray(destination, &destinationValues));
    AssertPETSc(VecGetArrayRead(source, &sourceValues));

    for (PetscInt i = 0; i < size; i += blockSize)
    {
        for (PetscInt cmptI = 0; cmptI < blockSize; ++cmptI)
        {
            const bool isPressure = cmptI == blockSize - 1;
            destinationValues[i + cmptI] =
                isPressure == pressurePart ? sourceValues[i + cmptI] : 0.0;
        }
    }

    AssertPETSc(VecRestoreArrayRead(source, &sourceValues));
    AssertPETSc(VecRestoreArray(destination, &destinationValues));
}


word blockName(const label rowBlock, const label colBlock)
{
    if (rowBlock == 3 && colBlock == 3)
    {
        return "P_pp";
    }

    if (rowBlock == 3)
    {
        return "P_pD";
    }

    if (colBlock == 3)
    {
        return "P_Dp";
    }

    return "P_DD";
}


void reportMixedBlockStats(const Mat P, const label blockSize)
{
    PetscInt rows = 0;
    PetscInt cols = 0;
    PetscInt rowStart = 0;
    PetscInt rowEnd = 0;
    AssertPETSc(MatGetSize(P, &rows, &cols));
    AssertPETSc(MatGetOwnershipRange(P, &rowStart, &rowEnd));

    for (label rowGroup = 0; rowGroup < 2; ++rowGroup)
    {
        const label rowBlock = rowGroup == 0 ? 0 : blockSize - 1;
        const label rowBlockEnd = rowGroup == 0 ? blockSize - 1 : blockSize;

        for (label colGroup = 0; colGroup < 2; ++colGroup)
        {
            const label colBlock = colGroup == 0 ? 0 : blockSize - 1;
            const label colBlockEnd =
                colGroup == 0 ? blockSize - 1 : blockSize;
            PetscInt nonZero = 0;
            scalar sumSqr = 0.0;
            scalar maxAbs = 0.0;
            scalar minNonZero = GREAT;
            scalar minDiag = GREAT;
            scalar maxDiag = -GREAT;
            scalar minRowSum = GREAT;
            scalar maxRowSum = -GREAT;
            List<scalar> columnSums(cols, 0.0);

            for (PetscInt row = rowStart; row < rowEnd; ++row)
            {
                if
                (
                    row % blockSize < rowBlock
                 || row % blockSize >= rowBlockEnd
                )
                {
                    continue;
                }

                PetscInt nCols = 0;
                const PetscInt* colIndices = nullptr;
                const PetscScalar* values = nullptr;
                AssertPETSc
                (
                    MatGetRow(P, row, &nCols, &colIndices, &values)
                );

                scalar rowSum = 0.0;
                for (PetscInt entryI = 0; entryI < nCols; ++entryI)
                {
                    const PetscInt col = colIndices[entryI];
                    if
                    (
                        col % blockSize < colBlock
                     || col % blockSize >= colBlockEnd
                    )
                    {
                        continue;
                    }

                    const scalar value = PetscRealPart(values[entryI]);
                    rowSum += value;
                    columnSums[col] += value;
                    sumSqr += sqr(value);
                    maxAbs = max(maxAbs, mag(value));
                    if (mag(value) > SMALL)
                    {
                        ++nonZero;
                        minNonZero = min(minNonZero, mag(value));
                    }
                    if (row == col)
                    {
                        minDiag = min(minDiag, value);
                        maxDiag = max(maxDiag, value);
                    }
                }

                minRowSum = min(minRowSum, rowSum);
                maxRowSum = max(maxRowSum, rowSum);
                AssertPETSc
                (
                    MatRestoreRow
                    (
                        P, row, &nCols, &colIndices, &values
                    )
                );
            }

            scalar minColSum = GREAT;
            scalar maxColSum = -GREAT;
            forAll(columnSums, colI)
            {
                if
                (
                    colI % blockSize >= colBlock
                 && colI % blockSize < colBlockEnd
                )
                {
                    minColSum = min(minColSum, columnSums[colI]);
                    maxColSum = max(maxColSum, columnSums[colI]);
                }
            }

            Info<< blockName(rowBlock, colBlock)
                << " matrix statistics:" << nl
                << "    nonzero count = " << nonZero << nl
                << "    Frobenius norm = " << Foam::sqrt(sumSqr) << nl
                << "    max abs coefficient = " << maxAbs << nl
                << "    min nonzero coefficient = " << minNonZero << nl
                << "    diagonal range = [" << minDiag << ", " << maxDiag
                << "]" << nl
                << "    row-sum range = [" << minRowSum << ", "
                << maxRowSum << "]" << nl
                << "    column-sum range = [" << minColSum << ", "
                << maxColSum << "]" << endl;
        }
    }
}

}
#endif


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "pressureOnly",
        "Run only the mixed pressure-block audit"
    );
    #include "setRootCase.H"
    #include "createTime.H"

#ifndef USE_PETSC
    FatalErrorInFunction
        << "ArosticaJvPvAudit requires a PETSc-enabled build"
        << abort(FatalError);
#else
    const bool pressureOnly = args.optionFound("pressureOnly");

    autoPtr<solidModel> solid =
        solidModel::New(runTime, dynamicFvMesh::defaultRegion);

    solidModels::nonLinGeomTotalLagTotalDispSolid& petscSolid =
        refCast<solidModels::nonLinGeomTotalLagTotalDispSolid>(*solid);

    Vec x = nullptr;
    Vec xPlus = nullptr;
    Vec xMinus = nullptr;
    Vec residual = nullptr;
    Vec residualRepeat = nullptr;
    Vec residualPlus = nullptr;
    Vec residualMinus = nullptr;
    Vec jv = nullptr;
    Vec pv = nullptr;
    Vec difference = nullptr;
    Vec modeVector = nullptr;
    Mat P = nullptr;

    AssertPETSc(petscSolid.initialiseSolution(x));
    const IOdictionary solidProperties
    (
        IOobject
        (
            "solidProperties",
            runTime.constant(),
            petscSolid.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word solidCoeffsName
    (
        "arosticaNonLinearGeometryTotalLagrangianTotalDisplacementCoeffs"
    );
    const dictionary& solidCoeffs = solidProperties.subDict(solidCoeffsName);
    const scalar Vtot = gSum(petscSolid.mesh().V());
    const scalar L0 = Foam::cbrt(Vtot);
    const tmp<volScalarField> tTwoMu
    (
        2.0*petscSolid.mechanical().shearModulus()
    );
    scalar twoMuV = 0.0;
    forAll(tTwoMu(), cellI)
    {
        twoMuV += tTwoMu()[cellI]*petscSolid.mesh().V()[cellI];
    }
    reduce(twoMuV, sumOp<scalar>());
    const scalar twoMuRef = twoMuV/Vtot;
    const scalar pressureEqnScale =
        solidCoeffs.lookupOrDefault<scalar>("pressureScaleFactor", 1.0)
       *(
            solidCoeffs.lookupOrDefault<Switch>
            (
                "pressureScaleByTwoMu", true
            )
          ? twoMuRef
          : 1.0
        );
    scalar pressureUnknownScale = 1.0;
    if
    (
        solidCoeffs.lookupOrDefault<Switch>
        (
            "scaleMixedPetScFields", true
        )
    )
    {
        const word pressureUnknownScaleType =
            solidCoeffs.lookupOrDefault<word>
            (
                "pressureUnknownScale", "twoMu"
            );
        if
        (
            pressureUnknownScaleType == "twoMu"
         || pressureUnknownScaleType == "2mu"
        )
        {
            pressureUnknownScale = twoMuRef;
        }
        else if
        (
            pressureUnknownScaleType == "user"
         || pressureUnknownScaleType == "scalar"
        )
        {
            pressureUnknownScale = readScalar
            (
                solidCoeffs.lookup("pressureUnknownScaleValue")
            );
        }
    }
    const word pressureRowScaling = solidCoeffs.lookupOrDefault<word>
    (
        "pressureRowScaling", "legacy"
    );
    const tmp<volScalarField> tBulkModulus
    (
        petscSolid.mechanical().bulkModulus()
    );

    if (auditPrestressedState(solidProperties))
    {
        const scalar prestressedScale = solidCoeffs.lookupOrDefault<scalar>
        (
            "auditPrestressedScale", 1.0e-8
        );
        setPrestressedState(x, petscSolid.mesh(), prestressedScale);
        Info<< "Using deterministic affine prestressed audit state" << nl
            << "    coordinate scale = " << prestressedScale << nl
            << "    D = scale*(0.5*x + 0.1*y, -0.05*x + 0.4*y, "
            << "0.1*x + 0.15*z)" << endl;
    }
    AssertPETSc(VecDuplicate(x, &xPlus));
    AssertPETSc(VecDuplicate(x, &xMinus));
    AssertPETSc(VecDuplicate(x, &residual));
    AssertPETSc(VecDuplicate(x, &residualRepeat));
    AssertPETSc(VecDuplicate(x, &residualPlus));
    AssertPETSc(VecDuplicate(x, &residualMinus));
    AssertPETSc(VecDuplicate(x, &jv));
    AssertPETSc(VecDuplicate(x, &pv));
    AssertPETSc(VecDuplicate(x, &difference));
    AssertPETSc(VecDuplicate(x, &modeVector));
    AssertPETSc(petscSolid.initialiseJacobian(P));

    AssertPETSc(petscSolid.formResidual(residual, x));
    AssertPETSc(petscSolid.formResidual(residualRepeat, x));
    AssertPETSc(VecCopy(residualRepeat, difference));
    AssertPETSc(VecAXPY(difference, -1.0, residual));
    Info<< "Baseline repeated-residual difference = " << norm(difference)
        << endl;

    AssertPETSc(petscSolid.formJacobian(P, x));
    AssertPETSc(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
    AssertPETSc(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
    reportMixedBlockStats(P, 4);

    wordList historyFieldNames;
    historyNames(petscSolid.mesh(), historyFieldNames);
    scalarField historyReference(historyFieldNames.size(), 0.0);
    forAll(historyFieldNames, fieldI)
    {
        historyReference[fieldI] =
            historyChecksum(petscSolid.mesh(), historyFieldNames[fieldI]);
    }

    Info<< "Accepted history fields tracked = " << historyFieldNames << nl;
    forAll(historyFieldNames, fieldI)
    {
        Info<< "    " << historyFieldNames[fieldI]
            << " checksum = " << historyReference[fieldI] << nl;
    }

    scalarField epsilons(5);
    epsilons[0] = 1e-4;
    epsilons[1] = 1e-5;
    epsilons[2] = 1e-6;
    epsilons[3] = 1e-7;
    epsilons[4] = 1e-8;

    const label nModes = pressureOnly ? 0 : 20;
    scalar maximumHistoryChange = 0.0;
    scalar maximumResidualRepeatDifference = 0.0;

    Info<< "Arostica Jv/Pv action audit" << nl
        << "    central-difference epsilons = " << epsilons << nl
        << "    rows = " << petscSolid.mesh().nCells() << " cells" << endl;

    for (label mode = 0; mode < nModes; ++mode)
    {
        AssertPETSc(VecSet(modeVector, 0.0));
        AssertPETSc(VecSet(jv, 0.0));
        AssertPETSc(VecSet(pv, 0.0));
        setMode(modeVector, mode, petscSolid.mesh());
        const scalar modeNorm = norm(modeVector);
        if (modeNorm > VSMALL)
        {
            AssertPETSc(VecScale(modeVector, 1.0/modeNorm));
        }
        AssertPETSc(MatMult(P, modeVector, pv));

        forAll(epsilons, epsilonI)
        {
            const scalar epsilon = epsilons[epsilonI];

            AssertPETSc(VecCopy(x, xPlus));
            AssertPETSc(VecAXPY(xPlus, epsilon, modeVector));
            AssertPETSc(VecCopy(x, xMinus));
            AssertPETSc(VecAXPY(xMinus, -epsilon, modeVector));

            AssertPETSc(petscSolid.formResidual(residualPlus, xPlus));
            scalar historyChange = maximumHistoryDifference
            (
                petscSolid.mesh(), historyFieldNames, historyReference
            );
            maximumHistoryChange = max(maximumHistoryChange, historyChange);

            AssertPETSc(petscSolid.formResidual(residualMinus, xMinus));
            historyChange = maximumHistoryDifference
            (
                petscSolid.mesh(), historyFieldNames, historyReference
            );
            maximumHistoryChange = max(maximumHistoryChange, historyChange);

            AssertPETSc(VecCopy(residualPlus, jv));
            AssertPETSc(VecAXPY(jv, -1.0, residualMinus));
            AssertPETSc(VecScale(jv, 1.0/(2.0*epsilon)));

            AssertPETSc(VecCopy(jv, difference));
            AssertPETSc(VecAXPY(difference, -1.0, pv));

            const scalar jvNorm = norm(jv);
            const scalar diffNorm = norm(difference);
            scalarField jvComponents;
            scalarField pvComponents;
            componentNorms(jv, jvComponents);
            componentNorms(pv, pvComponents);
            const scalar jvDNorm = firstComponentsNorm(jvComponents);
            const scalar pvDNorm = firstComponentsNorm(pvComponents);

            Info<< "    mode " << mode << " (" << modeName(mode) << ')'
                << ", epsilon = " << epsilon
                << ": ||J_DD^FD v|| = " << jvDNorm
                << ", ||J_pD^FD v|| = " << jvComponents[3]
                << ", ||P_DD v|| = " << pvDNorm
                << ", ||P_pD v|| = " << pvComponents[3]
                << ", full relative difference = "
                << diffNorm/(jvNorm + VSMALL)
                << ", full cosine = " << dotCosine(jv, pv)
                << ", J components = " << jvComponents
                << ", P components = " << pvComponents << nl;
        }

        AssertPETSc(petscSolid.formResidual(residualRepeat, x));
        AssertPETSc(VecCopy(residualRepeat, difference));
        AssertPETSc(VecAXPY(difference, -1.0, residual));
        maximumResidualRepeatDifference = max
        (
            maximumResidualRepeatDifference,
            norm(difference)
        );
    }

    Vec jvBlock = nullptr;
    Vec pvBlock = nullptr;
    Vec blockDifference = nullptr;
    AssertPETSc(VecDuplicate(x, &jvBlock));
    AssertPETSc(VecDuplicate(x, &pvBlock));
    AssertPETSc(VecDuplicate(x, &blockDifference));

    const label displacementBlockModes[] =
    {
        0, 3, 4, 5, 10, 11, 12, 13
    };
    const label pressureBlockModes[] =
    {
        14, 6, 7, 8, 9, 15, 16
    };
    const label nDisplacementBlockModes =
        pressureOnly ? 0 : sizeof(displacementBlockModes)/sizeof(label);
    scalarField blockAuditEpsilons(2);
    blockAuditEpsilons[0] = 1e-6;
    blockAuditEpsilons[1] = 1e-7;

    Info<< "\nMixed pressure-block FD action audit" << nl
        << "    P_pD uses displacement-only perturbations" << nl
        << "    P_pp uses pressure-only perturbations" << endl;

    scalar maximumRequiredPppRelativeError = 0.0;
    scalar minimumRequiredPppCosine = GREAT;

    for (label modeI = 0;
         modeI < nDisplacementBlockModes;
         ++modeI)
    {
        const label mode = displacementBlockModes[modeI];
        AssertPETSc(VecSet(modeVector, 0.0));
        setMode(modeVector, mode, petscSolid.mesh());
        const scalar modeNorm = norm(modeVector);
        if (modeNorm > VSMALL)
        {
            AssertPETSc(VecScale(modeVector, 1.0/modeNorm));
        }
        AssertPETSc(MatMult(P, modeVector, pv));
        setBlockPart(pvBlock, pv, true);

        forAll(blockAuditEpsilons, epsilonI)
        {
            const scalar epsilon = blockAuditEpsilons[epsilonI];
            AssertPETSc(VecCopy(x, xPlus));
            AssertPETSc(VecAXPY(xPlus, epsilon, modeVector));
            AssertPETSc(VecCopy(x, xMinus));
            AssertPETSc(VecAXPY(xMinus, -epsilon, modeVector));
            AssertPETSc(petscSolid.formResidual(residualPlus, xPlus));
            AssertPETSc(petscSolid.formResidual(residualMinus, xMinus));
            AssertPETSc(VecCopy(residualPlus, jv));
            AssertPETSc(VecAXPY(jv, -1.0, residualMinus));
            AssertPETSc(VecScale(jv, 1.0/(2.0*epsilon)));
            setBlockPart(jvBlock, jv, true);
            AssertPETSc(VecCopy(jvBlock, blockDifference));
            AssertPETSc(VecAXPY(blockDifference, -1.0, pvBlock));

            Info<< "    P_pD mode " << mode << " (" << modeName(mode)
                << "), epsilon = " << epsilon
                << ": ||J_pD v|| = " << norm(jvBlock)
                << ", ||P_pD v|| = " << norm(pvBlock)
                << ", relative error = "
                << norm(blockDifference)/(norm(jvBlock) + VSMALL)
                << ", cosine = " << dotCosine(jvBlock, pvBlock)
                << ", best-fit scalar = "
                << bestFitScalar(jvBlock, pvBlock)
                << ", scaled error = "
                << bestFitError
                   (
                       jvBlock,
                       pvBlock,
                       bestFitScalar(jvBlock, pvBlock)
                   ) << endl;

            maximumHistoryChange = max
            (
                maximumHistoryChange,
                maximumHistoryDifference
                (
                    petscSolid.mesh(), historyFieldNames, historyReference
                )
            );
        }
    }

    for (label modeI = 0;
         modeI < sizeof(pressureBlockModes)/sizeof(label);
         ++modeI)
    {
        const label mode = pressureBlockModes[modeI];
        AssertPETSc(VecSet(modeVector, 0.0));
        setMode(modeVector, mode, petscSolid.mesh());
        const scalar modeNorm = norm(modeVector);
        if (modeNorm > VSMALL)
        {
            AssertPETSc(VecScale(modeVector, 1.0/modeNorm));
        }
        AssertPETSc(MatMult(P, modeVector, pv));
        setBlockPart(pvBlock, pv, true);

        forAll(blockAuditEpsilons, epsilonI)
        {
            const scalar epsilon = blockAuditEpsilons[epsilonI];
            AssertPETSc(VecCopy(x, xPlus));
            AssertPETSc(VecAXPY(xPlus, epsilon, modeVector));
            AssertPETSc(VecCopy(x, xMinus));
            AssertPETSc(VecAXPY(xMinus, -epsilon, modeVector));
            AssertPETSc(petscSolid.formResidual(residualPlus, xPlus));

            scalarField plusFaceStabilisation;
            List<scalarField> plusBoundaryStabilisation;
            if (epsilonI == 0)
            {
                const surfaceScalarField& faceStabilisation =
                    petscSolid.mesh().lookupObject<surfaceScalarField>
                    (
                        "faceStabilisation(p)"
                    );
                plusFaceStabilisation = faceStabilisation.internalField();
                plusBoundaryStabilisation.setSize
                (
                    faceStabilisation.boundaryField().size()
                );
                forAll(faceStabilisation.boundaryField(), patchI)
                {
                    plusBoundaryStabilisation[patchI] =
                        faceStabilisation.boundaryField()[patchI];
                }
            }

            AssertPETSc(petscSolid.formResidual(residualMinus, xMinus));

            if (epsilonI == 0)
            {
                const surfaceScalarField& faceStabilisation =
                    petscSolid.mesh().lookupObject<surfaceScalarField>
                    (
                        "faceStabilisation(p)"
                    );
                const surfaceScalarField& gamma =
                    petscSolid.mesh().lookupObject<surfaceScalarField>
                    (
                        "rAUf"
                    );
                const surfaceScalarField& magSf = petscSolid.mesh().magSf();
                scalar internalFaceActionSqr = 0.0;
                scalar boundaryFaceActionSqr = 0.0;
                forAll(plusFaceStabilisation, faceI)
                {
                    const scalar action =
                        gamma[faceI]*magSf[faceI]
                       *(
                            plusFaceStabilisation[faceI]
                          - faceStabilisation[faceI]
                        )/(2.0*epsilon);
                    internalFaceActionSqr += sqr(action);
                }
                forAll(faceStabilisation.boundaryField(), patchI)
                {
                    forAll(faceStabilisation.boundaryField()[patchI], faceI)
                    {
                        const scalar action =
                            gamma.boundaryField()[patchI][faceI]
                           *magSf.boundaryField()[patchI][faceI]
                           *
                            (
                                plusBoundaryStabilisation[patchI][faceI]
                              - faceStabilisation.boundaryField()[patchI][faceI]
                            )/(2.0*epsilon);
                        boundaryFaceActionSqr += sqr(action);
                    }
                }
                reduce(internalFaceActionSqr, sumOp<scalar>());
                reduce(boundaryFaceActionSqr, sumOp<scalar>());
                Info<< "        physical stabilisation face-flux action: "
                    << "internal norm=" << Foam::sqrt(internalFaceActionSqr)
                    << ", boundary norm="
                    << Foam::sqrt(boundaryFaceActionSqr) << endl;
            }

            AssertPETSc(VecCopy(residualPlus, jv));
            AssertPETSc(VecAXPY(jv, -1.0, residualMinus));
            AssertPETSc(VecScale(jv, 1.0/(2.0*epsilon)));
            setBlockPart(jvBlock, jv, true);
            AssertPETSc(VecCopy(jvBlock, blockDifference));
            AssertPETSc(VecAXPY(blockDifference, -1.0, pvBlock));

            const scalar relativeError =
                norm(blockDifference)/(norm(jvBlock) + VSMALL);
            const scalar cosine = dotCosine(jvBlock, pvBlock);

            if (epsilonI == 0)
            {
                scalarField vHat;
                scalarField trueFinal;
                scalarField approximateFinal;
                pressureComponents(modeVector, vHat);
                pressureComponents(jvBlock, trueFinal);
                pressureComponents(pvBlock, approximateFinal);

                scalarField trueLocal(vHat.size(), 0.0);
                scalarField trueStabilisation(vHat.size(), 0.0);
                scalarField trueTotal(vHat.size(), 0.0);
                scalarField approximateLocal(vHat.size(), 0.0);
                scalarField approximateStabilisation(vHat.size(), 0.0);
                scalarField approximateTotal(vHat.size(), 0.0);
                scalarField finalRowMultiplier(vHat.size(), 0.0);

                forAll(vHat, cellI)
                {
                    const scalar rowMultiplier =
                        pressureEqnScale
                       *petscSolid.mesh().V()[cellI]
                       *
                        (
                            pressureRowScaling == "volumeRmsForce"
                          ? sqr(L0)/Foam::sqrt
                            (
                                Vtot*petscSolid.mesh().V()[cellI]
                            )
                          : 1.0
                        );
                    finalRowMultiplier[cellI] = rowMultiplier;

                    trueLocal[cellI] =
                       -pressureUnknownScale
                       *vHat[cellI]/tBulkModulus()[cellI];
                    approximateLocal[cellI] = trueLocal[cellI];
                    trueTotal[cellI] =
                        trueFinal[cellI]/rowMultiplier;
                    approximateTotal[cellI] =
                        approximateFinal[cellI]/rowMultiplier;
                    trueStabilisation[cellI] =
                        trueTotal[cellI] - trueLocal[cellI];
                    approximateStabilisation[cellI] =
                        approximateTotal[cellI] - approximateLocal[cellI];
                }

                auto reportScaledStage =
                [&]
                (
                    const word& name,
                    const scalarField& multiplier
                )
                {
                    scalarField tl(trueLocal*multiplier);
                    scalarField ts(trueStabilisation*multiplier);
                    scalarField tt(trueTotal*multiplier);
                    scalarField pl(approximateLocal*multiplier);
                    scalarField ps(approximateStabilisation*multiplier);
                    scalarField pt(approximateTotal*multiplier);
                    reportComponentStage(name, tl, ts, tt, pl, ps, pt);
                };

                scalarField beforeUnknown(vHat.size(), 1.0/pressureUnknownScale);
                scalarField afterUnknown(vHat.size(), 1.0);
                scalarField afterVolume(petscSolid.mesh().V());
                scalarField afterEqnScale
                (
                    pressureEqnScale*petscSolid.mesh().V()
                );
                reportScaledStage("intensive-before-pUnknown", beforeUnknown);
                reportScaledStage("intensive-after-pUnknown", afterUnknown);
                reportScaledStage("after-cell-volume", afterVolume);
                reportScaledStage("after-pressure-equation-scale", afterEqnScale);
                reportScaledStage("final-PETSc-row", finalRowMultiplier);

                if
                (
                    fieldNorm(trueStabilisation)
                 <= 1e-12*(fieldNorm(trueLocal) + VSMALL)
                )
                {
                    Info<< "        stabilisation action is a constant-mode "
                        << "null within 1e-12 of the local action" << endl;
                }
                else
                {
                    scalar maxTrueStabilisation = 0.0;
                    forAll(trueStabilisation, cellI)
                    {
                        maxTrueStabilisation = max
                        (
                            maxTrueStabilisation,
                            mag(trueStabilisation[cellI])
                        );
                    }
                    reduce(maxTrueStabilisation, maxOp<scalar>());
                    const scalar ratioThreshold =
                        max(1e-8*maxTrueStabilisation, VSMALL);
                    scalar minRatio = GREAT;
                    scalar maxRatio = -GREAT;
                    scalar sumRatio = 0.0;
                    label nRatio = 0;
                    label signMismatch = 0;
                    forAll(trueStabilisation, cellI)
                    {
                        if (mag(trueStabilisation[cellI]) > ratioThreshold)
                        {
                            const scalar ratio =
                                approximateStabilisation[cellI]
                               /trueStabilisation[cellI];
                            minRatio = min(minRatio, ratio);
                            maxRatio = max(maxRatio, ratio);
                            sumRatio += ratio;
                            ++nRatio;
                            if (ratio < 0.0)
                            {
                                ++signMismatch;
                            }
                        }
                    }
                    reduce(minRatio, minOp<scalar>());
                    reduce(maxRatio, maxOp<scalar>());
                    reduce(sumRatio, sumOp<scalar>());
                    reduce(nRatio, sumOp<label>());
                    reduce(signMismatch, sumOp<label>());
                    Info<< "        cellwise P_stab/J_stab ratio: count="
                        << nRatio << ", min=" << minRatio
                        << ", max=" << maxRatio << ", mean="
                        << sumRatio/max(nRatio, label(1))
                        << ", sign-mismatches=" << signMismatch << endl;
                }
            }

            // Required exact-pressure modes: constant, smooth, localised and
            // deterministic random. Additional pressure modes remain useful
            // diagnostics but do not define this focused regression gate.
            if (mode == 14 || mode == 6 || mode == 7 || mode == 15)
            {
                maximumRequiredPppRelativeError = max
                (
                    maximumRequiredPppRelativeError,
                    relativeError
                );
                minimumRequiredPppCosine = min
                (
                    minimumRequiredPppCosine,
                    cosine
                );
            }

            Info<< "    P_pp mode " << mode << " (" << modeName(mode)
                << "), epsilon = " << epsilon
                << ": ||J_pp v|| = " << norm(jvBlock)
                << ", ||P_pp v|| = " << norm(pvBlock)
                << ", relative error = " << relativeError
                << ", cosine = " << cosine
                << ", best-fit scalar = "
                << bestFitScalar(jvBlock, pvBlock)
                << ", scaled error = "
                << bestFitError
                   (
                       jvBlock,
                       pvBlock,
                       bestFitScalar(jvBlock, pvBlock)
                   ) << endl;

            maximumHistoryChange = max
            (
                maximumHistoryChange,
                maximumHistoryDifference
                (
                    petscSolid.mesh(), historyFieldNames, historyReference
                )
            );
        }
    }

    AssertPETSc(VecDestroy(&jvBlock));
    AssertPETSc(VecDestroy(&pvBlock));
    AssertPETSc(VecDestroy(&blockDifference));

    if (!pressureOnly)
    {
    // Stage 2 audit: decompose the complete momentum residual action at the
    // proven central-difference epsilon.  This section deliberately uses a
    // separate set of outputs and does not change any production switch.
    Vec jvD = nullptr;
    Vec pvD = nullptr;
    Vec termAction = nullptr;
    Vec termSum = nullptr;
    AssertPETSc(VecDuplicate(x, &jvD));
    AssertPETSc(VecDuplicate(x, &pvD));
    AssertPETSc(VecDuplicate(x, &termAction));
    AssertPETSc(VecDuplicate(x, &termSum));

    const scalar decompositionEpsilon = 1e-6;
    const label decompositionModes[] = {10, 11, 12, 13, 19, 18, 17};
    const label nDecompositionModes = sizeof(decompositionModes)/sizeof(label);

    Info<< "\nStage 2 term-by-term J_DD decomposition"
        << ", epsilon = " << decompositionEpsilon << nl
        << "    terms are intensive momentum residual actions before volume "
        << "scaling; PETSc comparisons are extensive" << endl;

    scalar maximumTermSumError = 0.0;
    for (label decompositionModeI = 0;
         decompositionModeI < nDecompositionModes;
         ++decompositionModeI)
    {
        const label mode = decompositionModes[decompositionModeI];
        AssertPETSc(VecSet(modeVector, 0.0));
        setMode(modeVector, mode, petscSolid.mesh());
        const scalar modeNorm = norm(modeVector);
        if (modeNorm > VSMALL)
        {
            AssertPETSc(VecScale(modeVector, 1.0/modeNorm));
        }
        AssertPETSc(MatMult(P, modeVector, pv));

        AssertPETSc(VecCopy(x, xPlus));
        AssertPETSc(VecAXPY(xPlus, decompositionEpsilon, modeVector));
        AssertPETSc(VecCopy(x, xMinus));
        AssertPETSc(VecAXPY(xMinus, -decompositionEpsilon, modeVector));

        solidModels::nonLinGeomTotalLagTotalDispSolid::
            momentumResidualDecompositionData baseDecomposition;
        solidModels::nonLinGeomTotalLagTotalDispSolid::
            momentumResidualDecompositionData plusDecomposition;
        solidModels::nonLinGeomTotalLagTotalDispSolid::
            momentumResidualDecompositionData minusDecomposition;

        AssertPETSc(petscSolid.formResidual(residual, x));
        petscSolid.momentumResidualDecomposition(baseDecomposition);
        AssertPETSc(petscSolid.formResidual(residualPlus, xPlus));
        petscSolid.momentumResidualDecomposition(plusDecomposition);
        AssertPETSc(petscSolid.formResidual(residualMinus, xMinus));
        petscSolid.momentumResidualDecomposition(minusDecomposition);

        AssertPETSc(VecCopy(residualPlus, jv));
        AssertPETSc(VecAXPY(jv, -1.0, residualMinus));
        AssertPETSc(VecScale(jv, 1.0/(2.0*decompositionEpsilon)));
        setDisplacementPart(jvD, jv);
        setDisplacementPart(pvD, pv);
        AssertPETSc(VecSet(termSum, 0.0));

        Info<< "  mode " << mode << " (" << modeName(mode) << ')' << nl;
        for (label termI = 0;
             termI < plusDecomposition.names.size();
             ++termI)
        {
            setCentralTerm
            (
                termAction,
                plusDecomposition,
                minusDecomposition,
                termI,
                decompositionEpsilon,
                petscSolid.mesh()
            );
            AssertPETSc(VecAXPY(termSum, 1.0, termAction));

            const scalar actionNorm = norm(termAction);
            const scalar faceNorm = centralFaceActionNorm
            (
                plusDecomposition,
                minusDecomposition,
                termI,
                decompositionEpsilon
            );
            scalar internalFaceNorm = 0.0;
            scalar boundaryFaceNorm = 0.0;
            centralFaceCategoryNorms
            (
                plusDecomposition,
                minusDecomposition,
                termI,
                decompositionEpsilon,
                petscSolid.mesh(),
                internalFaceNorm,
                boundaryFaceNorm
            );
            Info<< "    " << plusDecomposition.names[termI]
                << ": ||J_term v|| = " << actionNorm
                << ", fraction = "
                << actionNorm/(norm(jvD) + VSMALL)
                << ", cosine = " << dotCosine(termAction, jvD)
                << ", face-action norm = " << faceNorm
                << ", owner-cell face norm = " << internalFaceNorm
                << ", neighbour-cell face norm = " << internalFaceNorm
                << ", internal-face norm = " << internalFaceNorm
                << ", boundary-face norm = " << boundaryFaceNorm << nl;
        }

        AssertPETSc(VecCopy(termSum, difference));
        AssertPETSc(VecAXPY(difference, -1.0, jvD));
        const scalar termSumError = norm(difference)/(norm(jvD) + VSMALL);
        maximumTermSumError = max(maximumTermSumError, termSumError);

        const scalar fit = bestFitScalar(jvD, pvD);
        AssertPETSc(VecCopy(jvD, difference));
        AssertPETSc(VecAXPY(difference, -1.0, pvD));
        Info<< "    term-sum reconstruction relative error = "
            << termSumError << nl
            << "    ||J_DD v|| = " << norm(jvD) << nl
            << "    ||P_DD v|| = " << norm(pvD) << nl
            << "    P_DD relative action error = "
            << norm(difference)/(norm(jvD) + VSMALL) << nl
            << "    P_DD action cosine = " << dotCosine(jvD, pvD) << nl
            << "    current P_DD best-fit scalar = " << fit << nl
            << "    current P_DD error after best scalar = "
            << bestFitError(jvD, pvD, fit) << nl
            << "    direct/interpolated face-path action ratio = "
            << centralFaceActionNorm
               (
                   plusDecomposition,
                   minusDecomposition,
                   3,
                   decompositionEpsilon
               )/(centralFaceActionNorm
                  (
                      plusDecomposition,
                      minusDecomposition,
                      0,
                      decompositionEpsilon
                  )
                + centralFaceActionNorm
                  (
                      plusDecomposition,
                      minusDecomposition,
                      1,
                      decompositionEpsilon
                  )
                + centralFaceActionNorm
                  (
                      plusDecomposition,
                      minusDecomposition,
                      2,
                      decompositionEpsilon
                  )
                + VSMALL) << nl;

        reportGeometrySplit
        (
            baseDecomposition,
            plusDecomposition,
            minusDecomposition,
            decompositionEpsilon,
            petscSolid.mesh()
        );
        reportLargestCellsAndFaces
        (
            jvD,
            pvD,
            plusDecomposition,
            minusDecomposition,
            decompositionEpsilon,
            petscSolid.mesh()
        );

        const scalar historyChange = maximumHistoryDifference
        (
            petscSolid.mesh(), historyFieldNames, historyReference
        );
        maximumHistoryChange = max(maximumHistoryChange, historyChange);
    }

    Info<< "Maximum Stage 2 term-sum reconstruction relative error = "
        << maximumTermSumError << endl;

    AssertPETSc(VecDestroy(&jvD));
    AssertPETSc(VecDestroy(&pvD));
    AssertPETSc(VecDestroy(&termAction));
    AssertPETSc(VecDestroy(&termSum));
    }

    Info<< "Maximum accepted-history checksum change = "
        << maximumHistoryChange << nl
        << "Maximum repeated-baseline residual difference = "
        << maximumResidualRepeatDifference << nl
        << "Maximum required P_pp relative action error = "
        << maximumRequiredPppRelativeError << nl
        << "Minimum required P_pp action cosine = "
        << minimumRequiredPppCosine << endl;

    if
    (
        pressureOnly
     &&
        (
            maximumRequiredPppRelativeError > 1e-8
         || minimumRequiredPppCosine < 0.99999999
         || maximumHistoryChange != 0.0
         || maximumResidualRepeatDifference != 0.0
        )
    )
    {
        FatalErrorInFunction
            << "Exact P_pp action regression failed: relative error "
            << maximumRequiredPppRelativeError << ", cosine "
            << minimumRequiredPppCosine << ", history change "
            << maximumHistoryChange << ", repeated residual difference "
            << maximumResidualRepeatDifference << abort(FatalError);
    }

    AssertPETSc(MatDestroy(&P));
    AssertPETSc(VecDestroy(&x));
    AssertPETSc(VecDestroy(&xPlus));
    AssertPETSc(VecDestroy(&xMinus));
    AssertPETSc(VecDestroy(&residual));
    AssertPETSc(VecDestroy(&residualRepeat));
    AssertPETSc(VecDestroy(&residualPlus));
    AssertPETSc(VecDestroy(&residualMinus));
    AssertPETSc(VecDestroy(&jv));
    AssertPETSc(VecDestroy(&pv));
    AssertPETSc(VecDestroy(&difference));
    AssertPETSc(VecDestroy(&modeVector));

    Info<< "PASS: Jv/Pv action audit completed" << endl;
#endif
    return 0;
}

// ************************************************************************* //
