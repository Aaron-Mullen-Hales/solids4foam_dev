/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

Application
    ArosticaDirectFaceMixedPressure

Description
    Verifies physical PETSc pressure unpacking, direct-face mixed-pressure
    completion, pressure linearity, Cauchy/nominal force equivalence, and
    replacement of constitutive force on solidTraction patches.

\*---------------------------------------------------------------------------*/

#ifdef USE_PETSC
    #include "petscErrorHandling.H"
#endif

#include "fvCFD.H"
#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "solidModel.H"
#include "solidTractionFvPatchVectorField.H"

#include <cmath>
#include <limits>

using namespace Foam;


namespace
{

typedef solidModels::nonLinGeomTotalLagTotalDispSolid Model;
typedef Model::momentumResidualDecompositionData Decomposition;


template<class Type>
void flattenSurfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    Field<Type>& flat
)
{
    const fvMesh& mesh = field.mesh();
    flat.setSize(mesh.nFaces());

    forAll(field, faceI)
    {
        flat[faceI] = field[faceI];
    }

    forAll(field.boundaryField(), patchI)
    {
        const Field<Type>& patchField = field.boundaryField()[patchI];
        forAll(patchField, faceI)
        {
            flat[mesh.boundary()[patchI].start() + faceI] =
                patchField[faceI];
        }
    }
}


void selectedStress
(
    const Decomposition& decomposition,
    symmTensorField& stress
)
{
    stress.setSize(decomposition.passiveFaceStress.size());

    forAll(stress, faceI)
    {
        stress[faceI] =
            decomposition.passiveFaceStress[faceI]
          + decomposition.viscousFaceStress[faceI]
          + decomposition.pressureFaceStress[faceI];
    }
}


scalar maxMagnitude(const symmTensorField& field)
{
    scalar result = 0.0;
    forAll(field, faceI)
    {
        result = max(result, mag(field[faceI]));
    }
    reduce(result, maxOp<scalar>());
    return result;
}


void evaluateAtPressure
(
    Model& model,
    volScalarField& p,
    const scalarField& pressureShape,
    const scalar multiplier,
    symmTensorField& stress,
    scalarField* pressureFace,
    Decomposition* retainedDecomposition
)
{
    forAll(p, cellI)
    {
        p[cellI] = multiplier*pressureShape[cellI];
    }
    p.correctBoundaryConditions();

    if (p.dimensions() != dimPressure)
    {
        FatalErrorInFunction
            << "The unpacked/updated p field does not have pressure dimensions: "
            << p.dimensions() << abort(FatalError);
    }

    Decomposition decomposition;
    model.momentumResidualDecomposition(decomposition);
    selectedStress(decomposition, stress);

    if (pressureFace)
    {
        const surfaceScalarField pFace(fvc::interpolate(p));
        flattenSurfaceField(pFace, *pressureFace);
    }

    if (retainedDecomposition)
    {
        *retainedDecomposition = decomposition;
    }
}

}


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "expectOmission",
        "Require and report the pre-correction direct-face pressure omission"
    );
    argList::addBoolOption
    (
        "checkInterpolatedStateOnly",
        "Check pressure after production residual formation and exit"
    );
    argList::addOption
    (
        "pressureScaleMode",
        "word",
        "Expected pressure unknown scale: twoMu or user"
    );
    argList::addOption
    (
        "expectedPressureScale",
        "scalar",
        "Expected test-only user pressure unknown scale"
    );
    argList::addOption
    (
        "expectedSolidType",
        "word",
        "Expected solid-model runtime type"
    );

    #include "setRootCase.H"
    #include "createTime.H"

#ifndef USE_PETSC
    FatalErrorInFunction
        << "ArosticaDirectFaceMixedPressure requires USE_PETSC"
        << abort(FatalError);
    return 1;
#else
    const bool expectOmission = args.optionFound("expectOmission");
    const bool checkInterpolatedStateOnly =
        args.optionFound("checkInterpolatedStateOnly");
    word pressureScaleMode("twoMu");
    args.optionReadIfPresent("pressureScaleMode", pressureScaleMode);
    word expectedSolidType
    (
        "arosticaNonLinearGeometryTotalLagrangianTotalDisplacement"
    );
    args.optionReadIfPresent("expectedSolidType", expectedSolidType);

    autoPtr<solidModel> solid =
        solidModel::New(runTime, dynamicFvMesh::defaultRegion);
    Model& model = refCast<Model>(*solid);
    fvMesh& mesh = model.mesh();

    if
    (
        model.type()
     != expectedSolidType
    )
    {
        FatalErrorInFunction
            << "Unexpected solid-model runtime type " << model.type()
            << abort(FatalError);
    }

    volScalarField& p = mesh.lookupObjectRef<volScalarField>("p");
    if (p.dimensions() != dimPressure)
    {
        FatalErrorInFunction
            << "p must have dimPressure, found " << p.dimensions()
            << abort(FatalError);
    }

    scalar expectedPressureScale = 0.0;
    if (pressureScaleMode == "twoMu")
    {
        const volScalarField twoMu(2.0*model.mechanical().shearModulus());
        scalar weightedTwoMu = 0.0;
        scalar volume = 0.0;
        forAll(twoMu, cellI)
        {
            weightedTwoMu += twoMu[cellI]*mesh.V()[cellI];
            volume += mesh.V()[cellI];
        }
        reduce(weightedTwoMu, sumOp<scalar>());
        reduce(volume, sumOp<scalar>());
        expectedPressureScale = weightedTwoMu/volume;
    }
    else if (pressureScaleMode == "user")
    {
        if (!args.optionReadIfPresent("expectedPressureScale", expectedPressureScale))
        {
            FatalErrorInFunction
                << "-expectedPressureScale is required for user mode"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown -pressureScaleMode " << pressureScaleMode
            << "; expected twoMu or user" << abort(FatalError);
    }

    Vec x = nullptr;
    Vec residual = nullptr;
    AssertPETSc(model.initialiseSolution(x));
    AssertPETSc(VecDuplicate(x, &residual));
    AssertPETSc(VecSet(x, 0.0));

    scalarField pHat(mesh.nCells());
    forAll(pHat, cellI)
    {
        const point& centre = mesh.C()[cellI];
        pHat[cellI] =
            0.25 + 0.05*centre.x() + 0.03*centre.y() + 0.02*centre.z();
    }
    model.InsertFieldComponents<scalar>(pHat, x, 3);
    AssertPETSc(model.formResidual(residual, x));

    scalar pressureUnpackError = 0.0;
    scalar pressureUnpackScale = 1.0;
    forAll(p, cellI)
    {
        const scalar expected = expectedPressureScale*pHat[cellI];
        pressureUnpackError =
            max(pressureUnpackError, mag(p[cellI] - expected));
        pressureUnpackScale =
            max(pressureUnpackScale, max(mag(p[cellI]), mag(expected)));
    }
    reduce(pressureUnpackError, maxOp<scalar>());
    reduce(pressureUnpackScale, maxOp<scalar>());

    const scalar machineEpsilon = std::numeric_limits<scalar>::epsilon();
    const scalar pressureUnpackTolerance =
        1.0e-12 + 100.0*machineEpsilon*pressureUnpackScale;

    if
    (
        pressureUnpackError > pressureUnpackTolerance
     || p.dimensions() != dimPressure
    )
    {
        FatalErrorInFunction
            << "Physical pressure unpacking failed: error="
            << pressureUnpackError << ", tolerance="
            << pressureUnpackTolerance << ", dimensions=" << p.dimensions()
            << abort(FatalError);
    }

    if (checkInterpolatedStateOnly)
    {
        Decomposition decomposition;
        model.momentumResidualDecomposition(decomposition);

        scalarField pressureFace;
        const surfaceScalarField pFace(fvc::interpolate(p));
        flattenSurfaceField(pFace, pressureFace);

        symmTensorField interpolationError
        (
            decomposition.pressureFaceStress.size()
        );
        symmTensorField pressureIdentity(interpolationError.size());
        forAll(interpolationError, faceI)
        {
            pressureIdentity[faceI] =
                pressureFace[faceI]*symmTensor::I;
            interpolationError[faceI] =
                decomposition.pressureFaceStress[faceI]
              + pressureIdentity[faceI];
        }

        const scalar error = maxMagnitude(interpolationError);
        const scalar scale =
            max(1.0, maxMagnitude(pressureIdentity));
        const scalar tolerance =
            1.0e-12 + 100.0*machineEpsilon*scale;
        if (error > tolerance)
        {
            FatalErrorInFunction
                << "Interpolated-cell pressure stress failed: error="
                << error << ", tolerance=" << tolerance
                << abort(FatalError);
        }

        AssertPETSc(VecDestroy(&residual));
        AssertPETSc(VecDestroy(&x));
        Info<< "INTERPOLATED_CELL_RESULT"
            << " E_sigma=" << error
            << " tolerance=" << tolerance << nl
            << "PASS: generic interpolatedCell production residual path"
            << endl;
        return 0;
    }

    AssertPETSc(VecDestroy(&residual));
    AssertPETSc(VecDestroy(&x));

    scalarField pressureShape(mesh.nCells());
    forAll(pressureShape, cellI)
    {
        const point& centre = mesh.C()[cellI];
        pressureShape[cellI] =
            1000.0
          + 250.0*(centre.x() + 2.0*centre.y() + 3.0*centre.z());
    }

    symmTensorField sigmaZero;
    symmTensorField sigmaP;
    symmTensorField sigmaTwoP;
    scalarField pFace;
    Decomposition pressureDecomposition;
    evaluateAtPressure
    (
        model, p, pressureShape, 0.0, sigmaZero, nullptr, nullptr
    );
    evaluateAtPressure
    (
        model,
        p,
        pressureShape,
        1.0,
        sigmaP,
        &pFace,
        &pressureDecomposition
    );
    evaluateAtPressure
    (
        model, p, pressureShape, 2.0, sigmaTwoP, nullptr, nullptr
    );

    symmTensorField deltaP(sigmaP.size());
    symmTensorField deltaTwoP(sigmaP.size());
    symmTensorField completionError(sigmaP.size());
    symmTensorField linearityError(sigmaP.size());
    symmTensorField pressureIdentity(sigmaP.size());

    forAll(sigmaP, faceI)
    {
        deltaP[faceI] = sigmaP[faceI] - sigmaZero[faceI];
        deltaTwoP[faceI] = sigmaTwoP[faceI] - sigmaZero[faceI];
        pressureIdentity[faceI] = pFace[faceI]*symmTensor::I;
        completionError[faceI] = deltaP[faceI] + pressureIdentity[faceI];
        linearityError[faceI] = deltaTwoP[faceI] - 2.0*deltaP[faceI];
    }

    const scalar completionNorm = maxMagnitude(completionError);
    const scalar deltaNorm = maxMagnitude(deltaP);
    const scalar twoDeltaNorm = maxMagnitude(deltaTwoP);
    const scalar pressureIdentityNorm = maxMagnitude(pressureIdentity);
    const scalar completionScale =
        max(1.0, max(deltaNorm, pressureIdentityNorm));
    const scalar completionTolerance =
        1.0e-12 + 100.0*machineEpsilon*completionScale;
    const scalar linearityNorm = maxMagnitude(linearityError);
    const scalar linearityScale =
        max(1.0, max(twoDeltaNorm, 2.0*deltaNorm));
    const scalar linearityTolerance =
        1.0e-12 + 100.0*machineEpsilon*linearityScale;

    if (expectOmission)
    {
        if
        (
            deltaNorm > completionTolerance
         || completionNorm < 0.99*pressureIdentityNorm
        )
        {
            FatalErrorInFunction
                << "Pre-correction omission was not demonstrated: "
                << "|deltaSigma|=" << deltaNorm
                << ", E_sigma=" << completionNorm
                << ", |p_f I|=" << pressureIdentityNorm
                << abort(FatalError);
        }
    }
    else if
    (
        completionNorm > completionTolerance
     || linearityNorm > linearityTolerance
    )
    {
        FatalErrorInFunction
            << "Direct-face pressure completion failed: E_sigma="
            << completionNorm << " (tol " << completionTolerance
            << "), E_2p=" << linearityNorm << " (tol "
            << linearityTolerance << ')' << abort(FatalError);
    }

    tensorField faceF;
    vectorField referenceSf;
    flattenSurfaceField
    (
        mesh.lookupObject<surfaceTensorField>("Ff_"),
        faceF
    );
    flattenSurfaceField(mesh.Sf(), referenceSf);

    scalar forceError = 0.0;
    scalar forceScale = 1.0;
    forAll(sigmaP, faceI)
    {
        const vector currentForce =
            pressureDecomposition.faceAreaVectors[faceI] & sigmaP[faceI];
        const tensor nominalStress =
            det(faceF[faceI])
           *(
                tensor(sigmaP[faceI])
              & inv(faceF[faceI]).T()
            );
        // P*N in index notation.  OpenFOAM's Tensor&Vector operator performs
        // this product directly; Vector&Tensor would apply P^T instead.
        const vector nominalForce = nominalStress & referenceSf[faceI];
        forceError = max(forceError, mag(currentForce - nominalForce));
        forceScale =
            max(forceScale, max(mag(currentForce), mag(nominalForce)));
    }
    reduce(forceError, maxOp<scalar>());
    reduce(forceScale, maxOp<scalar>());
    const scalar forceTolerance = 1.0e-12 + 1.0e-9*forceScale;

    if (forceError > forceTolerance)
    {
        FatalErrorInFunction
            << "Cauchy/current-area and nominal-force forms differ: E_F="
            << forceError << ", tolerance=" << forceTolerance
            << abort(FatalError);
    }

    scalar tractionBoundaryError = 0.0;
    scalar tractionBoundaryScale = 1.0;
    label tractionBoundaryFaces = 0;
    const volVectorField& D = model.solutionD();
    forAll(D.boundaryField(), patchI)
    {
        if (!isA<solidTractionFvPatchVectorField>(D.boundaryField()[patchI]))
        {
            continue;
        }

        const solidTractionFvPatchVectorField& tractionPatch =
            refCast<const solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            );
        forAll(D.boundaryField()[patchI], patchFaceI)
        {
            const label faceI =
                mesh.boundary()[patchI].start() + patchFaceI;
            const scalar selectedArea =
                tractionPatch.useUndeformedArea()
              ? mesh.boundary()[patchI].magSf()[patchFaceI]
              : mag(pressureDecomposition.faceAreaVectors[faceI]);
            const vector nCurrent =
                pressureDecomposition.faceAreaVectors[faceI]
               /(mag(pressureDecomposition.faceAreaVectors[faceI]) + VSMALL);
            const vector expectedForce =
                selectedArea
               *(
                    tractionPatch.traction()[patchFaceI]
                  - tractionPatch.pressure()[patchFaceI]*nCurrent
                );
            vector actualForce(vector::zero);
            for (label termI = 0; termI < 9; ++termI)
            {
                actualForce +=
                    pressureDecomposition.faceForces[termI][faceI];
            }
            tractionBoundaryError =
                max(tractionBoundaryError, mag(actualForce - expectedForce));
            tractionBoundaryScale =
                max
                (
                    tractionBoundaryScale,
                    max(mag(actualForce), mag(expectedForce))
                );
            ++tractionBoundaryFaces;
        }
    }
    reduce(tractionBoundaryError, maxOp<scalar>());
    reduce(tractionBoundaryScale, maxOp<scalar>());
    reduce(tractionBoundaryFaces, sumOp<label>());
    const scalar tractionBoundaryTolerance =
        1.0e-12 + 1.0e-9*tractionBoundaryScale;

    if
    (
        tractionBoundaryFaces == 0
     || tractionBoundaryError > tractionBoundaryTolerance
    )
    {
        FatalErrorInFunction
            << "solidTraction replacement check failed: faces="
            << tractionBoundaryFaces << ", error=" << tractionBoundaryError
            << ", tolerance=" << tractionBoundaryTolerance
            << abort(FatalError);
    }

    Info<< "PRESSURE_PHYSICAL_FIELD_PROOF"
        << " mode=" << pressureScaleMode
        << " expectedScale=" << expectedPressureScale
        << " E_p=" << pressureUnpackError
        << " tolerance=" << pressureUnpackTolerance
        << " dimensions=" << p.dimensions() << nl
        << "DIRECT_FACE_PRESSURE_RESULT"
        << " expected=" << (expectOmission ? "omission" : "complete")
        << " deltaSigma=" << deltaNorm
        << " E_sigma=" << completionNorm
        << " tolerance=" << completionTolerance
        << " pressureIdentity=" << pressureIdentityNorm
        << " E_2p=" << linearityNorm
        << " E_2p_tolerance=" << linearityTolerance << nl
        << "FORCE_EQUIVALENCE"
        << " E_F=" << forceError
        << " tolerance=" << forceTolerance << nl
        << "TRACTION_REPLACEMENT"
        << " faces=" << tractionBoundaryFaces
        << " E_F_BC=" << tractionBoundaryError
        << " tolerance=" << tractionBoundaryTolerance << nl
        << "PASS: Aróstica direct-face mixed-pressure regression" << endl;

    return 0;
#endif
}


// ************************************************************************* //
