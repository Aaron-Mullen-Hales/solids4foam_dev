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
#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "solidModel.H"
#include "solidTractionFvPatchVectorField.H"

#include <cmath>

using namespace Foam;


#ifdef USE_PETSC
namespace
{

struct trialSnapshot
{
    vectorField displacement;
    vectorField velocity;
    vectorField traction;
};


trialSnapshot captureTrial
(
    const volVectorField& D,
    const label patchI
)
{
    const solidTractionFvPatchVectorField& patch =
        refCast<const solidTractionFvPatchVectorField>
        (
            D.boundaryField()[patchI]
        );
    const tmp<volVectorField> tDdot(fvc::ddt(D));

    trialSnapshot result
    {
        vectorField(D.boundaryField()[patchI]),
        vectorField(tDdot().boundaryField()[patchI]),
        vectorField(patch.traction())
    };

    return result;
}


scalar maximumDifference
(
    const trialSnapshot& a,
    const trialSnapshot& b
)
{
    scalar result = 0.0;

    forAll(a.displacement, faceI)
    {
        result = max(result, mag(a.displacement[faceI] - b.displacement[faceI]));
        result = max(result, mag(a.velocity[faceI] - b.velocity[faceI]));
        result = max(result, mag(a.traction[faceI] - b.traction[faceI]));
    }

    return result;
}


scalar maximumDifference(Vec a, Vec b)
{
    PetscInt size = 0;
    const PetscScalar* aValues = nullptr;
    const PetscScalar* bValues = nullptr;
    AssertPETSc(VecGetLocalSize(a, &size));
    AssertPETSc(VecGetArrayRead(a, &aValues));
    AssertPETSc(VecGetArrayRead(b, &bValues));

    scalar result = 0.0;
    for (PetscInt i = 0; i < size; ++i)
    {
        result = max(result, mag(PetscRealPart(aValues[i] - bValues[i])));
    }

    AssertPETSc(VecRestoreArrayRead(a, &aValues));
    AssertPETSc(VecRestoreArrayRead(b, &bValues));
    return result;
}


scalar maximumDifference
(
    const volSymmTensorField& a,
    const volSymmTensorField& b
)
{
    scalar result = 0.0;

    forAll(a, cellI)
    {
        result = max(result, mag(a[cellI] - b[cellI]));
    }

    forAll(a.boundaryField(), patchI)
    {
        const symmTensorField& aPatch = a.boundaryField()[patchI];
        const symmTensorField& bPatch = b.boundaryField()[patchI];
        forAll(aPatch, faceI)
        {
            result = max(result, mag(aPatch[faceI] - bPatch[faceI]));
        }
    }

    return result;
}


scalar maximumDifference
(
    const surfaceSymmTensorField& a,
    const surfaceSymmTensorField& b
)
{
    scalar result = 0.0;

    forAll(a, faceI)
    {
        result = max(result, mag(a[faceI] - b[faceI]));
    }

    forAll(a.boundaryField(), patchI)
    {
        const symmTensorField& aPatch = a.boundaryField()[patchI];
        const symmTensorField& bPatch = b.boundaryField()[patchI];
        forAll(aPatch, faceI)
        {
            result = max(result, mag(aPatch[faceI] - bPatch[faceI]));
        }
    }

    return result;
}


void perturbDisplacement(Vec x, const scalar scale = 1.0)
{
    PetscInt size = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(x, &size));
    AssertPETSc(VecGetBlockSize(x, &blockSize));
    AssertPETSc(VecGetArray(x, &values));

    if (blockSize < 3)
    {
        FatalErrorInFunction
            << "Unexpected PETSc block size " << blockSize
            << abort(FatalError);
    }

    for (PetscInt i = 0; i < size; i += blockSize)
    {
        values[i] += scale*0.017;
        values[i + 1] -= scale*0.013;
        values[i + 2] += scale*0.029;
    }

    AssertPETSc(VecRestoreArray(x, &values));
}


scalar vectorNorm(Vec x)
{
    PetscReal result = 0.0;
    AssertPETSc(VecNorm(x, NORM_2, &result));
    return result;
}

}
#endif




int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    autoPtr<solidModel> solid =
        solidModel::New(runTime, dynamicFvMesh::defaultRegion);

    volVectorField& D = solid->D();
    D.correctBoundaryConditions();

    label tractionPatchCount = 0;
    label tractionPatchI = -1;

    forAll(D.boundaryField(), patchI)
    {
        if (isA<solidTractionFvPatchVectorField>(D.boundaryField()[patchI]))
        {
            ++tractionPatchCount;
            tractionPatchI = patchI;

            const solidTractionFvPatchVectorField& bc =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            if (!bc.useUndeformedArea() || max(mag(bc.pressure())) > SMALL)
            {
                FatalErrorInFunction
                    << "Aróstica boundary has invalid area or pressure state"
                    << abort(FatalError);
            }
        }
    }

    if (tractionPatchCount != 1)
    {
        FatalErrorInFunction
            << "Expected one Aróstica traction patch, found "
            << tractionPatchCount << abort(FatalError);
    }

    const solidTractionFvPatchVectorField& tractionPatch =
        refCast<const solidTractionFvPatchVectorField>
        (
            D.boundaryField()[tractionPatchI]
        );
    const vectorField patchDisplacementA
    (
        D.boundaryField()[tractionPatchI]
    );
    const vectorField tractionA(tractionPatch.traction());
    const tmp<volVectorField> tDdotA(fvc::ddt(D));
    const vectorField patchVelocityA
    (
        tDdotA().boundaryField()[tractionPatchI]
    );

    D.primitiveFieldRef() = vector(0.017, -0.013, 0.029);
    D.correctBoundaryConditions();
    const vectorField patchDisplacementB
    (
        D.boundaryField()[tractionPatchI]
    );
    const solidTractionFvPatchVectorField& tractionPatchB =
        refCast<const solidTractionFvPatchVectorField>
        (
            D.boundaryField()[tractionPatchI]
        );
    const vectorField tractionB(tractionPatchB.traction());
    const tmp<volVectorField> tDdotB(fvc::ddt(D));
    const vectorField patchVelocityB
    (
        tDdotB().boundaryField()[tractionPatchI]
    );

    D.primitiveFieldRef() = vector::zero;
    D.correctBoundaryConditions();
    const solidTractionFvPatchVectorField& tractionPatchAagain =
        refCast<const solidTractionFvPatchVectorField>
        (
            D.boundaryField()[tractionPatchI]
        );
    const vectorField patchDisplacementAagain
    (
        D.boundaryField()[tractionPatchI]
    );
    const vectorField tractionAagain(tractionPatchAagain.traction());
    const tmp<volVectorField> tDdotAagain(fvc::ddt(D));
    const vectorField patchVelocityAagain
    (
        tDdotAagain().boundaryField()[tractionPatchI]
    );

    scalar maximumTrialDifference = 0.0;
    scalar maximumTrialSensitivity = 0.0;
    forAll(patchDisplacementA, faceI)
    {
        maximumTrialDifference =
            max
            (
                maximumTrialDifference,
                mag(patchDisplacementA[faceI] - patchDisplacementAagain[faceI])
            );
        maximumTrialDifference =
            max
            (
                maximumTrialDifference,
                mag(patchVelocityA[faceI] - patchVelocityAagain[faceI])
            );
        maximumTrialDifference =
            max
            (
                maximumTrialDifference,
                mag(tractionA[faceI] - tractionAagain[faceI])
            );
        maximumTrialSensitivity =
            max
            (
                maximumTrialSensitivity,
                mag(patchDisplacementA[faceI] - patchDisplacementB[faceI])
            );
        maximumTrialSensitivity =
            max
            (
                maximumTrialSensitivity,
                mag(patchVelocityA[faceI] - patchVelocityB[faceI])
            );
        maximumTrialSensitivity =
            max
            (
                maximumTrialSensitivity,
                mag(tractionA[faceI] - tractionB[faceI])
            );
    }

    if (maximumTrialDifference > 1e-12 || maximumTrialSensitivity <= 1e-12)
    {
        FatalErrorInFunction
            << "Trial boundary sequence failed: maximum repeat difference = "
            << maximumTrialDifference << ", maximum sensitivity = "
            << maximumTrialSensitivity << abort(FatalError);
    }

    Info<< "PASS: first-evaluation boundary trial sequence" << nl
        << "    maximum repeated-state difference = "
        << maximumTrialDifference << nl
        << "    maximum A/B sensitivity = " << maximumTrialSensitivity
        << endl;

#ifdef USE_PETSC
    solidModels::nonLinGeomTotalLagTotalDispSolid& petscSolid =
        refCast<solidModels::nonLinGeomTotalLagTotalDispSolid>(*solid);

    Vec xA = nullptr;
    Vec xB = nullptr;
    Vec fA = nullptr;
    Vec fB = nullptr;
    Vec fAagain = nullptr;
    Vec fBagain = nullptr;
    Vec xEps = nullptr;
    Vec fEps = nullptr;
    Vec jv = nullptr;
    Vec jvReference = nullptr;
    AssertPETSc(petscSolid.initialiseSolution(xA));
    AssertPETSc(VecDuplicate(xA, &xB));
    AssertPETSc(VecDuplicate(xA, &fA));
    AssertPETSc(VecDuplicate(xA, &fB));
    AssertPETSc(VecDuplicate(xA, &fAagain));
    AssertPETSc(VecDuplicate(xA, &fBagain));
    AssertPETSc(VecDuplicate(xA, &xEps));
    AssertPETSc(VecDuplicate(xA, &fEps));
    AssertPETSc(VecDuplicate(xA, &jv));
    AssertPETSc(VecDuplicate(xA, &jvReference));

    const fvMesh& mesh = solid->mesh();
    const bool hasViscoHistory =
        mesh.foundObject<volSymmTensorField>
        (
            "ArosticaEOld_myocardium"
        );
    scalar maximumProductionHistoryChange = 0.0;
    scalar maximumMffdHistoryChange = 0.0;

    autoPtr<volSymmTensorField> eOldSnapshotPtr;
    autoPtr<volSymmTensorField> eOldOldSnapshotPtr;
    autoPtr<surfaceSymmTensorField> efOldSnapshotPtr;
    autoPtr<surfaceSymmTensorField> efOldOldSnapshotPtr;

    const volSymmTensorField* eOldPtr = nullptr;
    const volSymmTensorField* eOldOldPtr = nullptr;
    const surfaceSymmTensorField* efOldPtr = nullptr;
    const surfaceSymmTensorField* efOldOldPtr = nullptr;

    if (hasViscoHistory)
    {
        eOldPtr = &mesh.lookupObject<volSymmTensorField>
        (
            "ArosticaEOld_myocardium"
        );
        eOldOldPtr = &mesh.lookupObject<volSymmTensorField>
        (
            "ArosticaEOldOld_myocardium"
        );
        efOldPtr = &mesh.lookupObject<surfaceSymmTensorField>
        (
            "ArosticaEfOld_myocardium"
        );
        efOldOldPtr = &mesh.lookupObject<surfaceSymmTensorField>
        (
            "ArosticaEfOldOld_myocardium"
        );

        eOldSnapshotPtr.reset(new volSymmTensorField(*eOldPtr));
        eOldOldSnapshotPtr.reset(new volSymmTensorField(*eOldOldPtr));
        efOldSnapshotPtr.reset(new surfaceSymmTensorField(*efOldPtr));
        efOldOldSnapshotPtr.reset(new surfaceSymmTensorField(*efOldOldPtr));
    }

    auto updateHistoryDifference = [&]()
    {
        if (hasViscoHistory)
        {
            maximumProductionHistoryChange = max
            (
                maximumProductionHistoryChange,
                max
                (
                    maximumDifference
                    (
                        *eOldPtr,
                        *eOldSnapshotPtr
                    ),
                    maximumDifference
                    (
                        *eOldOldPtr,
                        *eOldOldSnapshotPtr
                    )
                )
            );
            maximumProductionHistoryChange = max
            (
                maximumProductionHistoryChange,
                max
                (
                    maximumDifference
                    (
                        *efOldPtr,
                        *efOldSnapshotPtr
                    ),
                    maximumDifference
                    (
                        *efOldOldPtr,
                        *efOldOldSnapshotPtr
                    )
                )
            );
        }
    };

    AssertPETSc(petscSolid.formResidual(fA, xA));
    updateHistoryDifference();
    const trialSnapshot productionA = captureTrial(D, tractionPatchI);
    AssertPETSc(VecCopy(xA, xB));
    perturbDisplacement(xB);
    AssertPETSc(petscSolid.formResidual(fB, xB));
    updateHistoryDifference();
    const trialSnapshot productionB = captureTrial(D, tractionPatchI);
    AssertPETSc(petscSolid.formResidual(fAagain, xA));
    updateHistoryDifference();
    const trialSnapshot productionAagain = captureTrial(D, tractionPatchI);
    AssertPETSc(petscSolid.formResidual(fBagain, xB));
    updateHistoryDifference();
    const trialSnapshot productionBagain = captureTrial(D, tractionPatchI);

    const scalar productionTrialDifference =
        max
        (
            max
            (
                maximumDifference(productionA, productionAagain),
                maximumDifference(fA, fAagain)
            ),
            max
            (
                maximumDifference(productionB, productionBagain),
                maximumDifference(fB, fBagain)
            )
        );
    const scalar productionSensitivity =
        max
        (
            maximumDifference(productionA, productionB),
            maximumDifference(fA, fB)
        );

    if
    (
        productionTrialDifference > 1e-12
     || productionSensitivity <= 1e-12
    )
    {
        FatalErrorInFunction
            << "PETSc production trial sequence failed: repeat difference = "
            << productionTrialDifference << ", sensitivity = "
            << productionSensitivity << abort(FatalError);
    }

    Info<< "PASS: PETSc production trial sequence" << nl
        << "    maximum repeated-state difference = "
        << productionTrialDifference << nl
        << "    maximum residual/traction sensitivity = "
        << productionSensitivity << nl
        << "    maximum accepted-history change = "
        << maximumProductionHistoryChange << endl;

    if
    (
        hasViscoHistory
     && maximumProductionHistoryChange > 1e-12
    )
    {
        FatalErrorInFunction
            << "Production residual changed accepted viscoelastic history; "
            << "maximum difference = " << maximumProductionHistoryChange
            << abort(FatalError);
    }

    // PETSc MFFD-style finite differences over a smooth logistic constitutive
    // state.  The smallest perturbation is retained as the reference and
    // every evaluation is checked against the same accepted history.
    const scalarField epsilons
    (
        List<scalar>({1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2})
    );
    scalar referenceNorm = 0.0;
    scalar maximumMffdReferenceDifference = 0.0;

    Info<< "PETSc production MFFD-style sweep" << nl;
    forAll(epsilons, epsI)
    {
        const scalar epsilon = epsilons[epsI];
        AssertPETSc(VecCopy(xA, xEps));
        perturbDisplacement(xEps, epsilon);
        AssertPETSc(petscSolid.formResidual(fEps, xEps));
        updateHistoryDifference();
        maximumMffdHistoryChange = max
        (
            maximumMffdHistoryChange,
            maximumProductionHistoryChange
        );
        AssertPETSc(VecCopy(fEps, jv));
        AssertPETSc(VecAXPY(jv, -1.0, fA));
        AssertPETSc(VecScale(jv, 1.0/epsilon));

        const scalar jvNorm = vectorNorm(jv);
        if (epsI == 0)
        {
            referenceNorm = jvNorm;
            AssertPETSc(VecCopy(jv, jvReference));
        }
        else
        {
            AssertPETSc(VecCopy(jv, fAagain));
            AssertPETSc(VecAXPY(fAagain, -1.0, jvReference));
            maximumMffdReferenceDifference = max
            (
                maximumMffdReferenceDifference,
                vectorNorm(fAagain)
            );
        }

        Info<< "    epsilon = " << epsilon
            << ", ||Jv_FD||2 = " << jvNorm << nl;
    }

    Info<< "    reference ||Jv||2 = " << referenceNorm << nl
        << "    maximum coarse/reference difference = "
        << maximumMffdReferenceDifference << nl
        << "    maximum MFFD history change = "
        << maximumMffdHistoryChange << endl;

    AssertPETSc(VecDestroy(&xA));
    AssertPETSc(VecDestroy(&xB));
    AssertPETSc(VecDestroy(&fA));
    AssertPETSc(VecDestroy(&fB));
    AssertPETSc(VecDestroy(&fAagain));
    AssertPETSc(VecDestroy(&fBagain));
    AssertPETSc(VecDestroy(&xEps));
    AssertPETSc(VecDestroy(&fEps));
    AssertPETSc(VecDestroy(&jv));
    AssertPETSc(VecDestroy(&jvReference));
#endif

    // This enters the production residual, where the total-Lagrangian model
    // calls enforceTractionBoundaries() for the runtime-selected patch.
    if (!solid->evolve())
    {
        FatalErrorInFunction
            << "The one-cell total-Lagrangian smoke solve did not converge"
            << abort(FatalError);
    }

    Info<< "PASS: Aróstica boundary runtime smoke" << nl
        << "    total-Lagrangian traction enforcement path exercised" << nl
        << "    formula/reference-area errors are checked by the companion "
        << "unit test" << endl;

    return 0;
}

// ************************************************************************* //
