/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Application
    ArosticaGlobalJVersusPDiagnosis

Description
    Action-only, block-aware MFFD/finite-difference audit for a mixed
    total-Lagrangian solid model.  The diagnostic uses only the public PETSc
    residual, Jacobian, and solution interfaces.  It never assembles a full
    finite-difference matrix and never advances accepted material state.

\*---------------------------------------------------------------------------*/

#ifdef USE_PETSC
    #include "petscErrorHandling.H"
    #include <petscmat.h>
    #include <petscksp.h>
#endif

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "nonLinGeomTotalLagTotalDispSolid.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <execinfo.h>
#include <iomanip>
#include <limits>
#include <signal.h>
#include <string>
#include <unistd.h>
#include <vector>

using namespace Foam;

#ifdef USE_PETSC
namespace
{

typedef solidModels::nonLinGeomTotalLagTotalDispSolid Model;

const scalar Ceps = 100.0;
const scalar Cref = 100.0;
const scalar Cfloor = 1.0e-12;
const scalar RelativePlateauTolerance = 1.0e-3;
const scalar CosineAcceptance = 0.999999;
const scalar RelativeAcceptance = 1.0e-3;


struct BlockNorms
{
    scalar full;
    scalar displacement;
    scalar pressure;
};


struct MffdContext
{
    Model* model;
    label mode;
    const volVectorField* baselineD;
    const pointVectorField* baselinePointD;
    const struct DiagnosticStateSnapshot* baselineState;
};


struct DiagnosticStateSnapshot
{
    tmp<surfaceTensorField> Ff;
    tmp<volTensorField> F;
    tmp<surfaceTensorField> relFf;
    tmp<volTensorField> relF;
    tmp<volScalarField> p;
    tmp<volSymmTensorField> sigma;
    tmp<volSymmTensorField> E;
    tmp<volSymmTensorField> Edot;
    tmp<volSymmTensorField> Sviscous;
    tmp<volSymmTensorField> sigmaViscous;
    tmp<volSymmTensorField> EOld;
    tmp<volSymmTensorField> EOldOld;
    tmp<surfaceSymmTensorField> Ef;
    tmp<surfaceSymmTensorField> Edotf;
    tmp<surfaceSymmTensorField> Sviscousf;
    tmp<surfaceSymmTensorField> sigmaViscousf;
    tmp<surfaceSymmTensorField> EfOld;
    tmp<surfaceSymmTensorField> EfOldOld;
};


struct MffdResult
{
    BlockNorms action;
    scalar h;
    scalar repeat;
};


struct ActionRecord
{
    label mode;
    scalar epsilon;
    BlockNorms fd;
    BlockNorms mffd;
    BlockNorms preconditioner;
    BlockNorms residualScale;
    scalar mffdRepeat;
    scalar relativeError[3];
    scalar absoluteError[3];
    scalar cosine[3];
    scalar fdAdjacentChange[3];
    label discrepancyCell[3];
};


typedef std::uint64_t StateHash;


struct StateHashes
{
    StateHash primary;
    StateHash deformation;
    StateHash history;
    StateHash boundary;
    StateHash time;
    StateHash solution;
    StateHash residual;
};


const StateHash FnvOffset = UINT64_C(1469598103934665603);
const StateHash FnvPrime = UINT64_C(1099511628211);


void addHashBytes
(
    StateHash& hash,
    const void* data,
    const std::size_t size
)
{
    const unsigned char* bytes =
        static_cast<const unsigned char*>(data);
    for (std::size_t byteI = 0; byteI < size; ++byteI)
    {
        hash ^= StateHash(bytes[byteI]);
        hash *= FnvPrime;
    }
}


template<class Type>
void addHashValue(StateHash& hash, const Type& value)
{
    addHashBytes(hash, &value, sizeof(Type));
}


template<class FieldType>
void addPrimitiveFieldHash(StateHash& hash, const FieldType& field)
{
    const label size = field.size();
    addHashValue(hash, size);
    forAll(field, valueI)
    {
        addHashValue(hash, field[valueI]);
    }
}


template<class Type>
void addPrimitiveFieldHash
(
    StateHash& hash,
    const pointPatchField<Type>& field
)
{
    const tmp<Field<Type>> values = field.patchInternalField();
    addPrimitiveFieldHash(hash, values());
}


template<class GeoField>
void addGeometricFieldHash(StateHash& hash, const GeoField& field)
{
    addPrimitiveFieldHash(hash, field.primitiveField());
    forAll(field.boundaryField(), patchI)
    {
        addPrimitiveFieldHash(hash, field.boundaryField()[patchI]);
    }
    const label fieldTimeIndex = field.timeIndex();
    addHashValue(hash, fieldTimeIndex);
}


StateHash globalHash(const StateHash localHash)
{
    unsigned long long local = localHash;
    unsigned long long global = 0;
    MPI_Allreduce
    (
        &local,
        &global,
        1,
        MPI_UNSIGNED_LONG_LONG,
        MPI_BXOR,
        PETSC_COMM_WORLD
    );
    return StateHash(global);
}


StateHash vectorHash(const Vec vector)
{
    StateHash hash = FnvOffset;
    if (!vector)
    {
        return globalHash(hash);
    }

    PetscInt localBegin = 0;
    PetscInt localEnd = 0;
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetOwnershipRange(vector, &localBegin, &localEnd));
    AssertPETSc(VecGetArrayRead(vector, &values));
    addHashValue(hash, localBegin);
    addHashValue(hash, localEnd);
    for (PetscInt valueI = 0; valueI < localEnd - localBegin; ++valueI)
    {
        const PetscReal value = PetscRealPart(values[valueI]);
        addHashValue(hash, value);
    }
    AssertPETSc(VecRestoreArrayRead(vector, &values));
    return globalHash(hash);
}


StateHashes stateHashes
(
    const Model& model,
    const Vec solution,
    const Vec residual
)
{
    StateHashes hashes;
    hashes.primary = FnvOffset;
    hashes.deformation = FnvOffset;
    hashes.history = FnvOffset;
    hashes.boundary = FnvOffset;
    hashes.time = FnvOffset;

    const volVectorField& D = model.D();
    const pointVectorField& pointD = model.pointD();
    const volScalarField& p =
        model.mesh().lookupObject<volScalarField>("p");
    addGeometricFieldHash(hashes.primary, D);
    addGeometricFieldHash(hashes.primary, D.oldTime());
    addGeometricFieldHash(hashes.primary, D.oldTime().oldTime());
    addGeometricFieldHash(hashes.primary, pointD);
    addGeometricFieldHash(hashes.primary, pointD.oldTime());
    addGeometricFieldHash(hashes.primary, pointD.oldTime().oldTime());
    addGeometricFieldHash(hashes.primary, p);
    addGeometricFieldHash(hashes.primary, p.oldTime());
    addGeometricFieldHash(hashes.primary, model.sigma());
    addGeometricFieldHash(hashes.primary, model.sigma().oldTime());

    const surfaceTensorField* Ff =
        model.mesh().findObject<surfaceTensorField>("Ff_");
    const volTensorField* F =
        model.mesh().findObject<volTensorField>("F_");
    const surfaceTensorField* relFf =
        model.mesh().findObject<surfaceTensorField>("relFf_");
    const volTensorField* relF =
        model.mesh().findObject<volTensorField>("relF_");
    const bool haveFf = Ff;
    const bool haveF = F;
    const bool haveRelFf = relFf;
    const bool haveRelF = relF;
    addHashValue(hashes.deformation, haveFf);
    addHashValue(hashes.deformation, haveF);
    addHashValue(hashes.deformation, haveRelFf);
    addHashValue(hashes.deformation, haveRelF);
    if (Ff)
    {
        addGeometricFieldHash(hashes.deformation, *Ff);
        addGeometricFieldHash(hashes.deformation, Ff->oldTime());
    }
    if (F)
    {
        addGeometricFieldHash(hashes.deformation, *F);
    }
    if (relFf)
    {
        addGeometricFieldHash(hashes.deformation, *relFf);
    }
    if (relF)
    {
        addGeometricFieldHash(hashes.deformation, *relF);
    }

    static const char* const cellHistoryNames[] =
    {
        "ArosticaE_myocardium",
        "ArosticaEdot_myocardium",
        "ArosticaSviscous_myocardium",
        "ArosticasigmaViscous_myocardium",
        "ArosticaEOld_myocardium",
        "ArosticaEOldOld_myocardium"
    };
    for
    (
        std::size_t fieldI = 0;
        fieldI < sizeof(cellHistoryNames)/sizeof(cellHistoryNames[0]);
        ++fieldI
    )
    {
        addGeometricFieldHash
        (
            hashes.history,
            model.mesh().lookupObject<volSymmTensorField>
            (cellHistoryNames[fieldI])
        );
    }

    static const char* const faceHistoryNames[] =
    {
        "ArosticaEf_myocardium",
        "ArosticaEdotf_myocardium",
        "ArosticaSviscousf_myocardium",
        "ArosticasigmaViscousf_myocardium",
        "ArosticaEfOld_myocardium",
        "ArosticaEfOldOld_myocardium"
    };
    for
    (
        std::size_t fieldI = 0;
        fieldI < sizeof(faceHistoryNames)/sizeof(faceHistoryNames[0]);
        ++fieldI
    )
    {
        addGeometricFieldHash
        (
            hashes.history,
            model.mesh().lookupObject<surfaceSymmTensorField>
            (faceHistoryNames[fieldI])
        );
    }

    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& patch = D.boundaryField()[patchI];
        addPrimitiveFieldHash(hashes.boundary, patch);
        const bool updated = patch.updated();
        addHashValue(hashes.boundary, updated);
        if (isA<fixedGradientFvPatchVectorField>(patch))
        {
            addPrimitiveFieldHash
            (
                hashes.boundary,
                refCast<const fixedGradientFvPatchVectorField>
                (patch).gradient()
            );
        }
    }

    const label timeIndex = model.mesh().time().timeIndex();
    const scalar timeValue = model.mesh().time().value();
    const scalar deltaT = model.mesh().time().deltaTValue();
    addHashValue(hashes.time, timeIndex);
    addHashValue(hashes.time, timeValue);
    addHashValue(hashes.time, deltaT);
    const label dOldTimeIndex = D.oldTime().timeIndex();
    const label dOldOldTimeIndex = D.oldTime().oldTime().timeIndex();
    const label pointDOldTimeIndex = pointD.oldTime().timeIndex();
    const label pointDOldOldTimeIndex =
        pointD.oldTime().oldTime().timeIndex();
    addHashValue(hashes.time, dOldTimeIndex);
    addHashValue(hashes.time, dOldOldTimeIndex);
    addHashValue(hashes.time, pointDOldTimeIndex);
    addHashValue(hashes.time, pointDOldOldTimeIndex);

    hashes.primary = globalHash(hashes.primary);
    hashes.deformation = globalHash(hashes.deformation);
    hashes.history = globalHash(hashes.history);
    hashes.boundary = globalHash(hashes.boundary);
    hashes.time = globalHash(hashes.time);
    hashes.solution = vectorHash(solution);
    hashes.residual = vectorHash(residual);
    return hashes;
}


Foam::string stateHashString(const StateHash hash)
{
    return Foam::string
    (
        std::to_string(static_cast<unsigned long long>(hash))
    );
}


void reportStateHashes
(
    const Model& model,
    const Vec solution,
    const Vec residual,
    const label mode,
    const word& stage
)
{
    const StateHashes hashes =
        stateHashes(model, solution, residual);
    Info<< "STATE_HASH mode=" << mode
        << " stage=" << stage
        << " primary=" << stateHashString(hashes.primary)
        << " deformation=" << stateHashString(hashes.deformation)
        << " history=" << stateHashString(hashes.history)
        << " boundary=" << stateHashString(hashes.boundary)
        << " solution=" << stateHashString(hashes.solution)
        << " residual=" << stateHashString(hashes.residual)
        << " time=" << stateHashString(hashes.time)
        << " timeIndex=" << model.mesh().time().timeIndex()
        << endl;
}


word modeName(const label mode);


void restoreDiagnosticBoundaryGradients
(
    Model& model,
    const volVectorField* baselineD
)
{
    if (!baselineD)
    {
        return;
    }

    forAll(model.D().boundaryFieldRef(), patchI)
    {
        const fvPatchField<vector>& baselinePatch =
            baselineD->boundaryField()[patchI];
        fvPatchField<vector>& trialPatch =
            model.D().boundaryFieldRef()[patchI];
        if
        (
            isA<fixedGradientFvPatchVectorField>(baselinePatch)
         && isA<fixedGradientFvPatchVectorField>(trialPatch)
        )
        {
            refCast<fixedGradientFvPatchVectorField>(trialPatch).gradient() =
                refCast<const fixedGradientFvPatchVectorField>
                (baselinePatch).gradient();
            // Keep this explicitly restored value for the next BC evaluate.
            trialPatch.setUpdated(true);
        }
    }
}


void diagnosticSignalHandler(const int signalNumber)
{
    void* frames[64];
    const int frameCount = backtrace(frames, 64);
    const char header[] = "DIAGNOSTIC_SIGNAL_BACKTRACE\n";
    write(STDERR_FILENO, header, sizeof(header) - 1);
    backtrace_symbols_fd(frames, frameCount, STDERR_FILENO);
    _exit(128 + signalNumber);
}


void prepareStateForResidual
(
    Model& model,
    const Vec x,
    const label mode,
    const word& stage,
    const volVectorField* baselineD,
    const pointVectorField* baselinePointD,
    const DiagnosticStateSnapshot* baselineState
)
{
    reportStateHashes
    (
        model,
        x,
        nullptr,
        mode,
        word(stage + "_entry")
    );
    if (baselineD && baselinePointD)
    {
        model.D() = *baselineD;
        model.D().oldTime() = baselineD->oldTime();
        model.D().oldTime().oldTime() = baselineD->oldTime().oldTime();
        model.pointD() = *baselinePointD;
        model.pointD().oldTime() = baselinePointD->oldTime();
        model.pointD().oldTime().oldTime() =
            baselinePointD->oldTime().oldTime();

        restoreDiagnosticBoundaryGradients(model, baselineD);
    }
    if (baselineState)
    {
        model.mesh().lookupObjectRef<volScalarField>("p") =
            baselineState->p();
        // solidTraction::updateCoeffs() uses the current Cauchy stress to
        // reconstruct D's fixed gradient. Restore the accepted stress before
        // any trial boundary correction so a preceding plus state cannot
        // contaminate the following base/minus callback.
        model.sigma() = baselineState->sigma();
    }
    // PETSc calls reuse one OpenFOAM field object for every trial.  Some
    // mixed/traction patches retain their updated() flag after reconstructing
    // a prior trial, which makes fvc::snGrad(D) reuse the prior fixed-gradient
    // state even after D and pointD have been restored.  Invalidate only the
    // diagnostic trial boundary fields; this does not alter production code
    // or the accepted parent case.
    forAll(model.D().boundaryFieldRef(), patchI)
    {
        model.D().boundaryFieldRef()[patchI].setUpdated(false);
    }
    Info<< "CHECKPOINT mode=" << mode << " " << stage
        << " displacement-unpack" << endl;
    vectorField& D = model.D();
    model.ExtractFieldComponents<vector>
    (
        x,
        D,
        0,
        model.twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );
    Info<< "CHECKPOINT mode=" << mode << " " << stage
        << " displacement-unpack-done" << endl;
    model.D().correctBoundaryConditions();
    restoreDiagnosticBoundaryGradients(model, baselineD);
    Info<< "CHECKPOINT mode=" << mode << " " << stage
        << " boundary-correction-done" << endl;
    model.mechanical().grad(model.D(), model.gradD());
    Info<< "CHECKPOINT mode=" << mode << " " << stage
        << " cell-gradient-done" << endl;
    model.D().correctBoundaryConditions();
    restoreDiagnosticBoundaryGradients(model, baselineD);
    model.mechanical().interpolate(model.D(), model.pointD());
    model.pointD().correctBoundaryConditions();

    if (baselineState)
    {
        surfaceTensorField& Ff =
            model.mesh().lookupObjectRef<surfaceTensorField>("Ff_");
        Ff = baselineState->Ff();
        Ff.oldTime() = baselineState->Ff().oldTime();
        model.mesh().lookupObjectRef<volTensorField>("F_") =
            baselineState->F();
        model.mesh().lookupObjectRef<surfaceTensorField>("relFf_") =
            baselineState->relFf();
        model.mesh().lookupObjectRef<volTensorField>("relF_") =
            baselineState->relF();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticaE_myocardium") = baselineState->E();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticaEdot_myocardium") = baselineState->Edot();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticaSviscous_myocardium") = baselineState->Sviscous();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticasigmaViscous_myocardium") =
            baselineState->sigmaViscous();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticaEOld_myocardium") = baselineState->EOld();
        model.mesh().lookupObjectRef<volSymmTensorField>
            ("ArosticaEOldOld_myocardium") = baselineState->EOldOld();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticaEfOld_myocardium") = baselineState->EfOld();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticaEfOldOld_myocardium") = baselineState->EfOldOld();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticaEf_myocardium") = baselineState->Ef();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticaEdotf_myocardium") = baselineState->Edotf();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticaSviscousf_myocardium") = baselineState->Sviscousf();
        model.mesh().lookupObjectRef<surfaceSymmTensorField>
            ("ArosticasigmaViscousf_myocardium") =
            baselineState->sigmaViscousf();
    }

    // formResidual's directConstitutive path obtains the face deformation
    // gradient through the mechanical-law update.  Refresh that same public
    // state before every trial so a preceding plus/minus evaluation cannot
    // leave a stale Ff field behind.
    surfaceTensorField gradDf
    (
        IOobject
        (
            "grad(D)f",
            model.mesh().time().timeName(),
            model.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        model.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        calculatedFvPatchTensorField::typeName
    );
    model.mechanical().grad(model.D(), model.pointD(), gradDf);
    scalar minFaceJ = VGREAT;
    scalar maxFaceJ = -VGREAT;
    label invalidFaceJ = 0;
    scalar maxPreparedD = 0.0;
    scalar maxPreparedPointD = 0.0;
    scalar maxPreparedBoundaryD = 0.0;
    scalar maxPreparedGradDf = 0.0;
    label maxPreparedGradDfFace = -1;
    label maxPreparedGradDfPatch = -1;
    label maxPreparedGradDfPatchFace = -1;
    forAll(model.D().primitiveField(), cellI)
    {
        maxPreparedD = max(maxPreparedD, mag(model.D()[cellI]));
    }
    forAll(model.pointD().primitiveField(), pointI)
    {
        maxPreparedPointD =
            max(maxPreparedPointD, mag(model.pointD()[pointI]));
    }
    forAll(model.D().boundaryField(), patchI)
    {
        forAll(model.D().boundaryField()[patchI], faceI)
        {
            maxPreparedBoundaryD = max
            (
                maxPreparedBoundaryD,
                mag(model.D().boundaryField()[patchI][faceI])
            );
        }
    }
    forAll(gradDf, faceI)
    {
        const scalar magnitude = mag(gradDf[faceI]);
        if (magnitude > maxPreparedGradDf)
        {
            maxPreparedGradDf = magnitude;
            maxPreparedGradDfFace = faceI;
        }
    }
    forAll(gradDf.boundaryField(), patchI)
    {
        forAll(gradDf.boundaryField()[patchI], faceI)
        {
            const scalar magnitude =
                mag(gradDf.boundaryField()[patchI][faceI]);
            if (magnitude > maxPreparedGradDf)
            {
                maxPreparedGradDf = magnitude;
                maxPreparedGradDfFace =
                    model.mesh().boundary()[patchI].start() + faceI;
                maxPreparedGradDfPatch = patchI;
                maxPreparedGradDfPatchFace = faceI;
            }
        }
    }
    reduce(maxPreparedD, maxOp<scalar>());
    reduce(maxPreparedPointD, maxOp<scalar>());
    reduce(maxPreparedBoundaryD, maxOp<scalar>());
    const surfaceScalarField faceJ
    (
        IOobject
        (
            "diagnosticFaceJ",
            model.mesh().time().timeName(),
            model.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(I + gradDf.T())
    );
    forAll(faceJ, faceI)
    {
        const scalar J = faceJ[faceI];
        minFaceJ = min(minFaceJ, J);
        maxFaceJ = max(maxFaceJ, J);
        invalidFaceJ += (!std::isfinite(J) || J <= 0);
    }
    forAll(faceJ.boundaryField(), patchI)
    {
        const scalarField& patchJ = faceJ.boundaryField()[patchI];
        forAll(patchJ, faceI)
        {
            minFaceJ = min(minFaceJ, patchJ[faceI]);
            maxFaceJ = max(maxFaceJ, patchJ[faceI]);
            invalidFaceJ +=
                (!std::isfinite(patchJ[faceI]) || patchJ[faceI] <= 0);
        }
    }
    reduce(minFaceJ, minOp<scalar>());
    reduce(maxFaceJ, maxOp<scalar>());
    reduce(invalidFaceJ, sumOp<label>());
    Info<< "PREPARED_FACE_STATE mode=" << mode << " stage=" << stage
        << " minJ=" << minFaceJ
        << " maxJ=" << maxFaceJ
        << " maxD=" << maxPreparedD
        << " maxBoundaryD=" << maxPreparedBoundaryD
        << " maxPointD=" << maxPreparedPointD
        << " maxGradDf=" << maxPreparedGradDf
        << " maxGradDfFace=" << maxPreparedGradDfFace
        << " maxGradDfPatch=" << maxPreparedGradDfPatch
        << " maxGradDfPatchFace=" << maxPreparedGradDfPatchFace
        << " invalidJ=" << invalidFaceJ << endl;
    if (maxPreparedGradDfFace >= 0)
    {
        const face& diagnosticFace =
            model.mesh().faces()[maxPreparedGradDfFace];
        Info<< "PREPARED_FACE_DETAIL mode=" << mode
            << " stage=" << stage
            << " face=" << maxPreparedGradDfFace
            << " owner=" << model.mesh().faceOwner()[maxPreparedGradDfFace]
            << " points=" << diagnosticFace << endl;
        forAll(diagnosticFace, pointI)
        {
            const label pointLabel = diagnosticFace[pointI];
            Info<< "  point=" << pointLabel
                << " value=" << model.pointD().primitiveField()[pointLabel]
                << endl;
        }
        if (maxPreparedGradDfFace < model.mesh().nInternalFaces())
        {
            Info<< "  neighbour="
                << model.mesh().faceNeighbour()[maxPreparedGradDfFace]
                << endl;
        }
        if (maxPreparedGradDfPatch >= 0)
        {
            const label patchI = maxPreparedGradDfPatch;
            const label patchFaceI = maxPreparedGradDfPatchFace;
            Info<< "  gradDfPatch="
                << gradDf.boundaryField()[patchI][patchFaceI]
                << " DPatch=" << model.D().boundaryField()[patchI][patchFaceI]
                << " ownerD="
                << model.D().primitiveField()
                   [model.mesh().faceOwner()[maxPreparedGradDfFace]]
                << endl;
        }
        else
        {
            Info<< "  gradDf=" << gradDf[maxPreparedGradDfFace] << endl;
        }
    }
    Info<< "CHECKPOINT mode=" << mode << " " << stage
        << " face-gradient-and-point-state-ready" << endl;
    reportStateHashes
    (
        model,
        x,
        nullptr,
        mode,
        word(stage + "_prepared")
    );
}


PetscErrorCode mffdResidual(void* context, Vec x, Vec f)
{
    MffdContext* ctx = static_cast<MffdContext*>(context);
    prepareStateForResidual
    (
        *ctx->model, x, ctx->mode, "mffd-residual-preparation"
        , ctx->baselineD, ctx->baselinePointD, ctx->baselineState
    );
    Info<< "CHECKPOINT mode=" << ctx->mode
        << " mffd-residual before-formResidual" << endl;
    const label status = ctx->model->formResidual(f, x);
    reportStateHashes
    (
        *ctx->model,
        x,
        f,
        ctx->mode,
        "mffd_residual_after"
    );
    Info<< "CHECKPOINT mode=" << ctx->mode
        << " mffd-residual after-formResidual status=" << status << endl;
    return status == 0
         ? PETSC_SUCCESS
         : PETSC_ERR_USER;
}


scalar vecNorm(Vec v)
{
    PetscReal value = 0.0;
    AssertPETSc(VecNorm(v, NORM_2, &value));
    return value;
}


scalar dotCosine(Vec a, Vec b)
{
    PetscScalar dot = 0.0;
    AssertPETSc(VecDot(a, b, &dot));
    return PetscRealPart(dot)/(vecNorm(a)*vecNorm(b) + VSMALL);
}


scalar blockDotCosine(Vec a, Vec b, const label component)
{
    PetscInt localSize = 0;
    PetscInt blockSize = 0;
    const PetscScalar* aValues = nullptr;
    const PetscScalar* bValues = nullptr;
    PetscReal localDot = 0.0;
    PetscReal localA = 0.0;
    PetscReal localB = 0.0;
    AssertPETSc(VecGetLocalSize(a, &localSize));
    AssertPETSc(VecGetBlockSize(a, &blockSize));
    AssertPETSc(VecGetArrayRead(a, &aValues));
    AssertPETSc(VecGetArrayRead(b, &bValues));
    for (PetscInt i = 0; i < localSize; ++i)
    {
        if
        (
            component == 0
         || (component == 1 && i % blockSize < 3)
         || (component == 2 && i % blockSize == 3)
        )
        {
            const PetscReal av = PetscRealPart(aValues[i]);
            const PetscReal bv = PetscRealPart(bValues[i]);
            localDot += av*bv;
            localA += av*av;
            localB += bv*bv;
        }
    }
    AssertPETSc(VecRestoreArrayRead(b, &bValues));
    AssertPETSc(VecRestoreArrayRead(a, &aValues));
    PetscReal dot = 0.0;
    PetscReal normA = 0.0;
    PetscReal normB = 0.0;
    MPI_Allreduce(&localDot, &dot, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&localA, &normA, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&localB, &normB, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    return dot/(Foam::sqrt(normA)*Foam::sqrt(normB) + VSMALL);
}


BlockNorms blockNorms(Vec v)
{
    PetscInt localSize = 0;
    PetscInt blockSize = 0;
    const PetscScalar* values = nullptr;
    PetscReal localD = 0.0;
    PetscReal localP = 0.0;
    PetscReal localFull = 0.0;

    AssertPETSc(VecGetLocalSize(v, &localSize));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    AssertPETSc(VecGetArrayRead(v, &values));

    for (PetscInt i = 0; i < localSize; ++i)
    {
        const PetscReal value = PetscRealPart(values[i]);
        localFull += value*value;
        if (blockSize > 0 && i % blockSize < 3)
        {
            localD += value*value;
        }
        else if (blockSize > 3 && i % blockSize == 3)
        {
            localP += value*value;
        }
    }
    AssertPETSc(VecRestoreArrayRead(v, &values));

    PetscReal full = 0.0;
    PetscReal displacement = 0.0;
    PetscReal pressure = 0.0;
    MPI_Allreduce(&localFull, &full, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce
    (
        &localD, &displacement, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localP, &pressure, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD
    );

    BlockNorms result;
    result.full = Foam::sqrt(full);
    result.displacement = Foam::sqrt(displacement);
    result.pressure = Foam::sqrt(pressure);
    return result;
}


label largestDiscrepancyCell(Vec a, Vec b, const label component)
{
    PetscInt localSize = 0;
    PetscInt blockSize = 0;
    const PetscScalar* aValues = nullptr;
    const PetscScalar* bValues = nullptr;
    scalar localMaximum = -1.0;
    label localCell = -1;

    AssertPETSc(VecGetLocalSize(a, &localSize));
    AssertPETSc(VecGetBlockSize(a, &blockSize));
    AssertPETSc(VecGetArrayRead(a, &aValues));
    AssertPETSc(VecGetArrayRead(b, &bValues));
    for (PetscInt i = 0; i < localSize; i += max(blockSize, PetscInt(1)))
    {
        scalar sum = 0.0;
        for (PetscInt componentI = 0; componentI < blockSize; ++componentI)
        {
            if
            (
                component == 0
             || (component == 1 && componentI < 3)
             || (component == 2 && componentI == 3)
            )
            {
                sum += sqr
                (
                    PetscRealPart
                    (
                        aValues[i + componentI] - bValues[i + componentI]
                    )
                );
            }
        }
        const scalar magnitude = Foam::sqrt(sum);
        if (magnitude > localMaximum)
        {
            localMaximum = magnitude;
            localCell = i/blockSize;
        }
    }
    AssertPETSc(VecRestoreArrayRead(b, &bValues));
    AssertPETSc(VecRestoreArrayRead(a, &aValues));

    // The reported cell is deliberately local.  A parallel diagnostic must
    // not invent a global cell map from an inaccessible distributed vector.
    return localCell;
}


void setBlockMode(Vec v, const fvMesh& mesh, const label mode)
{
    PetscInt localSize = 0;
    PetscInt blockSize = 0;
    PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(v, &localSize));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    AssertPETSc(VecGetArray(v, &values));

    const volVectorField* f0Ptr = mesh.findObject<volVectorField>("f0");
    const volVectorField* s0Ptr = mesh.findObject<volVectorField>("s0");

    Info<< "CHECKPOINT mode=" << mode << " fibre-lookup"
        << " f0=" << (f0Ptr ? "present" : "absent")
        << " s0=" << (s0Ptr ? "present" : "absent") << endl;

    if (blockSize < 4)
    {
        FatalErrorInFunction
            << "Expected a four-component PETSc cell block, got blockSize="
            << blockSize << abort(FatalError);
    }

    if (mode == 6 && !f0Ptr)
    {
        FatalErrorInFunction
            << "fibreDirection requires the registered volVectorField f0"
            << abort(FatalError);
    }

    Info<< "CHECKPOINT mode=" << mode << " fibre-field-class"
        << " class=" << (f0Ptr ? f0Ptr->type() : word("absent"))
        << " cells=" << mesh.nCells()
        << " internal=" << (f0Ptr ? f0Ptr->primitiveField().size() : 0)
        << " boundaryPatches=" << mesh.boundary().size() << endl;

    for (PetscInt localI = 0; localI < localSize; localI += blockSize)
    {
        const label cellI = localI/blockSize;
        const vector c = mesh.C()[cellI];
        vector d(0, 0, 0);
        scalar pressure = 0.0;

        switch (mode)
        {
            case 0: d = vector(1, -0.7, 0.4); break;
            case 1: pressure = 1.0; break;
            case 2: d = vector(1, 0.3, -0.5); pressure = 0.8; break;
            case 3: d = vector(1, 0, 0); break;
            case 4: d = vector(-c.y(), c.x(), 0); break;
            case 5: d = vector(c.x(), c.y(), c.z()); break;
            case 6:
                if (f0Ptr) d = f0Ptr->primitiveField()[cellI];
                break;
            case 7:
                if (s0Ptr) d = s0Ptr->primitiveField()[cellI];
                break;
            case 8: pressure = (cellI % 2) ? -1.0 : 1.0; break;
            case 9:
                pressure = scalar((37*cellI + 17*(cellI % 11) + 5) % 101 - 50)
                          /50.0;
                break;
            default: break;
        }

        values[localI] = d.x();
        values[localI + 1] = d.y();
        values[localI + 2] = d.z();
        if (blockSize > 3)
        {
            values[localI + 3] = pressure;
        }
    }
    AssertPETSc(VecRestoreArray(v, &values));
}


void reportModeVector
(
    Vec v,
    const fvMesh& mesh,
    const label mode,
    const scalar normalisation
)
{
    PetscInt localSize = 0;
    PetscInt globalSize = 0;
    PetscInt blockSize = 0;
    PetscInt ownershipBegin = 0;
    PetscInt ownershipEnd = 0;
    PetscInt localDisplacementEntries = 0;
    PetscInt localPressureEntries = 0;
    PetscInt localZeroEntries = 0;
    PetscInt localNonFiniteEntries = 0;
    PetscReal localL1 = 0.0;
    PetscReal localL2Sqr = 0.0;
    PetscReal localLinf = 0.0;
    PetscReal localMinimum = VGREAT;
    PetscReal localMaximum = -VGREAT;
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(v, &localSize));
    AssertPETSc(VecGetSize(v, &globalSize));
    AssertPETSc(VecGetBlockSize(v, &blockSize));
    AssertPETSc(VecGetOwnershipRange(v, &ownershipBegin, &ownershipEnd));
    AssertPETSc(VecGetArrayRead(v, &values));
    for (PetscInt i = 0; i < localSize; ++i)
    {
        const PetscReal value = PetscRealPart(values[i]);
        const bool finite = std::isfinite(value);
        localL1 += mag(value);
        localL2Sqr += value*value;
        localLinf = max(localLinf, mag(value));
        localMinimum = min(localMinimum, value);
        localMaximum = max(localMaximum, value);
        localZeroEntries += (value == 0.0);
        localNonFiniteEntries += !finite;
        if (blockSize > 0 && i % blockSize < 3)
        {
            ++localDisplacementEntries;
        }
        else if (blockSize > 3 && i % blockSize == 3)
        {
            ++localPressureEntries;
        }
    }
    AssertPETSc(VecRestoreArrayRead(v, &values));

    PetscReal L1 = 0.0;
    PetscReal L2Sqr = 0.0;
    PetscReal Linf = 0.0;
    PetscReal minimum = 0.0;
    PetscReal maximum = 0.0;
    PetscInt displacementEntries = 0;
    PetscInt pressureEntries = 0;
    PetscInt zeroEntries = 0;
    PetscInt nonFiniteEntries = 0;
    MPI_Allreduce(&localL1, &L1, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&localL2Sqr, &L2Sqr, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(&localLinf, &Linf, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce
    (
        &localMinimum, &minimum, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localMaximum, &maximum, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localDisplacementEntries, &displacementEntries, 1, MPIU_INT,
        MPI_SUM, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localPressureEntries, &pressureEntries, 1, MPIU_INT,
        MPI_SUM, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localZeroEntries, &zeroEntries, 1, MPIU_INT,
        MPI_SUM, PETSC_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localNonFiniteEntries, &nonFiniteEntries, 1, MPIU_INT,
        MPI_SUM, PETSC_COMM_WORLD
    );

    label boundaryCells = 0;
    if (mode == 6)
    {
        boolList boundaryCell(mesh.nCells(), false);
        forAll(mesh.boundary(), patchI)
        {
            const labelUList& cells = mesh.boundary()[patchI].faceCells();
            forAll(cells, faceI)
            {
                boundaryCell[cells[faceI]] = true;
            }
        }
        forAll(boundaryCell, cellI)
        {
            boundaryCells += boundaryCell[cellI];
        }
        reduce(boundaryCells, sumOp<label>());
    }

    Info<< "MODE_VECTOR mode=" << mode << " (" << modeName(mode) << ')'
        << " localSize=" << localSize
        << " globalSize=" << globalSize
        << " blockSize=" << blockSize
        << " displacementEntries=" << displacementEntries
        << " pressureEntries=" << pressureEntries
        << " ownership=" << ownershipBegin << ':' << ownershipEnd
        << " L1=" << L1
        << " L2=" << Foam::sqrt(L2Sqr)
        << " Linf=" << Linf
        << " min=" << minimum
        << " max=" << maximum
        << " zeroEntries=" << zeroEntries
        << " nonFiniteEntries=" << nonFiniteEntries
        << " normalisation=" << normalisation
        << " boundaryAdjacentCells=" << boundaryCells << endl;

    if (mode == 6)
    {
        const volVectorField* f0Ptr =
            mesh.findObject<volVectorField>("f0");
        scalar minimumMagnitude = VGREAT;
        scalar maximumMagnitude = -VGREAT;
        label minimumCell = -1;
        label maximumCell = -1;
        label invalidMagnitude = 0;
        if (f0Ptr)
        {
            forAll(f0Ptr->primitiveField(), cellI)
            {
                const scalar magnitude = mag
                (
                    f0Ptr->primitiveField()[cellI]
                );
                if (magnitude < minimumMagnitude)
                {
                    minimumMagnitude = magnitude;
                    minimumCell = cellI;
                }
                if (magnitude > maximumMagnitude)
                {
                    maximumMagnitude = magnitude;
                    maximumCell = cellI;
                }
                invalidMagnitude +=
                    (!std::isfinite(magnitude) || magnitude <= VSMALL);
            }
        }
        reduce(minimumMagnitude, minOp<scalar>());
        reduce(maximumMagnitude, maxOp<scalar>());
        reduce(invalidMagnitude, sumOp<label>());
        Info<< "FIBRE_STATS cells=" << mesh.nCells()
            << " minMag=" << minimumMagnitude
            << " minCellLocal=" << minimumCell
            << " maxMag=" << maximumMagnitude
            << " maxCellLocal=" << maximumCell
            << " invalidMagnitudeCells=" << invalidMagnitude << endl;
    }
}


void reportPerturbedState
(
    const Model& model,
    Vec base,
    Vec candidate,
    const label mode,
    const scalar epsilon,
    const word& stage
)
{
    Vec difference = nullptr;
    AssertPETSc(VecDuplicate(base, &difference));
    AssertPETSc(VecCopy(candidate, difference));
    AssertPETSc(VecAXPY(difference, -1.0, base));
    BlockNorms delta = blockNorms(difference);
    PetscReal localMaximumDeltaD = 0.0;
    PetscInt localSize = 0;
    PetscInt blockSize = 0;
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetLocalSize(difference, &localSize));
    AssertPETSc(VecGetBlockSize(difference, &blockSize));
    AssertPETSc(VecGetArrayRead(difference, &values));
    for (PetscInt i = 0; i < localSize; ++i)
    {
        if (blockSize > 0 && i % blockSize < 3)
        {
            localMaximumDeltaD = max
            (
                localMaximumDeltaD,
                mag(PetscRealPart(values[i]))
            );
        }
    }
    AssertPETSc(VecRestoreArrayRead(difference, &values));
    PetscReal maximumDeltaD = 0.0;
    MPI_Allreduce
    (
        &localMaximumDeltaD, &maximumDeltaD, 1, MPIU_REAL,
        MPI_MAX, PETSC_COMM_WORLD
    );
    AssertPETSc(VecDestroy(&difference));

    const volVectorField& D = model.D();
    const pointVectorField& pointD = model.pointD();
    scalar maxD = 0.0;
    scalar maxPointD = 0.0;
    forAll(D.primitiveField(), cellI)
    {
        maxD = max(maxD, mag(D.primitiveField()[cellI]));
    }
    forAll(pointD.primitiveField(), pointI)
    {
        maxPointD = max(maxPointD, mag(pointD.primitiveField()[pointI]));
    }
    reduce(maxD, maxOp<scalar>());
    reduce(maxPointD, maxOp<scalar>());

    const volScalarField* JPtr = model.mesh().findObject<volScalarField>("J");
    scalar minJ = VGREAT;
    scalar maxJ = -VGREAT;
    label invalidJ = 0;
    if (JPtr)
    {
        forAll(JPtr->primitiveField(), cellI)
        {
            const scalar J = JPtr->primitiveField()[cellI];
            minJ = min(minJ, J);
            maxJ = max(maxJ, J);
            invalidJ += (!std::isfinite(J) || J <= 0);
        }
    }
    reduce(minJ, minOp<scalar>());
    reduce(maxJ, maxOp<scalar>());
    reduce(invalidJ, sumOp<label>());

    Info<< "PERTURBED_STATE mode=" << mode
        << " epsilon=" << epsilon
        << " stage=" << stage
        << " deltaD_l2=" << delta.displacement
        << " maxDeltaD=" << maximumDeltaD
        << " maxD=" << maxD
        << " maxPointD=" << maxPointD
        << " minJ=" << minJ
        << " maxJ=" << maxJ
        << " invalidJ=" << invalidJ
        << " finite=" << (invalidJ == 0 ? "yes" : "no") << endl;
}


word modeName(const label mode)
{
    static const word names[] =
    {
        "uniformDisplacement",
        "uniformPressure",
        "coupledDisplacementPressure",
        "translation",
        "rotation",
        "dilation",
        "fibreDirection",
        "sheetDirection",
        "checkerboardPressure",
        "deterministicPressure"
    };
    return mode >= 0 && mode < 10 ? names[mode] : "unknown";
}


BlockNorms maximumResidualScale
(
    const BlockNorms& base,
    const BlockNorms& plus,
    const BlockNorms& minus
)
{
    BlockNorms result;
    result.full = max(base.full, max(plus.full, minus.full));
    result.displacement = max
    (
        base.displacement,
        max(plus.displacement, minus.displacement)
    );
    result.pressure = max(base.pressure, max(plus.pressure, minus.pressure));
    return result;
}


void centralFDActionWithBase
(
    Model& model,
    Vec base,
    const BlockNorms& baseResidualScale,
    Vec direction,
    const label mode,
    const scalar epsilon,
    Vec action,
    BlockNorms& residualScale,
    const volVectorField* baselineD,
    const pointVectorField* baselinePointD
    , const DiagnosticStateSnapshot* baselineState
)
{
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " allocate-vectors" << endl;
    Vec plus = nullptr;
    Vec minus = nullptr;
    Vec residualPlus = nullptr;
    Vec residualMinus = nullptr;
    AssertPETSc(VecDuplicate(base, &plus));
    AssertPETSc(VecDuplicate(base, &minus));
    AssertPETSc(VecDuplicate(base, &residualPlus));
    AssertPETSc(VecDuplicate(base, &residualMinus));
    AssertPETSc(VecCopy(base, plus));
    AssertPETSc(VecAXPY(plus, epsilon, direction));
    AssertPETSc(VecCopy(base, minus));
    AssertPETSc(VecAXPY(minus, -epsilon, direction));
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " constructed-plus-minus" << endl;
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " restore-before-plus" << endl;
    prepareStateForResidual
    (
        model, base, mode, "central-restore-before-plus"
        , baselineD, baselinePointD, baselineState
    );
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " restore-after-plus" << endl;
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " plus-before-formResidual" << endl;
    prepareStateForResidual
    (
        model, plus, mode, "central-plus"
        , baselineD, baselinePointD, baselineState
    );
    AssertPETSc(model.formResidual(residualPlus, plus));
    reportStateHashes
    (
        model,
        plus,
        residualPlus,
        mode,
        "central_plus_after"
    );
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " plus-after-formResidual" << endl;
    reportPerturbedState(model, base, plus, mode, epsilon, "plus-after");
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " restore-before-minus" << endl;
    prepareStateForResidual
    (
        model, base, mode, "central-restore-before-minus"
        , baselineD, baselinePointD, baselineState
    );
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " restore-after-minus" << endl;
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " minus-before-formResidual" << endl;
    prepareStateForResidual
    (
        model, minus, mode, "central-minus"
        , baselineD, baselinePointD, baselineState
    );
    AssertPETSc(model.formResidual(residualMinus, minus));
    reportStateHashes
    (
        model,
        minus,
        residualMinus,
        mode,
        "central_minus_after"
    );
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " minus-after-formResidual" << endl;
    reportPerturbedState(model, base, minus, mode, epsilon, "minus-after");
    Info<< "CHECKPOINT central-FD mode=" << mode << " epsilon=" << epsilon
        << " restore-after-minus" << endl;
    prepareStateForResidual
    (
        model, base, mode, "central-restore-after-minus"
        , baselineD, baselinePointD, baselineState
    );

    residualScale = maximumResidualScale
    (
        baseResidualScale,
        blockNorms(residualPlus),
        blockNorms(residualMinus)
    );
    AssertPETSc(VecCopy(residualPlus, action));
    AssertPETSc(VecAXPY(action, -1.0, residualMinus));
    AssertPETSc(VecScale(action, 1.0/(2.0*epsilon)));

    AssertPETSc(VecDestroy(&plus));
    AssertPETSc(VecDestroy(&minus));
    AssertPETSc(VecDestroy(&residualPlus));
    AssertPETSc(VecDestroy(&residualMinus));
}


MffdResult mffdAction
(
    Mat mffd,
    Vec direction,
    Vec first,
    Vec second
)
{
    Info<< "CHECKPOINT MFFD before-first-MatMult" << endl;
    AssertPETSc(MatMult(mffd, direction, first));
    Info<< "CHECKPOINT MFFD after-first-MatMult" << endl;
    Info<< "CHECKPOINT MFFD before-second-MatMult" << endl;
    AssertPETSc(MatMult(mffd, direction, second));
    Info<< "CHECKPOINT MFFD after-second-MatMult" << endl;
    scalar h = 0.0;
    AssertPETSc(MatMFFDGetH(mffd, &h));

    Vec difference = nullptr;
    AssertPETSc(VecDuplicate(first, &difference));
    AssertPETSc(VecCopy(first, difference));
    AssertPETSc(VecAXPY(difference, -1.0, second));

    MffdResult result;
    result.action = blockNorms(first);
    result.h = h;
    result.repeat = vecNorm(difference)/(vecNorm(first) + VSMALL);
    AssertPETSc(VecDestroy(&difference));
    return result;
}


scalar median(std::vector<scalar> values)
{
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size()/2;
    if (values.size() % 2) return values[middle];
    return 0.5*(values[middle - 1] + values[middle]);
}


scalar actionThreshold
(
    const scalar epsilon,
    const scalar residualScale,
    const scalar actionReference
)
{
    const scalar epsilonMachine = std::numeric_limits<scalar>::epsilon();
    return max
    (
        Ceps*epsilonMachine*residualScale/epsilon,
        max(Cref*epsilonMachine*actionReference, Cfloor*actionReference)
    );
}


scalar actionComponent(const BlockNorms& norms, const label component)
{
    return component == 0
         ? norms.full
         : component == 1 ? norms.displacement : norms.pressure;
}


word componentName(const label component)
{
    return component == 0 ? "complete" : component == 1 ? "momentum" : "pressure";
}


word baseClassification
(
    const scalar fd,
    const scalar mffd,
    const scalar tau
)
{
    const bool fdResolved = fd > tau;
    const bool mffdResolved = mffd > tau;
    if (!fdResolved && !mffdResolved) return "both-actions-numerically-zero";
    if (!fdResolved && mffdResolved) return "FD-zero/MFFD-nonzero";
    if (fdResolved && !mffdResolved) return "MFFD-zero/FD-nonzero";
    return "resolved";
}


void reportThreshold
(
    const ActionRecord& record,
    const label component,
    const scalar reference
)
{
    const scalar tau = actionThreshold
    (
        record.epsilon,
        actionComponent(record.residualScale, component),
        reference
    );
    Info<< "THRESHOLD mode=" << record.mode
        << " (" << modeName(record.mode) << ")"
        << " epsilon=" << record.epsilon
        << " action=" << componentName(component)
        << " Rscale=" << actionComponent(record.residualScale, component)
        << " Aref=" << reference
        << " tau=" << tau
        << " constants=(Ceps " << Ceps << ", Cref " << Cref
        << ", Cfloor " << Cfloor << ')' << nl;
}


bool suitablePair
(
    const ActionRecord& first,
    const ActionRecord& second,
    const label component,
    const scalar reference
)
{
    const scalar tauFirst = actionThreshold
    (
        first.epsilon,
        actionComponent(first.residualScale, component),
        reference
    );
    const scalar tauSecond = actionThreshold
    (
        second.epsilon,
        actionComponent(second.residualScale, component),
        reference
    );
    const scalar firstFD = actionComponent(first.fd, component);
    const scalar secondFD = actionComponent(second.fd, component);
    return
        firstFD > tauFirst
     && secondFD > tauSecond
     && actionComponent(first.mffd, component) > tauFirst
     && actionComponent(second.mffd, component) > tauSecond
     && first.fdAdjacentChange[component] <= RelativePlateauTolerance
     && second.fdAdjacentChange[component] <= RelativePlateauTolerance
     && first.discrepancyCell[component] == second.discrepancyCell[component];
}


void reportInterpretation
(
    const word& stateLabel,
    const bool pressureDefect,
    const bool anyDiscrepancy
)
{
    Info<< "INTERPRETATION ";
    if (stateLabel == "accepted-exact" || stateLabel == "accepted-complete")
    {
        if (pressureDefect)
        {
            Info<< "PRESSURE-BLOCK MFFD/FD DEFECT";
        }
        else if (anyDiscrepancy)
        {
            Info<< "MFFD/FD DISCREPANCY AFTER BLOCK-AWARE ACCEPTANCE";
        }
        else
        {
            Info<< "MFFD/FD ACTIONS RESOLVED AND CONSISTENT";
        }
    }
    else
    {
        Info<< "LIFECYCLE " << stateLabel
            << "; MFFD/FD RESULTS DIAGNOSTIC ONLY";
    }
    Info<< endl;
}


void reportTermDecomposition(Model& model)
{
    Model::momentumResidualDecompositionData decomposition;
    model.momentumResidualDecomposition(decomposition);
    Info<< "TERM_DECOMPOSITION current-state" << nl
        << "  terms=" << decomposition.names << endl;
    forAll(decomposition.names, termI)
    {
        scalar cellSqr = 0.0;
        scalar faceSqr = 0.0;
        forAll(decomposition.cellActions[termI], cellI)
        {
            cellSqr += magSqr(decomposition.cellActions[termI][cellI]);
        }
        forAll(decomposition.faceForces[termI], faceI)
        {
            faceSqr += magSqr(decomposition.faceForces[termI][faceI]);
        }
        reduce(cellSqr, sumOp<scalar>());
        reduce(faceSqr, sumOp<scalar>());
        Info<< "  TERM_DECOMPOSITION term=" << decomposition.names[termI]
            << " cellNorm=" << Foam::sqrt(cellSqr)
            << " faceNorm=" << Foam::sqrt(faceSqr) << endl;
    }
}

}
#endif


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "actionAudit",
        "Run the action-only block-aware MFFD/central-FD audit"
    );
    argList::addOption
    (
        "stateLabel",
        "word",
        "Lifecycle state: pre-step, accepted-incomplete, accepted-complete, or unavailable"
    );
    argList::addOption
    (
        "mode",
        "label",
        "Run only one block mode (0-9), for example -mode 6"
    );
    argList::addBoolOption
    (
        "skipMffd",
        "Skip the PETSc MFFD action (diagnostic isolation only)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

#ifndef USE_PETSC
    FatalErrorInFunction
        << "ArosticaGlobalJVersusPDiagnosis requires USE_PETSC"
        << abort(FatalError);
    return 1;
#else
    word stateLabel("unavailable");
    args.optionReadIfPresent("stateLabel", stateLabel);
    label modeFilter = -1;
    const bool skipMffd = args.optionFound("skipMffd");
    if (args.optionFound("mode"))
    {
        args.optionReadIfPresent("mode", modeFilter);
        if (modeFilter < 0 || modeFilter >= 10)
        {
            FatalErrorInFunction
                << "-mode must be between 0 and 9, got " << modeFilter
                << abort(FatalError);
        }
        Info<< "MODE FILTER " << modeFilter << " ("
            << modeName(modeFilter) << ')' << endl;
    }
    if (!args.optionFound("actionAudit"))
    {
        Info<< "No explicit action switch supplied; running -actionAudit" << nl;
    }

    autoPtr<solidModel> solid =
        solidModel::New(runTime, dynamicFvMesh::defaultRegion);
    Model& model = refCast<Model>(*solid);
    signal(SIGILL, diagnosticSignalHandler);

    Vec x = nullptr;
    Vec baseResidual = nullptr;
    Vec residualRepeat = nullptr;
    AssertPETSc(model.initialiseSolution(x));
    AssertPETSc(VecDuplicate(x, &baseResidual));
    AssertPETSc(VecDuplicate(x, &residualRepeat));
    reportStateHashes
    (
        model,
        x,
        nullptr,
        -1,
        "initialisation_residual_before"
    );
    AssertPETSc(model.formResidual(baseResidual, x));
    reportStateHashes
    (
        model,
        x,
        baseResidual,
        -1,
        "initialisation_residual_after"
    );
    reportStateHashes
    (
        model,
        x,
        nullptr,
        -1,
        "initialisation_repeat_before"
    );
    AssertPETSc(model.formResidual(residualRepeat, x));
    reportStateHashes
    (
        model,
        x,
        residualRepeat,
        -1,
        "initialisation_repeat_after"
    );

    tmp<volVectorField> baselineD = model.D().clone();
    tmp<pointVectorField> baselinePointD = model.pointD().clone();
    DiagnosticStateSnapshot baselineState;
    baselineState.Ff =
        model.mesh().lookupObject<surfaceTensorField>("Ff_").clone();
    baselineState.F =
        model.mesh().lookupObject<volTensorField>("F_").clone();
    baselineState.relFf =
        model.mesh().lookupObject<surfaceTensorField>("relFf_").clone();
    baselineState.relF =
        model.mesh().lookupObject<volTensorField>("relF_").clone();
    baselineState.p =
        model.mesh().lookupObject<volScalarField>("p").clone();
    baselineState.sigma = model.sigma().clone();
    baselineState.E =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticaE_myocardium").clone();
    baselineState.Edot =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticaEdot_myocardium").clone();
    baselineState.Sviscous =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticaSviscous_myocardium").clone();
    baselineState.sigmaViscous =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticasigmaViscous_myocardium").clone();
    baselineState.EOld =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticaEOld_myocardium").clone();
    baselineState.EOldOld =
        model.mesh().lookupObject<volSymmTensorField>
        ("ArosticaEOldOld_myocardium").clone();
    baselineState.EfOld =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticaEfOld_myocardium").clone();
    baselineState.EfOldOld =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticaEfOldOld_myocardium").clone();
    baselineState.Ef =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticaEf_myocardium").clone();
    baselineState.Edotf =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticaEdotf_myocardium").clone();
    baselineState.Sviscousf =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticaSviscousf_myocardium").clone();
    baselineState.sigmaViscousf =
        model.mesh().lookupObject<surfaceSymmTensorField>
        ("ArosticasigmaViscousf_myocardium").clone();
    reportStateHashes
    (
        model,
        x,
        baseResidual,
        -1,
        "accepted_baseline"
    );

    prepareStateForResidual
    (
        model,
        x,
        -1,
        "base-residual",
        &baselineD(),
        &baselinePointD(),
        &baselineState
    );
    reportStateHashes
    (
        model,
        x,
        nullptr,
        -1,
        "base_residual_before"
    );
    AssertPETSc(model.formResidual(baseResidual, x));
    reportStateHashes
    (
        model,
        x,
        baseResidual,
        -1,
        "base_residual_after"
    );
    prepareStateForResidual
    (
        model,
        x,
        -1,
        "repeated-base-residual",
        &baselineD(),
        &baselinePointD(),
        &baselineState
    );
    reportStateHashes
    (
        model,
        x,
        nullptr,
        -1,
        "repeat_residual_before"
    );
    AssertPETSc(model.formResidual(residualRepeat, x));
    reportStateHashes
    (
        model,
        x,
        residualRepeat,
        -1,
        "repeat_residual_after"
    );
    AssertPETSc(VecAXPY(residualRepeat, -1.0, baseResidual));
    reportStateHashes
    (
        model,
        x,
        residualRepeat,
        -1,
        "repeat_residual_difference"
    );

    const BlockNorms baseResidualScale = blockNorms(baseResidual);
    Info<< "MFFD/FD BLOCK-AWARE ACTION AUDIT" << nl
        << "  cells = " << model.mesh().nCells() << nl
        << "  residual repeat norm = " << vecNorm(residualRepeat) << nl
        << "  lifecycle state = " << stateLabel << nl
        << "  accepted-state primary comparison = "
        << ((stateLabel == "accepted-exact" || stateLabel == "accepted-complete")
            ? "allowed" : "not allowed") << nl
        << "  Krylov vector semantics: exact Arnoldi/Krylov basis vectors are not "
        << "exposed by this public API; retained vectors are standalone-"
        << "reconstructed search directions, plus unpreconditioned residuals "
        << "and solution updates where replay is requested." << nl
        << "  constants: Ceps=" << Ceps << " Cref=" << Cref
        << " Cfloor=" << Cfloor << endl;

    Mat P = nullptr;
    AssertPETSc(model.initialiseJacobian(P));
    AssertPETSc(model.formJacobian(P, x));
    AssertPETSc(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
    AssertPETSc(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
    reportTermDecomposition(model);

    MffdContext context;
    context.model = &model;
    context.mode = -1;
    context.baselineD = &baselineD();
    context.baselinePointD = &baselinePointD();
    context.baselineState = &baselineState;
    PetscInt n = 0;
    AssertPETSc(VecGetSize(x, &n));
    Mat mffd = nullptr;
    AssertPETSc
    (
        MatCreateMFFD(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, &mffd)
    );
    AssertPETSc(MatMFFDSetFunction(mffd, mffdResidual, &context));
    AssertPETSc(MatMFFDSetType(mffd, MATMFFD_DS));
    AssertPETSc(MatMFFDSetBase(mffd, x, baseResidual));

    scalarField epsilons(5);
    epsilons[0] = 1.0e-4;
    epsilons[1] = 1.0e-5;
    epsilons[2] = 1.0e-6;
    epsilons[3] = 1.0e-7;
    epsilons[4] = 1.0e-8;
    Info<< "  central-FD epsilons = " << epsilons << nl
        << "  acceptance pair = " << epsilons[2] << " and " << epsilons[3]
        << nl << "  plateau tolerance = " << RelativePlateauTolerance << nl
        << "  relative acceptance = " << RelativeAcceptance
        << ", cosine acceptance = " << CosineAcceptance << endl;

    std::vector<scalar> fullReferences;
    std::vector<scalar> displacementReferences;
    std::vector<scalar> pressureReferences;
    std::vector<ActionRecord> records;

    for (label mode = 0; mode < 10; ++mode)
    {
        if (modeFilter >= 0 && mode != modeFilter)
        {
            continue;
        }
        Vec direction = nullptr;
        Vec mffdVector = nullptr;
        Vec mffdRepeatVector = nullptr;
        Vec pVector = nullptr;
        Vec fdVector = nullptr;
        Vec previousFD = nullptr;
        AssertPETSc(VecDuplicate(x, &direction));
        AssertPETSc(VecDuplicate(x, &mffdVector));
        AssertPETSc(VecDuplicate(x, &mffdRepeatVector));
        AssertPETSc(VecDuplicate(x, &pVector));
        AssertPETSc(VecDuplicate(x, &fdVector));
        AssertPETSc(VecDuplicate(x, &previousFD));
        AssertPETSc(VecSet(direction, 0.0));
        setBlockMode(direction, model.mesh(), mode);
        const scalar directionNorm = vecNorm(direction);
        if (directionNorm > VSMALL)
        {
            AssertPETSc(VecScale(direction, 1.0/directionNorm));
        }
        reportModeVector
        (
            direction,
            model.mesh(),
            mode,
            directionNorm
        );

        AssertPETSc(MatMult(P, direction, pVector));
        context.mode = mode;
        const MffdResult mffdResult =
            skipMffd
          ? MffdResult{BlockNorms{0.0, 0.0, 0.0}, 0.0, 0.0}
          : mffdAction(mffd, direction, mffdVector, mffdRepeatVector);
        const BlockNorms mffdNorms = blockNorms(mffdVector);
        const BlockNorms pNorms = blockNorms(pVector);
        Info<< "MODE " << mode << " " << modeName(mode)
            << " vector=standalone-reconstructed-search-direction"
            << " MFFD_h=" << mffdResult.h
            << " MFFD_repeat=" << mffdResult.repeat << endl;

        forAll(epsilons, epsilonI)
        {
            ActionRecord record;
            record.mode = mode;
            record.epsilon = epsilons[epsilonI];
            centralFDActionWithBase
            (
                model,
                x,
                baseResidualScale,
                direction,
                mode,
                record.epsilon,
                fdVector,
                record.residualScale,
                &baselineD(),
                &baselinePointD(),
                &baselineState
            );
            record.fd = blockNorms(fdVector);
            record.mffd = mffdNorms;
            record.preconditioner = pNorms;
            record.mffdRepeat = mffdResult.repeat;

            if (epsilonI == 0)
            {
                record.fdAdjacentChange[0] = GREAT;
                record.fdAdjacentChange[1] = GREAT;
                record.fdAdjacentChange[2] = GREAT;
            }
            else
            {
                Vec difference = nullptr;
                AssertPETSc(VecDuplicate(fdVector, &difference));
                AssertPETSc(VecCopy(fdVector, difference));
                AssertPETSc(VecAXPY(difference, -1.0, previousFD));
                const scalar denominator = vecNorm(fdVector) + VSMALL;
                const scalar differenceNorm = vecNorm(difference);
                const BlockNorms differenceBlocks = blockNorms(difference);
                record.fdAdjacentChange[0] =
                    differenceNorm/denominator;
                record.fdAdjacentChange[1] =
                    differenceBlocks.displacement/(record.fd.displacement + VSMALL);
                record.fdAdjacentChange[2] =
                    differenceBlocks.pressure/(record.fd.pressure + VSMALL);
                AssertPETSc(VecDestroy(&difference));
            }
            AssertPETSc(VecCopy(fdVector, previousFD));

            Vec difference = nullptr;
            AssertPETSc(VecDuplicate(fdVector, &difference));
            AssertPETSc(VecCopy(mffdVector, difference));
            AssertPETSc(VecAXPY(difference, -1.0, fdVector));
            const BlockNorms differenceBlocks = blockNorms(difference);
            const scalar differenceNorm = vecNorm(difference);
            const scalar fdFull = record.fd.full + VSMALL;
            const scalar fdD = record.fd.displacement + VSMALL;
            const scalar fdP = record.fd.pressure + VSMALL;
            record.absoluteError[0] = differenceNorm;
            record.absoluteError[1] = differenceBlocks.displacement;
            record.absoluteError[2] = differenceBlocks.pressure;
            record.relativeError[0] = differenceNorm/fdFull;
            record.relativeError[1] = differenceBlocks.displacement/fdD;
            record.relativeError[2] = differenceBlocks.pressure/fdP;
            record.cosine[0] = dotCosine(fdVector, mffdVector);
            record.cosine[1] = blockDotCosine(fdVector, mffdVector, 1);
            record.cosine[2] = blockDotCosine(fdVector, mffdVector, 2);
            record.discrepancyCell[0] =
                largestDiscrepancyCell(fdVector, mffdVector, 0);
            record.discrepancyCell[1] =
                largestDiscrepancyCell(fdVector, mffdVector, 1);
            record.discrepancyCell[2] =
                largestDiscrepancyCell(fdVector, mffdVector, 2);
            AssertPETSc(VecDestroy(&difference));

            if (epsilonI == 2 || epsilonI == 3)
            {
                fullReferences.push_back(record.fd.full);
                if (mode != 1 && mode != 8 && mode != 9)
                {
                    displacementReferences.push_back(record.fd.displacement);
                }
                if (mode == 1 || mode == 8 || mode == 9)
                {
                    pressureReferences.push_back(record.fd.pressure);
                }
            }
            records.push_back(record);
        }

        AssertPETSc(VecDestroy(&direction));
        AssertPETSc(VecDestroy(&mffdVector));
        AssertPETSc(VecDestroy(&mffdRepeatVector));
        AssertPETSc(VecDestroy(&pVector));
        AssertPETSc(VecDestroy(&fdVector));
        AssertPETSc(VecDestroy(&previousFD));
    }

    if (modeFilter >= 0)
    {
        Info<< "PASS: standalone mode-only diagnostic completed for mode "
            << modeFilter << " (" << modeName(modeFilter) << ')' << endl;
        AssertPETSc(MatDestroy(&mffd));
        AssertPETSc(MatDestroy(&P));
        AssertPETSc(VecDestroy(&baseResidual));
        AssertPETSc(VecDestroy(&residualRepeat));
        AssertPETSc(VecDestroy(&x));
        return 0;
    }

    // The two middle epsilons are the documented plateau candidates.  Their
    // own block actions, not the complete mixed action, define Aref.
    scalar ArefFull = median(fullReferences);
    scalar ArefD = median(displacementReferences);
    scalar ArefP = median(pressureReferences);
    Info<< "ACTION REFERENCES (median plateau norms)" << nl
        << "  Aref_full = " << ArefFull << nl
        << "  Aref_D = " << ArefD << nl
        << "  Aref_p = " << ArefP << endl;

    bool pressureDefect = false;
    bool anyDiscrepancy = false;
    forAll(records, recordI)
    {
        ActionRecord& record = records[recordI];
        const scalar references[] = {ArefFull, ArefD, ArefP};
        for (label component = 0; component < 3; ++component)
        {
            reportThreshold(record, component, references[component]);
            const scalar tau = actionThreshold
            (
                record.epsilon,
                actionComponent(record.residualScale, component),
                references[component]
            );
            const word classification = baseClassification
            (
                actionComponent(record.fd, component),
                actionComponent(record.mffd, component),
                tau
            );
            const bool resolved =
                actionComponent(record.fd, component) > tau
             && actionComponent(record.mffd, component) > tau;
            const bool accepted = resolved
             && record.fdAdjacentChange[component] <= RelativePlateauTolerance
             && record.mffdRepeat <= 1.0e-12
             && record.relativeError[component] <= RelativeAcceptance
             && record.cosine[component] >= CosineAcceptance;
            if (resolved && !accepted) anyDiscrepancy = true;
            Info<< "RESULT mode=" << record.mode
                << " epsilon=" << record.epsilon
                << " action=" << componentName(component)
                << " FD=" << actionComponent(record.fd, component)
                << " MFFD=" << actionComponent(record.mffd, component)
                << " P=" << actionComponent(record.preconditioner, component)
                << " absError=" << record.absoluteError[component]
                << " relativeError=" << record.relativeError[component]
                << " cosine=" << record.cosine[component]
                << " FD_adjacent_change="
                << record.fdAdjacentChange[component]
                << " discrepancy_local_cell="
                << record.discrepancyCell[component]
                << " classification=" << classification
                << " acceptance=" << (accepted ? "pass" : "not-accepted")
                << endl;
        }

    }

    for (label mode = 0; mode < 10; ++mode)
    {
        if (modeFilter >= 0 && mode != modeFilter)
        {
            continue;
        }
        const ActionRecord& first = records[mode*epsilons.size() + 2];
        const ActionRecord& second = records[mode*epsilons.size() + 3];
        const scalar tauFirstP = actionThreshold
        (
            first.epsilon, first.residualScale.pressure, ArefP
        );
        const scalar tauSecondP = actionThreshold
        (
            second.epsilon, second.residualScale.pressure, ArefP
        );
        if
        (
            first.fd.pressure <= tauFirstP
         && second.fd.pressure <= tauSecondP
         && first.mffd.pressure > tauFirstP
         && second.mffd.pressure > tauSecondP
         && first.fdAdjacentChange[2] <= RelativePlateauTolerance
         && second.fdAdjacentChange[2] <= RelativePlateauTolerance
        )
        {
            pressureDefect = true;
        }
        Info<< "TWO-EPSILON mode=" << mode << " (" << modeName(mode) << ')';
        for (label component = 0; component < 3; ++component)
        {
            const scalar references[] = {ArefFull, ArefD, ArefP};
            const bool pair = suitablePair
            (
                first, second, component, references[component]
            );
            Info<< ' ' << componentName(component) << '='
                << (pair ? "suitable-and-unchanged" : "inconclusive")
                << " locations=" << first.discrepancyCell[component]
                << "/" << second.discrepancyCell[component];
        }
        Info<< endl;
    }

    reportInterpretation(stateLabel, pressureDefect, anyDiscrepancy);

    AssertPETSc(MatDestroy(&mffd));
    AssertPETSc(MatDestroy(&P));
    AssertPETSc(VecDestroy(&baseResidual));
    AssertPETSc(VecDestroy(&residualRepeat));
    AssertPETSc(VecDestroy(&x));
    Info<< "PASS: standalone block-aware MFFD/FD diagnosis completed" << endl;
    return 0;
#endif
}

// ************************************************************************* //
