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

#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fixedGradientFvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"

#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagTotalDispSolid, dictionary
);



// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

namespace
{

// Representability ceiling for the explicit cofactor/determinant arithmetic in
// admissibleDeformation() below. With M = max_ij |F_ij|, each cofactor is a
// difference of two products so |cof| <= 2 M^2, and the determinant is a sum of
// three F*cof terms so |det| <= 6 M^3. Requiring 6 M^3 <= VGREAT therefore
// guarantees that forming cof(F) and det(F) cannot overflow, so FOAM_SIGFPE
// cannot fire inside the guard itself.
//
// This is a floating-point representability bound, not a physical deformation
// limit: the ceiling is ~5e99, far beyond any deformation gradient that is not
// already non-finite.
const scalar maxTrialFComponent = std::pow(VGREAT/6.0, 1.0/3.0);


// Validate a trial deformation gradient and, if it is admissible, return its
// Jacobian and its genuine three-dimensional inverse.
//
// The inverse is formed explicitly as adj(F)/det(F) rather than by calling
// Foam::inv(). The field-level Foam::inv() applies Tensor::safeInv(), whose
// "2-D" branch is a *diagonal-magnitude heuristic*, not a rank or dimension
// test: when any diagonal component is small relative to the diagonal norm it
// adds one to that diagonal, inverts the shifted tensor and subtracts one
// again, returning a finite but incorrect surrogate. A proper 90-degree
// rotation about z, F = ((0 -1 0)(1 0 0)(0 0 1)), is a valid deformation with
// J = 1 for which safeInv() returns ((-0.5 0.5 0)(-0.5 -0.5 0)(0 0 1)) instead
// of the exact inverse F^T. safeInv() also returns the zero tensor whenever
// |det| < ROOTVSMALL. Neither substitution is acceptable for a deformation
// gradient, so this solid model does not use it.
//
// Admissibility rejects only states that are mathematically or numerically
// invalid. No physical deformation limit and no condition-number cap is
// imposed: an orientation-preserving F is accepted however ill-conditioned.
// Note that the max|F_ij| ceiling below is a numerical-arithmetic limit, so
// this is not literally "any finite F": a structured tensor with a component
// above ~5e99 can be mathematically invertible with representable cofactor
// arithmetic and will still be rejected. That conservatism is deliberate --
// the screen must be evaluated before any product is formed -- and it is a
// floating-point bound, not a physical or condition-number cap.
//
// On return, J is the computed determinant when it could be evaluated, and
// -VGREAT as a sentinel when F itself was rejected before that point.
inline bool admissibleDeformation
(
    const tensor& F,
    scalar& J,
    tensor& Finv
)
{
    J = -VGREAT;

    // (1) Every component finite, and small enough that the cofactor and
    //     determinant arithmetic below cannot overflow. std::isfinite performs
    //     no arithmetic and raises no floating-point exception.
    scalar maxF = 0;
    for (direction cmptI = 0; cmptI < tensor::nComponents; ++cmptI)
    {
        if (!std::isfinite(F[cmptI]))
        {
            return false;
        }

        maxF = max(maxF, mag(F[cmptI]));
    }

    if (maxF > maxTrialFComponent)
    {
        return false;
    }

    // (2) Cofactors, and the determinant expanded along the first row
    const scalar cofXX = F.yy()*F.zz() - F.yz()*F.zy();
    const scalar cofXY = F.yz()*F.zx() - F.yx()*F.zz();
    const scalar cofXZ = F.yx()*F.zy() - F.yy()*F.zx();
    const scalar cofYX = F.xz()*F.zy() - F.xy()*F.zz();
    const scalar cofYY = F.xx()*F.zz() - F.xz()*F.zx();
    const scalar cofYZ = F.xy()*F.zx() - F.xx()*F.zy();
    const scalar cofZX = F.xy()*F.yz() - F.xz()*F.yy();
    const scalar cofZY = F.xz()*F.yx() - F.xx()*F.yz();
    const scalar cofZZ = F.xx()*F.yy() - F.xy()*F.yx();

    J = F.xx()*cofXX + F.xy()*cofXY + F.xz()*cofXZ;

    // (3) A non-finite Jacobian, or one that is not orientation preserving.
    //     Strict inequality: J <= 0 is an admissibility condition, not a
    //     tolerance. J < 0 is a locally inverted deformation and the
    //     constitutive law's J^(-2/3) is undefined there.
    if (!std::isfinite(J) || J <= 0)
    {
        return false;
    }

    // (4) Representability of adj(F)/J. Dividing the largest cofactor by
    //     VGREAT cannot overflow, so this test is itself safe, and passing it
    //     guarantees every quotient below is bounded by VGREAT.
    scalar maxCof = mag(cofXX);
    maxCof = max(maxCof, mag(cofXY));
    maxCof = max(maxCof, mag(cofXZ));
    maxCof = max(maxCof, mag(cofYX));
    maxCof = max(maxCof, mag(cofYY));
    maxCof = max(maxCof, mag(cofYZ));
    maxCof = max(maxCof, mag(cofZX));
    maxCof = max(maxCof, mag(cofZY));
    maxCof = max(maxCof, mag(cofZZ));

    if (J < maxCof/VGREAT)
    {
        return false;
    }

    // (5) The genuine 3-D inverse: adj(F)/J, where adj(F) = cof(F)^T
    Finv = tensor
    (
        cofXX/J, cofYX/J, cofZX/J,
        cofXY/J, cofYY/J, cofZY/J,
        cofXZ/J, cofYZ/J, cofZZ/J
    );

    // (6) Belt and braces: the bound at (4) makes this unreachable, but a
    //     non-finite inverse must never be handed downstream
    for (direction cmptI = 0; cmptI < tensor::nComponents; ++cmptI)
    {
        if (!std::isfinite(Finv[cmptI]))
        {
            return false;
        }
    }

    return true;
}

} // End anonymous namespace



void nonLinGeomTotalLagTotalDispSolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    // D() = D().oldTime() + U()*runTime().deltaT();
    D() = D().oldTime() + U()*runTime().deltaT()
        + 0.5*sqr(runTime().deltaT())*A_;

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Jacobian and genuine 3-D inverse of the deformation gradient, with
    // validation (see checkTrialGeometry()). This is not a PETSc callback, so
    // an inadmissible state here is a hard error rather than a recoverable
    // trial rejection
    if (!checkTrialGeometry())
    {
        FatalErrorInFunction
            << "Inadmissible deformation gradient in the predictor: "
            << trialGeometryDescription() << exit(FatalError);
    }

    updateFaceKinematics();

    // Calculate the stress using run-time selectable mechanical law
    if (solvePressure())
    {
        mechanical().correctStressComponents(sigma(), sigmaPreserved_);
    }
    else
    {
        mechanical().correct(sigma());
    }

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        sigma() = dev(sigma()) + sigmaPreserved_ - p()*I;
    }
}


void nonLinGeomTotalLagTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& force,
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        vectorField& forceP = boundaryFieldRef(force)[patchI];

        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // traction.boundaryFieldRef()[patchI] =
            //     tracPatch.traction() - nPatch*tracPatch.pressure();
            if (tracPatch.useUndeformedArea())
            {
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfCurrentPatch;
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<slipFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Set shear traction to zero
            // traction.boundaryFieldRef()[patchI] =
                // sqr(nPatch) & traction.boundaryField()[patchI];
            forceP = sqr(nPatch) & force.boundaryField()[patchI];
        }
    }
}


bool nonLinGeomTotalLagTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    if (predictor_ && newTimeStep())
    {
        predict();
    }

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    int iCorr = 0;
    scalar currentResidualNorm = 0;
    scalar initialResidualNorm = 0;
    scalar deltaXNorm = 0;
    scalar xNorm = 0;
    const convergenceParameters convParam =
        readConvergenceParameters(solidModelDict());

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Undeformed unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Calculate deformed area vectors and normals
        const surfaceVectorField SfCurrent
        (
            fvc::interpolate(J_*Finv_.T()) & mesh().Sf()
        );
        const surfaceScalarField magSfCurrent(mag(SfCurrent));
        const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

        // Traction vectors at the faces
        surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

        // Add stabilisation to the traction
        // We add this before enforcing the traction condition as the stabilisation
        // is set to zero on traction boundaries
        momentumStabilisation().updateVector(D(), &gradD());
        traction += impKf_*momentumStabilisation().faceVector();

        // Calculate the force at the faces
        surfaceVectorField force(magSfCurrent*traction);

        // Enforce traction boundary conditions
        enforceTractionBoundaries(force, D(), nCurrent, magSfCurrent);

        // Momentum equation total displacement total Lagrangian form
#ifndef OPENFOAM_COM
        // Assemble the RHS in stages.
        // The equivalent chained tmp fvMatrix expression is stable on OpenFOAM.com.
        tmp<fvVectorMatrix> tRhsEqn
        (
            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
        );
        tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
        tmpRef(tRhsEqn) += fvc::div(force);
        tmpRef(tRhsEqn) += rho()*g();

        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == tRhsEqn
        );
#else
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(force)
          + rho()*g()
        );
#endif

        // Add damping
        if (dampingCoeff().value() > SMALL)
        {
            DEqn += dampingCoeff()*rho()*fvm::ddt(D());
        }

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system and store the residual
        currentResidualNorm = mag(DEqn.solve().initialResidual());

        // Norm of the solution correction
        deltaXNorm =
            sqrt
            (
                gSum
                (
                    magSqr
                    (
#ifdef OPENFOAM_NOT_EXTEND
                        D().primitiveField() - D().prevIter().primitiveField()
#else
                        D().internalField() - D().prevIter().internalField()
#endif
                    )
                )
            );

        // Norm of the solution
#ifdef OPENFOAM_NOT_EXTEND
        xNorm = sqrt(gSum(magSqr(D().primitiveField())));
#else
        xNorm = sqrt(gSum(magSqr(D().internalField())));
#endif

        // Store the initial residual
        if (iCorr == 0)
        {
            initialResidualNorm = currentResidualNorm;
            Info<< "Initial Residual Norm = " << initialResidualNorm << nl
                << "Initial Solution Norm = " << xNorm << endl;
        }

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Jacobian and genuine 3-D inverse, with validation: see
        // checkTrialGeometry()
        if (!checkTrialGeometry())
        {
            FatalErrorInFunction
                << "Inadmissible deformation gradient at iteration " << iCorr
                << ": " << trialGeometryDescription() << exit(FatalError);
        }

        updateFaceKinematics();

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
        !checkConvergence
        (
            currentResidualNorm,
            initialResidualNorm,
            deltaXNorm,
            xNorm,
            ++iCorr,
            convParam
        )
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool nonLinGeomTotalLagTotalDispSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Seed the PETSc solution vector from the predicted fields
        packSolution(foamPetscSnesHelper::solution());
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // foamPetscSnesHelper::solve() treats SNES_DIVERGED_FUNCTION_DOMAIN as
    // a soft, non-fatal outcome by design, so that a single rejected
    // trial iterate (formResidual() detecting a locally inverted
    // deformation) lets PETSc's line search back off gracefully within a
    // Newton solve. That is a different situation from the whole step
    // genuinely failing to converge to any valid state, which solve()'s
    // return value alone does not distinguish. Check the actual final
    // SNES convergence reason here so a fully diverged step is reported
    // loudly instead of the diverged solution being silently unpacked
    // and accepted as this time step's answer.
    SNESConvergedReason snesReason = SNES_CONVERGED_ITERATING;
    const PetscErrorCode reasonErr =
        SNESGetConvergedReason(snes(), &snesReason);

    if (reasonErr != 0)
    {
        // A genuine PETSc API failure, not a convergence outcome
        FatalErrorInFunction
            << "SNESGetConvergedReason failed with PETSc error code "
            << label(reasonErr) << exit(FatalError);
    }

    if (snesReason < 0)
    {
        // The nonlinear solve did not converge. A diverged PETSc vector is
        // never unpacked or accepted as this time step's answer, whatever
        // stopOnPetscError is set to; what the switch selects is only whether
        // that is fatal.
        //
        // Restore the last accepted state so the model is left in the
        // well-defined configuration it had at the start of the step rather
        // than at whatever trial the failed solve happened to evaluate last.
        volVectorField& D = const_cast<volVectorField&>(this->D());
        D = D.oldTime();

        if (solvePressure())
        {
            volScalarField& p = const_cast<volScalarField&>(this->p());
            p = p.oldTime();
        }

        packSolution(foamPetscSnesHelper::solution());

        // Restored state: refresh the face state too, so anything written for
        // this failed step is consistent with the fields actually held
        unpackSolution(foamPetscSnesHelper::solution(), true);

        // Refuse to shift the accepted history in updateTotalFields()
        solveFailed_ = true;

        // Feed the existing end-of-run convergence report, so that end() can
        // no longer claim "The momentum equation converged in all time-steps"
        maxIterReached()++;

        if (stopOnPetscError())
        {
            FatalErrorInFunction
                << "The PETSc SNES solver failed to converge this time step "
                << "(SNESConvergedReason = "
                << SNESConvergedReasons[snesReason]
                << "); refusing to accept a diverged solution as the answer "
                << "for time = " << mesh().time().timeName() << nl
                << "Set `stopOnPetscError` to `false` to report the failure "
                << "and continue instead of stopping here"
                << exit(FatalError);
        }

        WarningInFunction
            << "The PETSc SNES solver failed to converge this time step "
            << "(SNESConvergedReason = " << SNESConvergedReasons[snesReason]
            << ") at time = " << mesh().time().timeName() << nl
            << "    The diverged solution has NOT been accepted: the fields "
            << "have been restored to the last accepted state and the "
            << "accumulated history will not be advanced." << nl
            << "    evolve() returns false. Note that the solids4Foam "
            << "application does not currently act on that return value, so "
            << "the run will continue to the next time step from the restored "
            << "state." << endl;

        return false;
    }

    solveFailed_ = false;

    // Map the PETSc solution back into the D field (and p when active),
    // refreshing dependent kinematic fields and boundary conditions.
    //
    // This is the accepted state, and the only unpack that rebuilds the face
    // kinematics: it runs before the law's face constitutive fields are
    // recomputed, and the application calls updateTotalFields() after it, which
    // shifts Ef_ into EfOld_/EfOldOld_. The accepted face history is therefore
    // built from the accepted Ff
    unpackSolution(foamPetscSnesHelper::solution(), true);

    if (solvePressure())
    {
        // Update dpdt
        autoPtrRef(dpdtPtr_) = fvc::ddt(p());
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagTotalDispSolid::nonLinGeomTotalLagTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "D",
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::CELLS,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    // Finv_ and J_ are given harmless placeholders here and are populated in
    // the constructor body through checkTrialGeometry(), the single validated
    // path. They must NOT be initialised with Foam::inv(F_)/Foam::det(F_):
    // member initialisers run before the constructor body, so those would
    // evaluate the generic OpenFOAM inverse/determinant on an unvalidated F_.
    // For an F_ read from a time directory that is finite but extreme (e.g.
    // diag(1e200, 1e-200, 1)), Tensor::safeInv() squares the diagonal
    // components, overflowing to Inf and raising FE_OVERFLOW under FOAM_SIGFPE
    // before any check could run; and for a valid 90-degree rotation it would
    // return an incorrect surrogate inverse.
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    trialGeometryValid_(true),
    trialMinJ_(1.0),
    trialRejectionCount_(0),
    faceSubstitutionCount_(0),
    faceRefreshCount_(0),
    faceRefreshSkipCount_(0),
    warnedUncorrectedPatches_(false),
    solveFailed_(false),
    sigmaPreserved_
    (
        IOobject
        (
            "sigmaPreserved",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    A_
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::d2dt2(D())
    ),
    impK_
    (
        solvePressure()
      ? 2.0*mechanical().shearModulus()
      : mechanical().impK()
    ),
    impKf_(fvc::interpolate(impK_)),
    rImpK_(1.0/impK_),
    viscousImpKfPtr_(),
    viscousImpKfDeltaT_(0),
    viscousJacobian_
    (
        solidModelDict().lookupOrDefault<Switch>("viscousJacobian", true)
    ),
    faceForceTreatment_
    (
        solidModelDict().lookupOrDefault<word>
        (
            "faceForceTreatment", "interpolatedCell"
        )
    ),
    faceForceDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>("faceForceDiagnostics", false)
    ),
    directFaceFallbackCount_(0),
    rKappaPtr_(),
    dpdtPtr_(),
    rAUfTimeIndex_(-1),
    rAUfDeltaT_(0),
    scaleMixedPetScFields_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "scaleMixedPetScFields", true
        )
    ),
    pressureUnknownScaleType_
    (
        solidModelDict().lookupOrDefault<word>("pressureUnknownScale", "twoMu")
    ),
    pressureUnknownScale_(1.0),
    pressureScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("pressureScaleFactor", 1.0)
    ),
    pressureScaleByTwoMu_
    (
        solidModelDict().lookupOrDefault<Switch>("pressureScaleByTwoMu", true)
    ),
    twoMuRef_(1.0),
    pressureEqnScale_(pressureScaleFactor_),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    )
{
    DisRequired();

    // The Finv_ and J_ member initialisers above use Foam::inv()/Foam::det()
    // for convenience. Foam::inv() applies Tensor::safeInv(), which can return
    // a finite but incorrect surrogate (see admissibleDeformation()), so
    // recompute both here through the validated path. F_ is the identity unless
    // it was read from a time directory.
    if (!checkTrialGeometry())
    {
        FatalErrorInFunction
            << "Inadmissible initial deformation gradient: "
            << trialGeometryDescription() << exit(FatalError);
    }

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields.
    // updateFaceKinematics() must run first: correct() refreshes the law's
    // face viscous fields from Ff(), and Ff() is otherwise still the identity
    // it was created with (or a stale identity read from a time directory)
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        updateFaceKinematics();

        mechanical().correct(sigma());
    }

    Info<< "solvePressure = " << solvePressure() << endl;
    if (solvePressure())
    {
        if (solutionAlg() != solutionAlgorithm::PETSC_SNES)
        {
            FatalErrorInFunction
                << "The solution algorithm must be "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " when solvePressure is enabled" << abort(FatalError);
        }

        // Ensure p is created and the oldTime is stored
        p().oldTime();

        // Initialise dpdt field
        dpdtPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(p())
            )
        );

        // Use a volume-weighted average of 2*mu as the physical scale
        // of the pressure equation. The pressure-row residual and
        // Jacobian are then multiplied by
        // pressureEqnScale_ = pressureScaleFactor_ * twoMuRef_ so that
        // their natural magnitude is comparable to the momentum block.
        const volScalarField twoMu(2.0*mechanical().shearModulus());
        scalar twoMuV = 0;
        scalar Vtot = 0;
        forAll(twoMu, cellI)
        {
            const scalar Vc = mesh().V()[cellI];
            twoMuV += twoMu[cellI]*Vc;
            Vtot += Vc;
        }
        reduce(twoMuV, sumOp<scalar>());
        reduce(Vtot, sumOp<scalar>());
        twoMuRef_ = twoMuV/Vtot;
        pressureEqnScale_ =
            pressureScaleFactor_*(pressureScaleByTwoMu_ ? twoMuRef_ : 1.0);

        pressureUnknownScale_ = 1.0;
        if (scaleMixedPetScFields_)
        {
            if
            (
                pressureUnknownScaleType_ == "twoMu"
             || pressureUnknownScaleType_ == "2mu"
            )
            {
                pressureUnknownScale_ = twoMuRef_;
            }
            else if
            (
                pressureUnknownScaleType_ == "user"
             || pressureUnknownScaleType_ == "scalar"
            )
            {
                pressureUnknownScale_ =
                    readScalar
                    (
                        solidModelDict().lookup("pressureUnknownScaleValue")
                    );
            }
            else if (pressureUnknownScaleType_ == "none")
            {
                pressureUnknownScale_ = 1.0;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown pressureUnknownScale "
                    << pressureUnknownScaleType_ << nl
                    << "Valid options are twoMu, user, scalar, none"
                    << abort(FatalError);
            }

            if (pressureUnknownScale_ <= VSMALL)
            {
                FatalErrorInFunction
                    << "pressureUnknownScale must be positive, found "
                    << pressureUnknownScale_ << abort(FatalError);
            }
        }

        Info<< "pressureEqnScale = " << pressureEqnScale_
            << ", where pressureScaleFactor = " << pressureScaleFactor_
            << " and 2*mu = " << twoMuRef_ << endl;

        Info<< "PETSc pressure unknown scale = " << pressureUnknownScale_
            << " (scaleMixedPetScFields = " << scaleMixedPetScFields_
            << ", pressureUnknownScale = " << pressureUnknownScaleType_
            << ")" << endl;
    }

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();

        if (!checkTrialGeometry())
        {
            FatalErrorInFunction
                << "Inadmissible deformation gradient read on restart: "
                << trialGeometryDescription() << exit(FatalError);
        }

        // Reconstruct Ff from the restarted displacement rather than trusting
        // any identity Ff_ written to the time directory
        updateFaceKinematics();

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }

    // Check the gradScheme
    const word gradDScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + D().name() +')')
#else
        mesh().schemesDict().gradScheme("grad(" + D().name() +')')
#endif
    );

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        if (gradDScheme != "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithm" << abort(FatalError);
        }

        // Set extrapolateValue to true for solidTraction boundaries
        forAll(D().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
                        boundaryFieldRef(D())[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


void nonLinGeomTotalLagTotalDispSolid::makeRKappa() const
{
    if (rKappaPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rKappaPtr_.set(new volScalarField(1.0/mechanical().bulkModulus()));
}


const volScalarField& nonLinGeomTotalLagTotalDispSolid::rKappa() const
{
    if (rKappaPtr_.empty())
    {
        makeRKappa();
    }

    return rKappaPtr_();
}


Foam::string nonLinGeomTotalLagTotalDispSolid::trialGeometryDescription() const
{
    // trialMinJ_ holds the smallest determinant actually evaluated. When F
    // itself was rejected the determinant was never formed and the sentinel
    // -VGREAT is reported instead, so say what really failed rather than
    // printing a meaningless "min(J)"
    if (trialMinJ_ <= -VGREAT)
    {
        return string
        (
            "F has a non-finite component, or one whose magnitude is outside "
            "the range for which the cofactor/determinant arithmetic is "
            "representable"
        );
    }

    return "min(J) = " + ::Foam::name(trialMinJ_);
}


void nonLinGeomTotalLagTotalDispSolid::refreshFaceKinematics()
{
    updateFaceKinematics();
}


void nonLinGeomTotalLagTotalDispSolid::updateFaceKinematics()
{
    faceRefreshCount_++;

    // Vertex displacements, using the same gradient-corrected interpolation
    // this model already uses elsewhere so that pointD has a single definition
    mechanical().interpolate(D(), gradD(), pointD());

    // Face displacement gradient: fvc::fGrad(D, pointD), the established
    // solids4foam face-gradient operator, which fills internal and boundary
    // patch values
    mechanical().grad(D(), pointD(), gradDf_);

    // Make the traction-boundary face gradient purely kinematic.
    //
    // fvc::fGrad builds the in-plane part from the point values and adds the
    // normal part as n*fvc::snGrad(D) (fvcGradf.C:108), so the normal row of
    // the face gradient is whatever the boundary condition reports as snGrad.
    // On a fixedGradient displacement patch that is not a kinematic
    // derivative: tractionBoundarySnGrad() returns
    //
    //     (n & gradD) + (applied traction - internal traction)/impK
    //
    // The second term belongs to the finite-volume traction enforcement
    // iteration. It is proportional to a boundary force residual divided by the
    // arbitrary implicit stiffness impK (the representative-shear-modulus
    // surrogate, 550 Pa here), and it must not become constitutive deformation:
    // at t = 0.02 in Case B it made the endocardial face strain 0.134 against an
    // internal-face maximum of 0.0024, and that contaminated state was being
    // committed into the viscous face history EfOld_/EfOldOld_.
    //
    // The correction keeps the accurate point-based tangential reconstruction
    // and replaces only the normal row with the uncontaminated one-sided cell
    // gradient. With OpenFOAM's convention gradD_ij = d(D_j)/d(X_i), the normal
    // derivative is (n & G) and the normal row is n*(n & G), so
    //
    //     Gkin = Gf + n*(n & (Gc - Gf))
    //
    // gives exactly
    //
    //     n & Gkin           = n & Gc          (uncontaminated normal part)
    //     (I - n*n) & Gkin   = (I - n*n) & Gf  (tangential part untouched)
    //
    // Gc is gradD().boundaryField(), which leastSquaresS4fGrad extrapolates from
    // cell values with useBoundaryFaceValues = false
    // (leastSquaresS4fGrad.C:119), so it carries no boundary-condition term.
    //
    // Applied only where snGrad is an imposed gradient rather than a kinematic
    // one. fixedGradient is the precise capability test: solidTraction and every
    // BC derived from it, plus tractionPressureDisplacement, are
    // fixedGradientFvPatchVectorField, and they are exactly the classes that
    // obtain their gradient from tractionBoundarySnGrad(). fixedValue-type,
    // symmetry, empty/wedge and coupled patches are excluded, because their
    // snGrad IS a kinematic derivative and fvc::fGrad's value there is the
    // better one. The solidDirectionMixed family (normalDisplacementZeroShear,
    // fixedDisplacementZeroShear, displacementOrTraction) is traction-driven
    // only in the directions where its valueFraction vanishes; correcting the
    // whole normal row would corrupt its Dirichlet directions, so it is left
    // alone and reported once
    forAll(gradDf_.boundaryField(), patchI)
    {
        const fvPatchVectorField& pD = D().boundaryField()[patchI];

        if (!isA<fixedGradientFvPatchVectorField>(pD))
        {
            if
            (
                !warnedUncorrectedPatches_
             && isA<solidDirectionMixedFvPatchVectorField>(pD)
            )
            {
                warnedUncorrectedPatches_ = true;

                WarningInFunction
                    << "Patch " << mesh().boundary()[patchI].name()
                    << " is a solidDirectionMixed-family displacement condition."
                    << " Its face deformation gradient is traction-driven only"
                    << " in the directions where valueFraction vanishes, so the"
                    << " normal-row correction is not applied there and its"
                    << " face Ff may retain a traction/impK contribution."
                    << endl;
            }

            continue;
        }

        const vectorField n(mesh().boundary()[patchI].nf());
        const tensorField& Gc = gradD().boundaryField()[patchI];
        tensorField& Gf = gradDf_.boundaryFieldRef()[patchI];

        Gf += n*(n & (Gc - Gf));
    }

    // Final safety guard, retained from the previous revision but demoted: it is
    // NOT the mechanism that removes traction contamination (the normal-row
    // correction above is), and it is not a definition of authoritative face
    // kinematics. It exists only so that a face which is still not a deformation
    // cannot reach the constitutive routines, where J <= 0 is a FatalError.
    //
    // This branch is discontinuous: Ff, Jf, Ef, Edotf and any face stress jump
    // when a face crosses the admissibility boundary. It must therefore never be
    // treated as part of a smooth constitutive map, and a future exact face
    // tangent must either differentiate the active branch or refuse to claim
    // exactness while any substitution is present. A healthy Case B step
    // requires zero substitutions
    label nSubstituted = 0;
    scalar worstJ = VGREAT;
    labelList nSubPerPatch(mesh().boundary().size(), 0);

    {
        const labelUList& own = mesh().owner();
        const tensorField& gradDI = gradD().primitiveField();
        tensorField& gradDfI = gradDf_.primitiveFieldRef();

        forAll(gradDfI, faceI)
        {
            scalar Jf = -VGREAT;
            tensor Finvf(tensor::zero);

            if (!admissibleDeformation(I + gradDfI[faceI].T(), Jf, Finvf))
            {
                gradDfI[faceI] = gradDI[own[faceI]];
                nSubstituted++;
                worstJ = min(worstJ, Jf);
            }
        }
    }

    forAll(gradDf_.boundaryField(), patchI)
    {
        const tensorField& gradDP = gradD().boundaryField()[patchI];
        tensorField& gradDfP = gradDf_.boundaryFieldRef()[patchI];

        forAll(gradDfP, faceI)
        {
            scalar Jf = -VGREAT;
            tensor Finvf(tensor::zero);

            if (!admissibleDeformation(I + gradDfP[faceI].T(), Jf, Finvf))
            {
                gradDfP[faceI] = gradDP[faceI];
                nSubstituted++;
                worstJ = min(worstJ, Jf);
                nSubPerPatch[patchI]++;
            }
        }
    }

    reduce(nSubstituted, sumOp<label>());
    reduce(worstJ, minOp<scalar>());
    forAll(nSubPerPatch, patchI)
    {
        reduce(nSubPerPatch[patchI], sumOp<label>());
    }

    if (nSubstituted > 0)
    {
        // Never silent
        faceSubstitutionCount_ += nSubstituted;

        Info<< "    Face deformation gradient: " << nSubstituted
            << " face(s) still inadmissible after the kinematic boundary"
            << " correction (min J = " << worstJ
            << "); the cell displacement gradient was used there instead"
            << endl;

        forAll(nSubPerPatch, patchI)
        {
            if (nSubPerPatch[patchI] > 0)
            {
                Info<< "        patch " << mesh().boundary()[patchI].name()
                    << ": " << nSubPerPatch[patchI] << " face(s)" << endl;
            }
        }
    }

    // Hand it to the mechanical law: Ff() = I + gradDf^T
    mechanical().updateFf();
}


bool nonLinGeomTotalLagTotalDispSolid::checkTrialGeometry()
{
    // Validate F_ and, where admissible, fill J_ and Finv_ with the true
    // determinant and the genuine 3-D inverse. Foam::det() and Foam::inv() are
    // deliberately not used here: see admissibleDeformation() above.
    //
    // Internal AND boundary values are covered. Both are consumed downstream:
    // the mechanical law evaluates the constitutive response on every boundary
    // patch, fvc::interpolate(J_*Finv_.T()) & Sf forms the Nanson current-area
    // vector using the boundary values of J_ and Finv_ on boundary faces, and
    // tractionBoundarySnGrad() reads Finv_.boundaryField() directly.
    //
    // When this returns false the contents of J_ and Finv_ must not be used;
    // unpackSolution() returns immediately and formResidual()/formJacobian()
    // both bail out on the trialGeometryValid_ flag.
    bool valid = true;
    scalar minJ = VGREAT;

    // Internal cells
    {
        const tensorField& FI = F_.primitiveField();
        scalarField& JI = primitiveFieldRef(J_);
        tensorField& FinvI = primitiveFieldRef(Finv_);

        forAll(FI, cellI)
        {
            scalar Jc = -VGREAT;
            tensor Finvc(tensor::zero);

            if (admissibleDeformation(FI[cellI], Jc, Finvc))
            {
                JI[cellI] = Jc;
                FinvI[cellI] = Finvc;
            }
            else
            {
                valid = false;
            }

            minJ = min(minJ, Jc);
        }
    }

    // Boundary patches
    forAll(F_.boundaryField(), patchI)
    {
        const tensorField& FP = F_.boundaryField()[patchI];
        scalarField& JP = boundaryFieldRef(J_)[patchI];
        tensorField& FinvP = boundaryFieldRef(Finv_)[patchI];

        forAll(FP, faceI)
        {
            scalar Jc = -VGREAT;
            tensor Finvc(tensor::zero);

            if (admissibleDeformation(FP[faceI], Jc, Finvc))
            {
                JP[faceI] = Jc;
                FinvP[faceI] = Finvc;
            }
            else
            {
                valid = false;
            }

            minJ = min(minJ, Jc);
        }
    }

    reduce(valid, andOp<bool>());
    reduce(minJ, minOp<scalar>());

    trialMinJ_ = minJ;

    return valid;
}


void nonLinGeomTotalLagTotalDispSolid::updateRAUfIfStale()
{
    if (!solvePressure())
    {
        // rAUf is unused unless we are solving the mixed system
        return;
    }

    const label tIdx = runTime().timeIndex();
    const scalar dt = runTime().deltaT().value();

    // rAUf depends on (impKf_, rho, mesh, deltaT) only. It is fresh
    // when both the timeIndex and deltaT match what we cached
    if
    (
        rAUfTimeIndex_ >= 0
     && rAUfTimeIndex_ == tIdx
     && mag(dt - rAUfDeltaT_) <= SMALL*max(mag(dt), SMALL)
    )
    {
        return;
    }

    // Build the approximate momentum diagonal. fvm::laplacian and
    // fvm::d2dt2 read only the BC structure of D, not its values, so
    // the resulting diagonal is independent of the current Newton
    // iterate and any MFFD perturbation
    fvVectorMatrix approxMomJ
    (
        fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D())
    );
    approxMomJ.relax();
    rAUf() = -1.0/(fvc::interpolate(approxMomJ.A()));

    rAUfTimeIndex_ = tIdx;
    rAUfDeltaT_ = dt;
}


const surfaceScalarField*
nonLinGeomTotalLagTotalDispSolid::viscousImpKfIfAny()
{
    if (!viscousJacobian_)
    {
        return nullptr;
    }

    const scalar dt = runTime().deltaT().value();

    // Stale exactly when deltaT changes: independent of D, p and the Newton
    // iterate, so it survives every MFFD probe untouched
    if
    (
        viscousImpKfPtr_.valid()
     && mag(dt - viscousImpKfDeltaT_) <= SMALL*max(mag(dt), SMALL)
    )
    {
        return &viscousImpKfPtr_();
    }

    tmp<surfaceScalarField> tvisc(mechanical().viscousImpKf());
    viscousImpKfDeltaT_ = dt;

    // Rate-independent law or zero viscosity: nothing to add
    if (max(mag(tvisc())).value() <= SMALL)
    {
        viscousImpKfPtr_.clear();
        return nullptr;
    }

    viscousImpKfPtr_.reset(new surfaceScalarField(tvisc));

    return &viscousImpKfPtr_();
}


#ifdef USE_PETSC

void nonLinGeomTotalLagTotalDispSolid::unpackSolution
(
    const Vec x,
    const bool refreshFaceState
)
{
    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions on D
    D.correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Enforce the boundary conditions again for any conditions that use gradD
    D.correctBoundaryConditions();

    // Increment of displacement
    DD() = D - D.oldTime();

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Total deformation gradient
    F_ = I + gradD().T();

    // Validate the trial deformation gradient and, if it is admissible, fill
    // J_ and Finv_ with the determinant and the genuine 3-D inverse. Nothing
    // downstream of this point sees an unvalidated F_, J_ or Finv_, on either
    // internal cells or boundary patches. See checkTrialGeometry() and
    // admissibleDeformation().
    trialGeometryValid_ = checkTrialGeometry();

    if (!trialGeometryValid_)
    {
        // Leave Finv_, sigma and p at their previous values: they are not
        // read on this path. formResidual() and formJacobian() check the
        // flag and hand the rejection to PETSc's domain-error mechanism,
        // which tags the residual vector so the line search backtracks.
        trialRejectionCount_++;

        Info<< "    Rejecting SNES trial state: " << trialGeometryDescription()
            << " is not an admissible deformation (rejection "
            << trialRejectionCount_ << ")" << endl;

        return;
    }

    // Face kinematics are rebuilt only when something will consume them.
    //
    // The momentum/pressure residual is built entirely from cell quantities
    // (fvc::interpolate(sigma()) and fvc::interpolate(J_*Finv_.T())), so an
    // ordinary MFFD residual probe never reads Ff. Rebuilding it there costs a
    // volPointInterpolation plus an fvc::fGrad per probe, which over the ~21k
    // (mesh1) and ~53k (mesh3) KSP iterations of a 400-step run was 27-42% of
    // total wall time for a field the probe cannot observe.
    //
    // The face state must be current at exactly one point per accepted step:
    // immediately before the law's face constitutive fields are recomputed in
    // the accepted unpack, because updateTotalFields() then shifts Ef_ into
    // EfOld_/EfOldOld_. That ordering is preserved here -- the refresh runs
    // BEFORE mechanical().correctStressComponents() below, so the accepted Ef_
    // is built from the accepted Ff, not from a stale one.
    //
    // A future face tangent must call refreshFaceKinematics() explicitly at the
    // Jacobian state it linearises about; see that function's documentation
    if (refreshFaceState)
    {
        updateFaceKinematics();
    }
    else
    {
        faceRefreshSkipCount_++;
    }

    // Calculate the stress using run-time selectable mechanical law
    if (solvePressure())
    {
        mechanical().correctStressComponents(sigma(), sigmaPreserved_);
    }
    else
    {
        mechanical().correct(sigma());
    }

    if (solvePressure())
    {
        // Copy the scaled pressure unknown pHat from x into the
        // physical p field via p = pressureUnknownScale_ * pHat
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        scalarField pHat(pI.size(), 0.0);
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pHat, blockSize_ - 1
        );
        pI = pressureUnknownScale_*pHat;

        // Enforce the boundary conditions on p
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) + sigmaPreserved_ - p*I;
    }
}


void nonLinGeomTotalLagTotalDispSolid::packSolution(Vec x)
{
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        primitiveField(D()),
        x,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    if (solvePressure())
    {
        // Insert the scaled pressure unknown pHat = p/pressureUnknownScale_
        scalarField pHat(primitiveField(p()));
        pHat /= pressureUnknownScale_;
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pHat,
            x,
            blockSize_ - 1
        );
    }
}

#endif // USE_PETSC


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagTotalDispSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    // {
    //     // Not yet implmented, although coupledUnsLinGeomLinearElasticSolid
    //     // could be combined with PETSc to achieve this.. todo!
    //     return evolveImplicitCoupled();
    // }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    // else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    // {
    //     return evolveExplicit();
    // }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
               ]
            // << solidModel::solutionAlgorithmNames_
            //    [
            //        solidModel::solutionAlgorithm::EXPLICIT
            //    ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


#ifdef USE_PETSC

label nonLinGeomTotalLagTotalDispSolid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label nonLinGeomTotalLagTotalDispSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


void nonLinGeomTotalLagTotalDispSolid::reportFaceForceDifference
(
    const surfaceVectorField& forceA,
    const surfaceVectorField& forceB
) const
{
    // Largest relative difference between the interpolated-cell (A) and
    // direct-face (B) internal force, and where it occurs. Diagnostic only.
    scalar worstRel = 0;
    scalar worstMagA = 0;
    scalar worstMagB = 0;
    label worstFace = -1;
    label worstOwn = -1;
    label worstNei = -1;
    vector worstCf(vector::zero);

    const scalar scale = max(gMax(mag(forceA.primitiveField())), SMALL);
    const labelUList& own = mesh().owner();
    const labelUList& nei = mesh().neighbour();

    forAll(forceA.primitiveField(), faceI)
    {
        const scalar d = mag(forceA[faceI] - forceB[faceI]);
        const scalar rel = d/scale;

        if (rel > worstRel)
        {
            worstRel = rel;
            worstMagA = mag(forceA[faceI]);
            worstMagB = mag(forceB[faceI]);
            worstFace = faceI;
            worstOwn = own[faceI];
            worstNei = nei[faceI];
            worstCf = mesh().Cf()[faceI];
        }
    }

    scalar l2 = 0;
    forAll(forceA.primitiveField(), faceI)
    {
        l2 += magSqr(forceA[faceI] - forceB[faceI]);
    }
    label nFaces = forceA.primitiveField().size();
    reduce(l2, sumOp<scalar>());
    reduce(nFaces, sumOp<label>());
    l2 = Foam::sqrt(l2/max(nFaces, 1));

    reduce(worstRel, maxOp<scalar>());

    Info<< "    faceForce diagnostics: max rel diff " << worstRel
        << ", rms diff " << l2 << " N"
        << " (scale " << scale << " N)" << nl
        << "      worst internal face " << worstFace
        << "  owner " << worstOwn << "  neighbour " << worstNei
        << "  at " << worstCf << nl
        << "      |interpolatedCell| " << worstMagA
        << "   |directFaceState| " << worstMagB << endl;
}


label nonLinGeomTotalLagTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Copy x into the D field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    unpackSolution(x);

    // If this trial iterate is not an admissible deformation state,
    // unpackSolution() has already skipped inv(F_) and the mechanical law.
    // Report it as a recoverable domain rejection so PETSc tags the
    // residual vector and the line search backtracks, rather than
    // proceeding to use partially-stale fields below. This is an expected
    // event during a line search, not a programming error.
    if (!trialGeometryValid_)
    {
        return foamPetscSnesHelper::SNES_EVAL_DOMAIN_ERROR;
    }

    // Take a non-const reference to D for local use below
    volVectorField& D = const_cast<volVectorField&>(this->D());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above, its BCs corrected and the pressure component of sigma
        // updated; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        // Calculate the pressure gradient
        const volVectorField gradp(fvc::grad(p));

        // Re-calculate the pressure stabilisation parameter
        pressureStabilisation().updateScalar(p, &gradp);

        // Refresh rAUf (the positive face-interpolated reciprocal of
        // the approximate momentum equation diagonal -- the solid
        // analogue of rAUf in pressure-velocity coupling, units [m^2/Pa])
        // only when the mesh or deltaT have changed. The diagonal is
        // value-independent of D and p, so this is safe under PETSc
        // matrix-free finite-difference perturbations
        updateRAUfIfStale();

        // Calculate pressure equation residual. Keep the finite-strain
        // volumetric term -0.5*(J^2-1)/J as in the published nonlinear
        // total-Lagrangian formulation; its linearisation about D=0
        // yields -V*div(D) which is what InsertFvmDivUIntoPETScMatrix
        // assembles in formJacobian.
        scalarField pressureResidual
        (
          - p*rKappa()
          + pressureStabilisation().cellScalar(&rAUf(), true)
          - 0.5*(pow(J_, 2.0) - 1.0)/J_
        );

        // Make residual extensive
        pressureResidual *= mesh.V();

        // Apply the physical row scaling. pressureEqnScale_ already
        // bakes in both the user-facing pressureScaleFactor and the
        // 2*mu physical scale.
        if (pressureEqnScale_ != 1.0)
        {
            pressureResidual *= pressureEqnScale_;
        }

        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
        );
    }

    // ---- face state for the internal momentum force -----------------------
    //
    // Two constructions are available, selected by faceForceTreatment:
    //
    //   interpolatedCell (default): interpolate the CELL stress and the CELL
    //     J*Finv^T to the faces.  This is the original implementation and is
    //     preserved here byte-for-byte.
    //
    //   directFaceState (experimental): build the face state directly from
    //     the corrected face kinematics gradDf_, so no cell->face
    //     interpolation of either the stress or the geometry is involved.
    //
    // Only the internal face force differs.  Traction-boundary enforcement,
    // momentum stabilisation, the pressure residual, the constitutive law,
    // the boundary conditions, the time scheme and the preconditioner are
    // identical in both paths.

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    // Traction vectors at the faces
    //surfaceVectorField traction(n & fvc::interpolate(sigma()));
    surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

    const bool useDirectFace = (faceForceTreatment_ == "directFaceState");

    if
    (
        !useDirectFace
     && faceForceTreatment_ != "interpolatedCell"
    )
    {
        FatalErrorInFunction
            << "Unknown faceForceTreatment " << faceForceTreatment_ << nl
            << "Valid options are interpolatedCell and directFaceState"
            << abort(FatalError);
    }

    autoPtr<surfaceScalarField> magSfDirectPtr;
    autoPtr<surfaceVectorField> tractionDirectPtr;
    autoPtr<surfaceVectorField> forceDirectFacePtr;

    if (useDirectFace || faceForceDiagnostics_)
    {
        // Every PETSc residual probe must see face kinematics consistent with
        // the trial D it was given, so the face state is refreshed here.  This
        // reuses the existing routine, and therefore keeps the fixedGradient
        // traction-boundary normal-row correction and the counted
        // inadmissible-face owner-cell fallback exactly as they are.
        updateFaceKinematics();

        // Material face stress through the mechanical law's own surface
        // interface, i.e. the same corrected stress split used for the cells
        surfaceSymmTensorField sigmaFace
        (
            IOobject
            (
                "sigmaFaceDirect",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        );

        {
            surfaceSymmTensorField sigmaFacePreserved(sigmaFace);
            mechanical().correctStressComponents(sigmaFace, sigmaFacePreserved);
            sigmaFace = dev(sigmaFace) + sigmaFacePreserved;
        }

        // Mixed pressure, added as -I_h[p]*I exactly as the cell path does
        if (solvePressure())
        {
            sigmaFace -= fvc::interpolate(this->p())*symmTensor::I;
        }

        // Current face area vector from the DIRECT face deformation gradient,
        // Jf*Ff^-T & Sf, rather than from an interpolated cell J*Finv^T.
        // admissibleDeformation() is the same explicit screen the cell state
        // uses; where a face fails it, the owner-cell state is substituted and
        // the event is counted and reported rather than silently absorbed.
        surfaceVectorField SfDirect
        (
            IOobject
            (
                "SfDirect",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimArea, vector::zero)
        );

        // Sf carries OpenFOAM's "oriented" attribute, and every quantity
        // derived from it downstream (nCurrent, the momentum stabilisation
        // face vector, ...) inherits it.  SfDirect plays exactly the same
        // role, so it must be marked the same way or the arithmetic below
        // fails the oriented/unoriented type check.
        SfDirect.setOriented(true);

        label nFallback = 0;

        {
            const labelUList& own = mesh.owner();
            const tensorField& gradDfI = gradDf_.primitiveField();
            const vectorField& SfI = mesh.Sf().primitiveField();
            vectorField& SfDirectI = SfDirect.primitiveFieldRef();

            forAll(SfDirectI, faceI)
            {
                scalar Jf = -VGREAT;
                tensor Finvf(tensor::zero);

                if (admissibleDeformation(I + gradDfI[faceI].T(), Jf, Finvf))
                {
                    SfDirectI[faceI] = Jf*(Finvf.T() & SfI[faceI]);
                }
                else
                {
                    // Owner-cell fallback, preserved deliberately for this
                    // first diagnostic version rather than redesigned
                    SfDirectI[faceI] =
                        J_[own[faceI]]*(Finv_[own[faceI]].T() & SfI[faceI]);
                    nFallback++;
                }
            }
        }

        forAll(SfDirect.boundaryField(), patchI)
        {
            const tensorField& gradDfP = gradDf_.boundaryField()[patchI];
            const vectorField& SfP = mesh.Sf().boundaryField()[patchI];
            const scalarField& JP = J_.boundaryField()[patchI];
            const tensorField& FinvP = Finv_.boundaryField()[patchI];
            vectorField& SfDirectP = SfDirect.boundaryFieldRef()[patchI];

            forAll(SfDirectP, faceI)
            {
                scalar Jf = -VGREAT;
                tensor Finvf(tensor::zero);

                if (admissibleDeformation(I + gradDfP[faceI].T(), Jf, Finvf))
                {
                    SfDirectP[faceI] = Jf*(Finvf.T() & SfP[faceI]);
                }
                else
                {
                    SfDirectP[faceI] = JP[faceI]*(FinvP[faceI].T() & SfP[faceI]);
                    nFallback++;
                }
            }
        }

        reduce(nFallback, sumOp<label>());

        if (nFallback > 0)
        {
            // Never silent
            directFaceFallbackCount_ += nFallback;

            Info<< "    directFaceState: " << nFallback
                << " face(s) inadmissible; owner-cell state substituted"
                << " (cumulative " << directFaceFallbackCount_ << ')' << endl;
        }

        // Internal face force from the direct face state.  The stabilisation
        // is added below in exactly the same way as for the default path.
        // Build the force as |Sf_direct| * (n_direct & sigma_f), i.e. in
        // exactly the same (unoriented) form as the interpolated-cell path
        // uses.  This is numerically identical to SfDirect & sigmaFace, but
        // Sf carries OpenFOAM's "oriented" flag, and an oriented field cannot
        // be added to the unoriented stabilisation term below.
        magSfDirectPtr.reset(new surfaceScalarField("magSfDirect", mag(SfDirect)));
        const surfaceVectorField nDirect(SfDirect/magSfDirectPtr());

        // Direct-face MATERIAL traction; the stabilisation is added later, in
        // exactly the same place and form as for the default treatment
        tractionDirectPtr.reset
        (
            new surfaceVectorField("tractionDirect", nDirect & sigmaFace)
        );

        forceDirectFacePtr.reset
        (
            new surfaceVectorField
            (
                "forceDirectFace", magSfDirectPtr()*tractionDirectPtr()
            )
        );

        if (faceForceDiagnostics_)
        {
            // Reference for the comparison: the interpolated-cell force
            // WITHOUT stabilisation, which is what forceDirectFace also
            // excludes at this point
            reportFaceForceDifference
            (
                surfaceVectorField
                (
                    "forceInterpolatedCell", magSfCurrent*traction
                ),
                forceDirectFacePtr()
            );
        }
    }

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    momentumStabilisation().updateVector(D, &gradD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Calculate the force at the faces.  This is the ORIGINAL expression,
    // unchanged: `traction` already carries the stabilisation.
    surfaceVectorField force(magSfCurrent*traction);

    if (useDirectFace)
    {
        // Same construction, but with the direct face state and the direct
        // material traction.  operator== assigns values without re-checking
        // dimensions or OpenFOAM's "oriented" metadata, so `force` keeps the
        // attributes it just acquired from the expression above.
        force ==
            magSfDirectPtr()
           *(
                tractionDirectPtr()
              + impKf_*momentumStabilisation().faceVector()
            );
    }

    // Enforce traction boundary conditions (unchanged in both treatments)
    enforceTractionBoundaries(force, D, nCurrent, magSfCurrent);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
        // fvc::div(magSfCurrent*traction)
        fvc::div(force)
      + rho()
       *(
            g() - fvc::d2dt2(D) - dampingCoeff()*fvc::ddt(D)
        )
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    //residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    return 0;
}


label nonLinGeomTotalLagTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the D field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    unpackSolution(x);

    // Refuse to assemble from an inadmissible trial state. unpackSolution()
    // returned early, so Finv_, sigma and p still hold the values from the
    // previous evaluation; assembling here would silently mix them with the
    // current gradD/F_. Hand the rejection to PETSc instead. In a
    // NEWTONLS solve this is unreachable in practice, because the Jacobian
    // is evaluated at an iterate whose residual has already been accepted;
    // it is a guard, not a normal path.
    if (!trialGeometryValid_)
    {
        return foamPetscSnesHelper::SNES_EVAL_DOMAIN_ERROR;
    }

    // Take a non-const reference to D for local use below
    volVectorField& D = const_cast<volVectorField&>(this->D());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above and its BCs corrected; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        {
            // Refresh rAUf only when the mesh or deltaT have changed
            // (the diagonal of the approximate momentum Jacobian is
            // independent of D and p values)
            updateRAUfIfStale();

            fvScalarMatrix approxPressureJ
            (
              - pressureEqnScale_*pressureUnknownScale_*fvm::Sp(rKappa(), p)
              + pressureEqnScale_*pressureUnknownScale_
               *pressureStabilisation().scalarJacobian(p, &rAUf())
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );

            // Insert D-in-p equation coefficients matching the
            // linearisation of -0.5*(J^2-1)/J about D=0, which to
            // leading order equals -V*div(D).
            //
            // InsertFvmDivU's sign convention: scale=+1 assembles
            // `-V*div(U)`. So we pass `+pressureEqnScale_` to get
            // J_pD = -alpha*V*div(D).
            foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
            (
                p,
                D,
                jac,
                blockSize_ - 1,             // row offset (p row)
                0,                          // column offset (D columns)
                solidModel::twoD() ? 2 : 3, // number of D components
                pressureEqnScale_           // scale (helper returns -V*div with +1)
            );

            // Insert p-in-D term. InsertFvmGrad's updated sign
            // convention: scale=+1 assembles `-V*grad(p)`. Apply the
            // pressure-unknown scale so the column corresponds to the
            // scaled unknown pHat = p/pressureUnknownScale_.
            foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
            (
                p,
                jac,
                0,                          // row offset
                blockSize_ - 1,             // column offset
                solidModel::twoD() ? 2 : 3, // number of D components
                pressureUnknownScale_       // scale (helper returns -V*grad with +1)
            );
        }
    }

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        momentumStabilisation().vectorJacobian(D, &impKf_)
      - rho()*fvm::d2dt2(D)
    );

    // Viscous material stiffness (Change 6) -- PRECONDITIONER ONLY.
    //
    // The residual carries a Kelvin-Voigt viscous stress S_v = eta*Edot with
    // Edot = current*E(F) + frozen history and dE = symm(F^T dF). Its
    // derivative with respect to D therefore carries the factor eta*current,
    // and about F = I the viscous divergence is that of an isotropic rate
    // material with mu_v = eta/2, lambda_v = 0. Following the impK
    // convention already used above (impK = 2*mu + lambda), the matching
    // Laplacian diffusivity is eta*current, which is what the law returns
    // through mechanicalModel::viscousImpKf().
    //
    // At the benchmark's dt = 5e-5 with BDF2 this is eta*coefft/dt
    // ~ 100*1.5/5e-5 = 3e6 Pa, some four orders of magnitude larger than the
    // impK surrogate (~550 Pa) that was previously the only material term in
    // P_DD. The true operator is viscous-diffusion dominated at this
    // timestep, so without this term the preconditioner was solving an
    // essentially different problem from the one MFFD applies.
    //
    // The sign is positive, matching momentumStabilisation()'s Laplacian:
    // approxJ carries dR/dD and R contains +div(viscous stress).
    //
    // This is assembled into the preconditioner ONLY. formResidual() does
    // not consult it, so the converged solution is unchanged. Note in
    // particular that rAUf() is deliberately NOT rebuilt with this
    // stiffness: rAUf enters the pressure stabilisation in the residual
    // (see formResidual), so altering it would change the physics.
    if (const surfaceScalarField* viscKf = viscousImpKfIfAny())
    {
        approxJ += fvm::laplacian(*viscKf, D, "laplacian(DD,D)");
    }

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D);
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    // Convert fvMatrix matrix to PETSc matrix
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        approxJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
    );

    // Support spring/dashpot stiffness (Change 5B).
    //
    // enforceTractionBoundaries() replaces the boundary face force by
    //
    //     forceP = (traction - nCurrent*pressure)*magSf
    //
    // and fvc::div() adds a boundary face force to its OWNER cell, so after
    // residual *= V the owner cell receives exactly +forceP. For the Arostica
    // supports the traction is a function of the OWNER-CELL state (Change 1):
    //
    //     epicardium  t = -alpha (D_c.N) N - beta (Ddot_c.N) N
    //     base        t = -alpha  D_c      - beta  Ddot_c
    //
    // With Ddot_c = c_time*D_c + history and the history frozen during the
    // current Newton linearisation, the derivative with respect to the owner
    // cell's own displacement unknown is
    //
    //     dR_c/dD_c = -kEff*|Sf|*(N outer N)   (normal-only support)
    //                 -kEff*|Sf|*I             (vector support)
    //
    //     kEff = alpha + beta*c_time
    //
    // This is a pure cell-diagonal block: each support face couples the owner
    // cell to itself only, so no off-diagonal entries arise and contributions
    // from several support faces on one cell simply accumulate. The sign is
    // negative, matching the inertia term above, because approxJ carries
    // dR/dD rather than its negation.
    //
    // The patch supplies kEff and the normal-only flag through
    // solidTractionFvPatchVectorField::supportTangent(), so the time
    // coefficient is derived by the same class that builds the residual and no
    // case-specific boundary type is referenced here.
    {
        const label nCmpt = solidModel::twoD() ? 2 : 3;
        const label bs = blockSize_;

        List<scalar> values(bs*bs, 0.0);

        forAll(D.boundaryField(), patchI)
        {
            const fvPatchVectorField& pD = D.boundaryField()[patchI];

            if (!isA<solidTractionFvPatchVectorField>(pD))
            {
                continue;
            }

            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>(pD);

            scalarField kEff;
            bool normalOnly = false;

            if (!tracPatch.supportTangent(kEff, normalOnly))
            {
                continue;
            }

            // enforceTractionBoundaries() uses the reference area when
            // useUndeformedArea() is set, and the current area otherwise
            scalarField areaP(mesh().boundary()[patchI].magSf());

            if (!tracPatch.useUndeformedArea())
            {
                const surfaceVectorField SfCur
                (
                    fvc::interpolate(J_*Finv_.T()) & mesh().Sf()
                );
                areaP = mag(SfCur.boundaryField()[patchI]);
            }

            const vectorField N(mesh().boundary()[patchI].nf());
            const labelUList& faceCells =
                mesh().boundary()[patchI].faceCells();

            forAll(faceCells, faceI)
            {
                const scalar w = kEff[faceI]*areaP[faceI];
                const vector& n = N[faceI];

                for (label i = 0; i < nCmpt; ++i)
                {
                    for (label j = 0; j < nCmpt; ++j)
                    {
                        values[i*bs + j] =
                            normalOnly
                          ? -w*n[i]*n[j]
                          : (i == j ? -w : 0.0);
                    }
                }

                const label globalRow =
                    foamPetscSnesHelper::globalCells().toGlobal
                    (
                        faceCells[faceI]
                    );

                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalRow, 1, &globalRow,
                        values.cdata(), ADD_VALUES
                    )
                );
            }
        }
    }

    return 0;
}

#endif // USE_PETSC

void nonLinGeomTotalLagTotalDispSolid::updateTotalFields()
{
    if (solveFailed_)
    {
        // The last nonlinear solve did not converge and its result was not
        // accepted, so the accumulated history (in particular the mechanical
        // law's old-time state) must not be advanced as though it had been
        WarningInFunction
            << "Not advancing accumulated fields: the nonlinear solve for time "
            << mesh().time().timeName() << " did not converge" << endl;

        return;
    }

    solidModel::updateTotalFields();
}


bool nonLinGeomTotalLagTotalDispSolid::solutionFailed() const
{
    return solveFailed_;
}


bool nonLinGeomTotalLagTotalDispSolid::reportsJacobianDomainError() const
{
    // formJacobian() can return SNES_EVAL_DOMAIN_ERROR, so PETSc's
    // Jacobian-domain check must be enabled for this model
    return true;
}


void nonLinGeomTotalLagTotalDispSolid::end()
{
    // Lifecycle evidence: the face kinematics are rebuilt only where something
    // consumes them (accepted unpacks, predictor, segregated iterations,
    // construction/restart), not on the MFFD residual and Jacobian probes that
    // dominate the iteration count
    Info<< nl
        << "Face kinematics: " << faceRefreshCount_ << " refresh(es), "
        << faceRefreshSkipCount_
        << " residual/Jacobian probe(s) skipped, "
        << faceSubstitutionCount_ << " admissibility substitution(s)" << endl;

    Info<< "    faceForceTreatment = " << faceForceTreatment_
        << ", directFaceState inadmissible-face fallbacks: "
        << directFaceFallbackCount_ << endl;

    solidModel::end();
}


tmp<vectorField> nonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
