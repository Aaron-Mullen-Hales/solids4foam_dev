/*---------------------------------------------------------------------------* \
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
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"


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

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    
    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        sigma() = dev(sigma()) - p()*I;
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
#ifdef OPENFOAM_NOT_EXTEND
        vectorField& forceP = force.boundaryFieldRef()[patchI];
#else
        vectorField& forceP = force.boundaryField()[patchI];
#endif

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

	    //       traction.boundaryFieldRef()[patchI] =
	    //  tracPatch.traction() - nPatch*tracPatch.pressure();
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
        // To-do: add a stabilisation traction function to momentumStabilisation
        const scalar scaleFactor =
            readScalar(stabilisation().dict().lookup("scaleFactor"));
        const surfaceTensorField gradDf(fvc::interpolate(gradD()));
        traction += scaleFactor*impKf_*(fvc::snGrad(D()) - (n & gradDf));

        // Calculate the force at the faces
        surfaceVectorField force(magSfCurrent*traction);

        // Enforce traction boundary conditions
        enforceTractionBoundaries(force, D(), nCurrent, magSfCurrent);

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(force)
          + rho()*g()
        );

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

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

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

        // Use the segregated solver as a predictor
        //evolveImplicitSegregated();

        // Map the D field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            D().primitiveFieldRef(),
#else
            D().internalField(),
#endif
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD()
          ? makeList<label>({0,1})
          : makeList<label>({0,1,2})
        );

        if (solvePressure())
        {
            // Map the p field to the SNES solution vector
            foamPetscSnesHelper::InsertFieldComponents<scalar>
            (
                p(),
                foamPetscSnesHelper::solution(),
                blockSize_ -1 // Location of first component
            );
        }
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    vectorField& DI = D();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    D().correctBoundaryConditions();

    if (solvePressure())
    {
        // Map the PETSc solution to the p field
        // p is located in the last ("blockSize - 1") component
        scalarField& pI = p();
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            foamPetscSnesHelper::solution(),
            pI,
            blockSize_ - 1 // Location of p component
        );

        p().correctBoundaryConditions();

        // Update dpdt
        autoPtrRef(dpdtPtr_) = fvc::ddt(p());
    }

    // Update gradient of displacement
    mechanical().grad(D(), gradD());
    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    // Update total deformation gradient
    //F_ = I + gradD().T(); // or your preferred linearized version

    // Update Jacobian to match -tr(gradD) exactly
    //J_ = 1.0 + tr(gradD()); // linearized Jacobian
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


void nonLinGeomTotalLagTotalDispSolid::makePDiffusivity() const
{
    if (pDiffusivityPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const scalar pressureSmoothingCoeff
    (
        readScalar(solidModelDict().lookup("pressureSmoothingCoeff"))
    );

    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D())
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D());
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    pDiffusivityPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "pDiffusivity",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -pressureSmoothingCoeff*impKf_/fvc::interpolate(approxJ.A())
        )
    );
}


const surfaceScalarField& nonLinGeomTotalLagTotalDispSolid::pDiffusivity() const
{
    if (pDiffusivityPtr_.empty())
    {
        makePDiffusivity();
    }

    return autoPtrRef(pDiffusivityPtr_);
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
        inv(F_)
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
        det(F_)
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
    pDiffusivityPtr_(),
    dpdtPtr_(),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
     ),
    ds_
    (
        IOobject
	(
            "ds",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
	),
        mesh(),
        dimensionedScalar("ds", (dimForce/dimVolume)/dimVelocity, 1.0)
    )
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields
    mechanical().correct(sigma());

    D().correctBoundaryConditions();
    D().storePrevIter();
    mechanical().grad(D(), gradD());
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
	//p();
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
        Finv_ = inv(F_);
        J_ = det(F_);

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
#ifdef OPENFOAM_NOT_EXTEND
                        D().boundaryFieldRef()[patchI]
#else
                        D().boundaryField()[patchI]
#endif
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


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


label nonLinGeomTotalLagTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    static label iterCount = 0;
    const fvMesh& mesh = this->mesh();

    //adding a scaling term to try make it easier for the linear solver in these problems.
    //scalar unitsPressureScale = solidModelDict().lookupOrDefault<scalar>("unitsPressureScale", 1.0);
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

    // Enforce the boundary conditions
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

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    //Finv transpose
    //Finv_T_ = Finv_.T()

    // Jacobian of the deformation gradient
    J_ = det(F_);
    //J_ = I + tr(gradD());
    //J_ =volScalarField("J_", mesh, 1.0) + tr(gradD());

    // volScalarField J_linear
    //   (
    //    IOobject("J_linear", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    //    mesh,
    //    dimensionedScalar("one", dimless, 1.0)
    //    );

    // J_linear += tr(gradD());
    // J_ = J_linear;   // assign to your main J_

    // J_ = 1.0 + tr(gradD());
    // forAll(J_, cellI)
    //   {
    // 	const tensor& Fc = gradD()[cellI];
    // 	// Force exact equality with -tr(gradD())
    // 	J_[cellI] = 1.0 - (-Fc.xx() - Fc.yy() - Fc.zz());
    //   }

    // // Boundary patches
    // forAll(J_.boundaryField(), patchI)
    //   {
    // 	const tensorField& gradDb = gradD().boundaryField()[patchI];
    // 	scalarField& Jb = J_.boundaryFieldRef()[patchI];

    // 	forAll(gradDb, faceI)
    // 	  {
    // 	    const tensor& Fc = gradDb[faceI];
    // 	    Jb[faceI] = 1.0 - (-Fc.xx() - Fc.yy() - Fc.zz());
    // 	  }
    //   }

    // Info<< "max(detF_manual - 1): " << -(max(J_) - 1) << nl
    // 	<< "min(detF_manual - 1): " <<  -(min(J_) - 1.0 )<< nl
    //  << "max(-tr(gradD))    : " << max(-tr(gradD())) << nl
    //  << "min(-tr(gradD))    : " << min(-tr(gradD())) << endl;

    //    const volTensorField& F = F_;       // deformation gradient


    // NEED THIS TO PRINT OUT MANUAL DETERMINANT OF F (J_ = det(F_))
    // const volTensorField& F = F_;
    // // Manual determinant computation (per cell)
    // volScalarField detF_manual(
    // 			       IOobject("detF_manual", mesh.time().timeName(), mesh),
    // 			       mesh,
    // 			       dimensionedScalar("zero", dimless, 0.0)
    // 			       );


    // // --- Internal field
    // forAll(detF_manual, cellI)
    //   {
    // 	const tensor& Fc = F[cellI];

    // 	// Force all intermediates to double
    // 	const double f11 = static_cast<double>(Fc.xx());
    // 	const double f12 = static_cast<double>(Fc.xy());
    // 	const double f13 = static_cast<double>(Fc.xz());
    // 	const double f21 = static_cast<double>(Fc.yx());
    // 	const double f22 = static_cast<double>(Fc.yy());
    // 	const double f23 = static_cast<double>(Fc.yz());
    // 	const double f31 = static_cast<double>(Fc.zx());
    // 	const double f32 = static_cast<double>(Fc.zy());
    // 	const double f33 = static_cast<double>(Fc.zz());

    // 	// Compute det in double precision
    // 	const double detF_d =
    // 	  f11*(f22*f33 - f23*f32)
    // 	  - f12*(f21*f33 - f23*f31)
    // 	  + f13*(f21*f32 - f22*f31);

    // 	// Store back as scalar
    // 	detF_manual[cellI] = static_cast<scalar>(detF_d);
    //   }

    // // --- Boundary field
    // forAll(detF_manual.boundaryField(), patchI)
    //   {
    // 	const tensorField& Fb = F.boundaryField()[patchI];
    // 	scalarField& detFb = detF_manual.boundaryFieldRef()[patchI];

    // 	forAll(Fb, faceI)
    // 	  {
    // 	    const tensor& Fc = Fb[faceI];

    // 	    // Force all intermediates to double
    // 	    const double f11 = static_cast<double>(Fc.xx());
    // 	    const double f12 = static_cast<double>(Fc.xy());
    // 	    const double f13 = static_cast<double>(Fc.xz());
    // 	    const double f21 = static_cast<double>(Fc.yx());
    // 	    const double f22 = static_cast<double>(Fc.yy());
    // 	    const double f23 = static_cast<double>(Fc.yz());
    // 	    const double f31 = static_cast<double>(Fc.zx());
    // 	    const double f32 = static_cast<double>(Fc.zy());
    // 	    const double f33 = static_cast<double>(Fc.zz());

    // 	    // Compute det in double precision
    // 	    const double detF_d =
    // 	      f11*(f22*f33 - f23*f32)
    // 	      - f12*(f21*f33 - f23*f31)
    // 	      + f13*(f21*f32 - f22*f31);

    // 	    // Store back as scalar
    // 	    detFb[faceI] = static_cast<scalar>(detF_d);
    // 	  }
    //   }



    // // --- Internal field
    // forAll(detF_manual, cellI)
    //   {
    // 	const tensor& Fc = F[cellI];
    // 	// Linearized det(F) ~ 1 + tr(F - I) = 1 + tr(gradD)
    // 	double trGradD = static_cast<double>(Fc.xx() + Fc.yy() + Fc.zz() - 3.0);
    // 	detF_manual[cellI] = scalar(1.0 + trGradD);
    //   }

    // // --- Boundary patches
    // forAll(detF_manual.boundaryField(), patchI)
    //   {
    // 	const tensorField& Fb = F.boundaryField()[patchI];
    // 	scalarField& detFb = detF_manual.boundaryFieldRef()[patchI];

    // 	forAll(Fb, faceI)
    // 	  {
    // 	    const tensor& Fc = Fb[faceI];
    // 	    double trGradD = static_cast<double>(Fc.xx() + Fc.yy() + Fc.zz() - 3.0);
    // 	    detFb[faceI] = scalar(1.0 + trGradD);
    // 	  }
    //   }




    // detF_manual.correctBoundaryConditions();

    // Info<< "max(detF_manual - 1): " << max(detF_manual) - 1.0 << nl
    // 	<< "min(detF_manual - 1): " <<  min(detF_manual) - 1.0 << nl
    // 	<< "max(-tr(gradD))    : " << max(-tr(gradD())) << nl
    // 	<< "min(-tr(gradD))    : " << min(-tr(gradD())) << endl;

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

	//dimensionedScalar mu_("mu", [1 -1 -2 0 0 0 0] );
	const dimensionedScalar mu_("mu", solidModelDict());
	const dimensionedScalar lambda_("lambda", solidModelDict());
        // Enforce the boundary conditions
        p.correctBoundaryConditions();
        // Replace the pressure component of stress
	//sigma() = 2*symm(mu_ *gradD()) + lambda_*tr(gradD())*I;
	//sigma() = 2*symm(mu_ *gradD()) + lambda_*tr(gradD())*I;

	//Normal Way -- non Guccione
	sigma() = dev(sigma()) - p*I;

	//Attempting Guccione- which might need all of sigma:
	//sigma() = sigma() - p*I;
        //sigma() = dev(sigma()) - (p * unitsPressureScale)*I;
    }

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    // const surfaceVectorField SfCurrent
    // (
    //     mesh.Sf()
    // );
    // const surfaceScalarField magSfCurrent(mag(SfCurrent));
    // const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);
    //    const surfaceVectorField nCurrent(mesh.Sf()/mesh.magSf());

    // Traction vectors at the faces
    //surfaceVectorField traction(n & fvc::interpolate(sigma()));
    surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

    //fvc::div(J_*Finv_ & sigma(), "div(sigma)");

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Calculate the force at the faces
    //surfaceVectorField force(traction);
    surfaceVectorField force(magSfCurrent*traction);

    // Enforce traction boundary conditions
    enforceTractionBoundaries(force, D, nCurrent, magSfCurrent);
    //enforceTractionBoundaries(traction, D, nCurrent, magSfCurrent);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
     //fvc::div(mesh.magSf()*traction)
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

    if (solvePressure())
    {
      //p.correctBoundaryConditions();
        volScalarField& p = const_cast<volScalarField&>(this->p());

	//p.correctBoundaryConditions();
        // Divided by bulkModulus form
        const volScalarField kappa("kappa", mechanical().bulkModulus());
        const surfaceScalarField kappaf(fvc::interpolate(kappa));

	const dimensionedScalar omega("omega", solidModelDict());
        const dimensionedScalar omegaTau("omegaTau", solidModelDict());
        const word stabilisationType =
          solidModelDict().lookupOrDefault<word>("stabilisationType", "rhieChow");

	scalarField pressureResidual(mesh.nCells(), 0.0);
	if (stabilisationType == "rhieChow")
	  {
	    // Create the diffusivity field properly

	    //const volScalarField Dp(pDiffusivity()/impK_);
	    pressureResidual =
	    (
	     - p/kappa
	     + fvc::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
	     - fvc::div
	     (
	      (pDiffusivity()/kappaf)*mesh.Sf()
	      //(pDiffusivity()/kappaf)*SfCurrent
	       & fvc::interpolate(fvc::grad(p))
	     )
	  //-0.5*kappa*(pow(J_, 2.0) - 1.0)/(J_) //*(1e-6))
	     //- tr(gradD())

	     -0.5*(pow(J_, 2.0) - 1.0)/(J_)
	     //-0.5*unitsPressureScale*(pow(J_, 2.0) - 1.0)/(J_)
	     );

	}
	else if (stabilisationType == "HigherOrder")
	  {
	    surfaceScalarField beta
	      (
	       "beta",
	       omega/(impKf_*sqr(mesh.deltaCoeffs()))
	       );

	    volScalarField sixthOrderStabilisation
	      (
	       "sixthOrderStabilisation",
	       fvc::laplacian
	       (
		fvc::laplacian
		(
		 fvc::laplacian(beta, p)   // beta is face-based: OK here
		 )
		)
	       );
	    //Info<< "dims(stab)    = " << sixthOrderStabilisation.dimensions() << nl;
	    //Info<< "dims(p)    = " << p.dimensions() << nl;
	    //Info<< "dims(kappa)    = " << kappa.dimensions() << nl;
	    //Info<< "max|stab|    = " << gMax(sixthOrderStabilisation) << nl;
	    //Info<< "max|p| = " << gMax(p) << nl;
	    //Info<< "max|kappa| = " << gMax(kappa) << nl;


	    //const volScalarField Dp(pDiffusivity()/impK_);
	    pressureResidual =
	    (
	     - p/kappa
	     //+ omega*fvc::laplacian(omega*fvc::laplacian(omega*fvc::laplacian(omega*1.0/(impKf_*sqr(mesh.deltaCoeffs())), p, "laplacian(Dp,p)")))
	     +fvc::laplacian
	       (
		fvc::laplacian
		(
		 fvc::laplacian(beta, p)   // beta is face-based: OK here
		 )
		)
	     //sixthOrderStabilisation
	  //-0.5*kappa*(pow(J_, 2.0) - 1.0)/(J_) //*(1e-6))
	     //- tr(gradD())

	     -0.5*(pow(J_, 2.0) - 1.0)/(J_)
	     //-0.5*unitsPressureScale*(pow(J_, 2.0) - 1.0)/(J_)
	     );
	    //Info<< "dims(p/kappa) = " << (p/kappa).dimensions() << nl;

	}
	else if (stabilisationType == "oosterlee")
	{
	  if (iterCount < 100000000){
	    // Oosterlee formulation
	    //p.correctBoundaryConditions();

	    const word pressureConfig = solidModelDict().lookupOrDefault<word>("pressureConfig", "reference");

	    if (pressureConfig == "deformed") {
	      // Deformed configuration pressure
	      //const surfaceScalarField kappafCurrent = fvc::interpolate(kappa*J_);
	      // const surfaceScalarField kappafCurrent(fvc::interpolate(kappa*J_));
	      // const surfaceVectorField gradpCurrent(fvc::interpolate(Finv_ & fvc::grad(p)));
	      // const surfaceScalarField DpDeformed(
	      // 					  fvc::interpolate(J_) * omega/(impKf_*sqr(mesh.deltaCoeffs()))
	      // 					  );

	      // Use Nanson's formula to get deformed face areas
	      const surfaceTensorField FinvTf(fvc::interpolate(Finv_.T()));
	      const surfaceScalarField Jf(fvc::interpolate(J_));
	      const surfaceVectorField SfCurrent(Jf * (FinvTf & mesh.Sf()));
	      const surfaceScalarField magSfCurrent(mag(SfCurrent));

	      // Compute deformed cell distances (approximate) by scaling
	      // the mesh deltacoeffs by the ratio of deformed to reference
	      const surfaceScalarField deltaCurrent(
						    (magSfCurrent/mesh.magSf()) * mesh.deltaCoeffs()
						    );

	      // Stabilization based on deformed geometry
	      const surfaceScalarField DpDeformed(
						  omega/(impKf_ * sqr(deltaCurrent))
						  );
	      //const surfaceScalarField DpDeformed = fvc::interpolate(J_ * Finv_ & Finv_.T()) * omega/(impKf_*sqr(mesh.deltaCoeffs()));
	      //const surfaceVectorField gradpCurrent = fvc::interpolate(Finv_ & fvc::grad(p));

	      pressureResidual =
		(
		 - p/kappa
		 + fvc::laplacian(DpDeformed, p, "laplacian(DpDeformed,p)")
		 //+ fvc::div(omega * kappafCurrent * (mesh.Sf()/mesh.magSf()) & gradpCurrent)
		 -0.5*(pow(J_, 2.0) - 1.0)/(J_)
		 //- (J_ - 1.0)
		 );
	    }
	    else {
	      // Reference configuration (your current with simplified volume change)
	      pressureResidual =
		(
		 - p/kappa
		 + fvc::laplacian(omega/(impKf_*sqr(mesh.deltaCoeffs())), p, "laplacian(Dp,p)")
		 -0.5*(pow(J_, 2.0) - 1.0)/(J_)
		 //- (J_ - 1.0)  // Use this instead of complex form for better convergence
		 );
	    }
	  }


	    //old previous approach:


	    // pressureResidual =
	    //   (
	    //    //p.correctBoundaryConditions();
	    //    //lambdaf
	    //    //- p/lambda +
	    //    - p/kappa +
	    //    fvc::laplacian(omega/(impKf_*sqr(mesh.deltaCoeffs())), p, "laplacian(Dp,p)")
	    //    //              - fvc::div(((omega/sqr(mesh.deltaCoeffs())))*mesh.Sf() & fvc::interpolate(fvc::grad(p)))
	    //    // Different stabilisation term(s) here
	    //    //+ 1 - (1 + tr(gradD()))
	    //    //- tr(gradD()) - 0.5*(sqr(tr(gradD())) - tr( gradD() & gradD() ) )
	    //    -0.5*(pow(J_, 2.0) - 1.0)/(J_)
	    //    //trying taylor series approximation
	    //    //2nd order
	    //    //-(J_ - 1 - 0.5*sqr(J_ -1))
	    //    // 3rd-order
	    //    //-(J_-1 - 0.5*sqr(J_-1) + (1.0/3.0)*pow(J_-1,3))
	    //    // 4th-order
	    //    //-(J_-1 - 0.5*sqr(J_-1) + (1.0/3.0)*pow(J_-1,3) - 0.25*pow(J_-1,4))
	    //    //- log(J_)
	    //    //-0.5*unitsPressureScale*(pow(J_, 2.0) - 1.0)/(J_)
	    //    );

	  else
	    {
	      //p.correctBoundaryConditions();

	      pressureResidual =
              (
	      //lambdaf
	       //- p/lambda +
	       - p/kappa +
	       fvc::laplacian(omega/(impKf_*sqr(mesh.deltaCoeffs())), p, "laplacian(Dp,p)")
	       //              - fvc::div(((omega/sqr(mesh.deltaCoeffs())))*mesh.Sf() & fvc::interpolate(fvc::grad(p)))
	       // Different stabilisation term(s) here
	       //- tr(gradD())
	       //+ 1 - (1 + tr(gradD()))
	       //-(J_ -1)
	       //-(J_ - 1 - 0.5*sqr(J_ -1))
	       //- tr(gradD()) - 0.5*(sqr(tr(gradD())) - tr( gradD() & gradD() ) )
	       -0.5*(pow(J_, 2.0) - 1.0)/(J_)
	       //- log(J_)
		 );
	    }
	}
	else if (stabilisationType == "Mixed")
	  {
	    pressureResidual =
	      (
	       -p/kappa +
	       0.5*(fvc::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
             - fvc::div
             (
              (pDiffusivity()/kappaf)*mesh.Sf()
              //(pDiffusivity()/kappaf)*SfCurrent
               & fvc::interpolate(fvc::grad(p))
             )
		    ) + 0.5*(fvc::laplacian(omega/(impKf_*sqr(mesh.deltaCoeffs())), p, "laplacian(Dp,p)"))
	       -0.5*(pow(J_, 2.0) -1.0)/(J_)
              );
	  }
	else
        {
            FatalError
                << "stabilisationType unknown: " << stabilisationType
                << exit(FatalError);
        }
	//Info << "impkf value: " << impKf_ << endl;

        // scalarField pressureResidual
        // (
        //   - p/kappa
        //   + fvc::laplacian
        //     (
        //         omega/sqr(mesh.deltaCoeffs()/impKf_), p, "laplacian(Dp,p)"
        //     )
        //   - 0.5*(pow(J_, 2.0) - 1.0)/J_
        // );

        // Make residual extensive
        pressureResidual *= mesh.V();
	//scalar matrixScale = solidModelDict().lookupOrDefault<scalar>("matrixScale", 1.0);
	//pressureResidual *= matrixScale;

	//adding a switch to look at scale of residuals
	Switch printResiduals =
	solidModelDict().lookupOrDefault<Switch>("printResiduals", false);
	if (printResiduals){
	  Info << "Displacement Residual:\n" << average(mag(residual)) << endl;
	Info << "Pressure Residual:\n" << average(mag(pressureResidual)) << endl;
	// Copy the pressureResidual into the f field as the final equation
	}


        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
	 );
    }

    // if (solvePressure())
    //   {
    // 	const volVectorField& Dfield = this->D();
    // 	const volScalarField& pfield = this->p();

    // 	// Explicitly dereference the tmp<> returned by mag()
    // 	const volScalarField& magD = mag(Dfield)();  // note the extra ()
    // 	const volScalarField& magP = mag(pfield)();  // same here

    // 	// Compute global maxima
    // 	scalar Dmax = gMax(magD.internalField());
    // 	scalar pmax = gMax(magP.internalField());

    // 	if (Dmax > SMALL)
    // 	  {
    // 	    Info<< "p:D ratio = " << pmax / Dmax
    // 		<< " (pmax = " << pmax
    // 		<< ", Dmax = " << Dmax << ")" << endl;
    // 	  }

    // 	else
    // 	  {
    // 	    Info<< "Warning: Dmax ~ 0, cannot compute p:D ratio" << endl;
    // 	  }
    //   }

    //Info << "printing itercount" << iterCount << endl;
    //iterCount++;


    return 0;
}


label nonLinGeomTotalLagTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{

    //This gives freedom to switch between pa and mpa ect by setting the scaling in the dict
    //scalar unitsPressureScale = solidModelDict().lookupOrDefault<scalar>("unitsPressureScale", 1.0);

    //note we only wish to build the preconditioner once.
    if(preconditionerBuilt_)
    {
      Info<< "Skipping pre conditioner build as was done."<< nl <<endl;
      return 0;
    }
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

    // Enforce the boundary conditions
    D.correctBoundaryConditions();

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();
    }

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D, "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D)
    );

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

    if (solvePressure())
    {
        const volScalarField& p = this->p();

        const volScalarField kappa("kappa", mechanical().bulkModulus());
        //const volScalarField rKappa(1.0/mechanical().bulkModulus());
        const volScalarField rKappa(1.0/kappa);
	//Info << "printing rKappa" << rKappa <<endl;
        const surfaceScalarField kappaf(fvc::interpolate(kappa));
        const dimensionedScalar omega(solidModelDict().lookup("omega"));
	//scalar matrixScale = solidModelDict().lookupOrDefault<scalar>("matrixScale", 1.0);
        const word stabilisationType =
          solidModelDict().lookupOrDefault<word>("stabilisationType", "rhieChow");
        const dimensionedScalar omegaTau("omegaTau", solidModelDict());

	// LAMBDA CASE DID NOT WORK CAN PROBABLY DELETE THIS STUFF WHEN CLEANING UP.
	// tmp<volScalarField> lambdaField = mechanical().bulkModulus() - (2.0/3.0)*mechanical().shearModulus();
	// const volScalarField& lambda = lambdaField(); // extract the underlying volScalarField

	// const volScalarField rLambda(1.0/kappa);
	// // Step 2: interpolate → tmp<surfaceScalarField>
	// tmp<surfaceScalarField> lambdafTmp = fvc::interpolate(lambda);

	// // Step 3: extract the underlying surfaceScalarField
	// const surfaceScalarField& lambdaf = lambdafTmp();
        {
	  // Calculate pressure equation matrix
	  //const dimensionedScalar one("one", dimless, 1);
	  if(stabilisationType == "rhieChow")
	    {
	      fvScalarMatrix approxPressureJ
		(
		  - fvm::Sp(rKappa, p)
		 //- fvm::Sp(rKappa , p)//*matrixScale
		 + fvm::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
		 //+ fvm::laplacian(
		 //                omega/(impKf_*sqr(mesh().deltaCoeffs())),
		 //                p,
		 //                "jacobian-laplacian(rAU,p)"
		 //                )//*matrixScale


		 );
	      //approxPressureJ *= matrixScale;

	      // Insert the pressure equation
	      foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
		(
		 approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
		 );
	    }
	  else if (stabilisationType =="Mixed")
	    {
	      surfaceScalarField Dp(omega/(impKf_*sqr(mesh().deltaCoeffs())));
              fvScalarMatrix approxPressureJ
		(
		  - fvm::Sp(rKappa, p)
		 //- fvm::Sp(rKappa , p)//*matrixScale
		 + 0.5*fvm::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
		 +0.5*fvm::laplacian(Dp, p, "jacobian-laplacian(rAU,p)")


		 );
	      foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
		(
		 approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
		 );
	    }
	  else if (stabilisationType =="HigherOrder")
	    {
	      surfaceScalarField Dp((omegaTau/(impKf_*sqr(mesh().deltaCoeffs()))));
              fvScalarMatrix approxPressureJ
		(
		  - fvm::Sp(rKappa, p)
		  + fvm::laplacian(Dp, p, "jacobian-laplacian(rAU,p)")

		  
		 );
	      foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
		(
		 approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
		 );
	    }
	  else if(stabilisationType == "oosterlee")
	    {
	      const word pressureConfig = solidModelDict().lookupOrDefault<word>("pressureConfig", "reference");

	      // Calculate Dp based on configuration
	      //surfaceScalarField Dp = omega/(impKf_*sqr(mesh().deltaCoeffs()));

	      surfaceScalarField Dp(omega/(impKf_*sqr(mesh().deltaCoeffs())));
	      if (pressureConfig == "deformed")
		{
		  // Use the SAME Nanson-based approach as the residual
		  const surfaceTensorField FinvTf(fvc::interpolate(Finv_.T()));
		  const surfaceScalarField Jf(fvc::interpolate(J_));
		  const surfaceScalarField magSfCurrent(mag(Jf * (FinvTf & mesh().Sf())));  // mesh().Sf()

		  // Compute length scale ratio using same formula as residual
		  const surfaceScalarField lengthScaleRatio(
							    pow(magSfCurrent/mesh().magSf(), 1.0/(mesh().nGeometricD()-1))  // mesh().magSf()
							    );

		  const surfaceScalarField deltaCurrent(lengthScaleRatio / mesh().deltaCoeffs());  // mesh().deltaCoeffs()

		  // Update Dp with deformed geometry
		  Dp = omega/(impKf_ * sqr(1.0/deltaCurrent));
		}

	      fvScalarMatrix approxPressureJ
		(
		 - fvm::Sp(rKappa, p)
		 + fvm::laplacian(Dp, p, "jacobian-laplacian(rAU,p)")
		 );
	      //surfaceScalarField omegaCoeff = fvc::interpolate(omega / mesh().deltaCoeffs());

	      // Insert the pressure equation
	      //approxPressureJ *= matrixScale;
	      foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
		(
		 approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
		 );
	    }
	  else
	    {
	      FatalError
		<< "stabilisationType unknown: " << stabilisationType
		<< exit(FatalError);
	    }


            // fvScalarMatrix approxPressureJ
            // (
            //   - fvm::Sp(rKappa, p)
            //   + fvm::laplacian
            //     (
            //         omega/sqr(mesh().deltaCoeffs())/impKf_,
            //         p,
            //         "jacobian-laplacian(rAU,p)"
            //     )
            // );

            // // Insert the pressure equation
            // foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            // (
            //     approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            // );
        }



	//volScalarField& J = const_cast<volScalarField&>(this->J()); // or wherever your J is

	//f(J) = 0.5*(J^2 - 1)/J
	//df/dJ = (J^2 + 1)/(2*J^2)
	//Chain rule weight = df/dJ * J

	// volScalarField weight
	//   (
	//    IOobject("volumetricWeight", mesh().time().timeName(), mesh()),
	//    (sqr(J_) + scalar(1.0)) / (scalar(2.0) * J_)
	//    );




	// foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
	//   (
	//    p,
	//    D,
	//    jac,
	//    blockSize_ - 1,          // row offset for pressure
	//    0,                       // column offset for displacement
	//    solidModel::twoD() ? 2 : 3,
	//    false,                   // flipSign (keep same as before)
	//    1.0,
	//    &weight                 // pass pointer to weight // scaleFactor
	//    );



        // Insert D-in-p equation coeffs coming from tr(grad(D)) == div(D)
	// the below is the block of the jacobian corresponding to Jpu
        foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
        (
            p,
            D,
            jac,
            blockSize_ - 1,            // row offset
            0,                         // column offset
            solidModel::twoD() ? 2 : 3//, // number of scalar components of D
            //matrixScale*unitsPressureScale //adding a scaling option
        );
	// {
	//   // Build the correction term for dRp/dD (∂R_p/∂D)
	//   // This represents −(1 + 1/sqr(J_)) * Finv_T · grad(ΔD)
	//   // Equivalent to a Laplacian-like correction weighted by deformation

	//   //volTensorField gradDivCorr = (1.0 + 1.0/sqr(J_)) * (Finv_.T());
	//   volTensorField gradDivCorr = ((1.0 + 1.0/sqr(J_)) * (Finv_.T())).ref();

	//   fvVectorMatrix dRp_dDcorr(
	// 			    fvm::laplacian(gradDivCorr, D, "dRp_dDcorr")

	// 			    );
	//   label nDim = solidModel::twoD() ? 2 : 3;
	//   for (label cellI = 0; cellI < D.internalField().size(); ++cellI)
	//     {
	//       scalar Jinv2 = 1.0 + 1.0/sqr(J_[cellI]); // (1 + 1/J^2) factor
	//       for (label i = 0; i < nDim; ++i)
	// 	{
	// 	  // Compute the nonlinear coefficient for this cell
	// 	  scalar coeff = -Jinv2 * Finv_.T()[cellI][i][i]; // diagonal approximation

	// 	  // Insert into PETSc Jacobian at block (pressure row, displacement column)
	// 	  foamPetscSnesHelper::AddScalarToPETScBlock(
	// 						     jac,
	// 						     cellI,              // row = pressure cell index
	// 						     i,                  // column = displacement component
	// 						     blockSize_-1,        // row block offset (pressure)
	// 						     coeff
	// 						     );
	// 	}
	//     }

	  // for (label i = 0; i < nDim; ++i)
	  //   {
	  //     // Extract the i-th component of the tensor as a scalar coefficient
	  //     volScalarField diffusivity = gradDivCorr.component(i).ref();

	  //     // Extract i-th component of D as a scalar field
	  //     volScalarField D_i = D.component(i)();

	  //     // Construct scalar fvMatrix
	  //     fvMatrix<scalar> dRp_dDcorr_i(
	  // 				    fvm::laplacian(diffusivity, D_i, "dRp_dDcorr_" + Foam::name(i))
	  // 				    );

	  //     // Insert into PETSc
	  //     foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>(
	  // 								 dRp_dDcorr_i, jac, blockSize_-1, i, 1
	  // 								 );
	  //   }

	  // Insert this into the same ∂Rp/∂D block (row = pressure, column = displacement)
	  //	  foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>(
	  //							     dRp_dDcorr, jac, blockSize_-1, 0, 1
	  //							     );
	//      }

        // Insert p-in-D term
        // Insert "-grad(p)" (equivalent to "-div(p*I)") into the D equation
        foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
        (
            p,
            jac,
            0,                         // row offset
            blockSize_ - 1,            // column offset
            solidModel::twoD() ? 2 : 3//, // number of scalar equations to insert
	    //unitsPressureScale
        );
    }
    preconditionerBuilt_ = false;

    return 0;
}

#endif // USE_PETSC

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
