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
    nonLinGeomTotalLagTotalDispSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Philip Cardiff, UCD.

SourceFiles
    nonLinGeomTotalLagTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomTotalLagTotalDispSolid_H
#define nonLinGeomTotalLagTotalDispSolid_H

#include "solidModel.H"
#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class nonLinGeomTotalLagTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomTotalLagTotalDispSolid
:
    public solidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Jacobian departure from incompressibility
        volScalarField Jminus1_;

        //- Geometric Jacobian from reconstructed point displacement
        volScalarField Jgeom_;

        //- Difference between gradient and geometric Jacobians
        volScalarField JgradMinusJgeom_;

        //- Geometric Jacobian departure from incompressibility
        volScalarField JgeomMinus1_;

        //- Pressure compressibility term in the unscaled pressure residual
        volScalarField pressureCompressibilityTerm_;

        //- Pressure stabilisation term in the unscaled pressure residual
        volScalarField pressureStabilisationTerm_;

        //- Volumetric constraint term in the unscaled pressure residual
        volScalarField volumetricConstraintTerm_;

        //- Sum of unscaled pressure residual diagnostic terms
        volScalarField pressureConstraintResidual_;

        //- Inertia force density in the momentum residual
        volVectorField inertiaForceDensity_;

        //- Damping force density in the momentum residual
        volVectorField dampingForceDensity_;

        //- Use a supplied reference-volume body-force field when enabled
        Switch useBodyForceField_;

        //- Body force density in the momentum residual
        volVectorField bodyForceDensity_;

        //- Physical stress force density excluding momentum stabilisation
        volVectorField physicalStressForceDensity_;

        //- Momentum stabilisation force density
        volVectorField momentumStabForceDensity_;

        //- Mask for cells with Jgrad - 1 above threshold
        volScalarField highJgradMask_;

        //- Mask for cells with Jgrad - 1 below negative threshold
        volScalarField lowJgradMask_;

        //- Mask for cells with Jgeom - 1 above threshold
        volScalarField highJgeomMask_;

        //- Mask for cells with Jgeom - 1 below negative threshold
        volScalarField lowJgeomMask_;

        //- Mask for cells with large Jgrad - Jgeom difference
        volScalarField largeJgradMinusJgeomMask_;

        //- Acceleration
        volVectorField A_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        const surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Reciprocal of the bulk modulus
        //  Used in the pressure equation when solvePressure is enabled
        mutable autoPtr<volScalarField> rKappaPtr_;

        //- Time rate of pressure
        //  This is used by the optional linear predictor
        mutable autoPtr<volScalarField> dpdtPtr_;

        //- Time index when solidModel::rAUf() was last computed
        //  (-1 = never built)
        mutable label rAUfTimeIndex_;

        //- deltaT for which solidModel::rAUf() was last computed.
        //  rAUf depends only on impKf_, rho, mesh and deltaT (not on
        //  the values of D or p), so it is safe to reuse across PETSc
        //  matrix-free residual evaluations within a Newton solve
        mutable scalar rAUfDeltaT_;

        //- Store pressure in the PETSc solution vector as p/pScale
        const Switch scaleMixedPetScFields_;

        //- Scale type used for the PETSc pressure unknown
        const word pressureUnknownScaleType_;

        //- Pressure scale used for pHat = p/pressureUnknownScale_
        scalar pressureUnknownScale_;

        //- User-facing dimensionless scale factor applied to the
        //  pressure-equation residual and Jacobian rows relative to
        //  the momentum equations
        const scalar pressureScaleFactor_;

        //- If true, multiply pressureScaleFactor_ by a representative
        //  2*mu so that the pressure-equation residual has units
        //  comparable to the momentum equation
        const Switch pressureScaleByTwoMu_;

        //- Representative value of 2*mu used as the physical scale of
        //  the pressure equation. Computed once at construction
        scalar twoMuRef_;

        //- Effective scale applied to the pressure residual and the
        //  pressure-row Jacobian entries.
        //  pressureEqnScale_ = pressureScaleFactor_
        //                       * (pressureScaleByTwoMu_ ? twoMuRef_ : 1)
        scalar pressureEqnScale_;

        //- Pressure-row normalisation. "legacy" retains the established
        //  pressureEqnScale_*V residual rows. "volumeRmsForce" uses a
        //  volume-weighted RMS normalisation with a characteristic force
        //  pressureEqnScale_*L^2, without changing the equation root.
        const word pressureRowScaling_;

        //- Report momentum stabilisation diagnostics
        const Switch reportMomentumStabDiagnostics_;

        //- Residual-call interval for momentum stabilisation diagnostics
        const label momentumStabReportInterval_;

        //- Write momentum stabilisation diagnostic fields
        const Switch writeMomentumStabFields_;

        //- Boundary treatment for the momentum stabilisation face force
        const word momentumStabilisationBoundaryTreatment_;

        //- Face traction route: legacy interpolated cell stress or direct
        //  constitutive evaluation from the current face deformation.
        const word faceStressTreatment_;

        //- Optional passive material tangent contribution to the compact
        //  displacement preconditioner. This changes P only, never the
        //  physical residual.
        const Switch preconditionerMaterialTangent_;

        //- Optional Aróstica spring-dashpot tangent contribution to P only.
        const Switch preconditionerBoundaryTangent_;

        //- Optional direct-face current-trial viscous nominal-force tangent.
        const Switch preconditionerViscousTangent_;

        //- Optional direct-face passive nominal-force tangent.
        const Switch preconditionerPassiveNominalTangent_;

        //- Use the least-squares displacement gradient stencil when
        //  assembling the pressure-equation displacement block.
        const Switch preconditionerLeastSquaresPressureCoupling_;

        //- Include the derivative of the interpolated least-squares pressure
        //  gradient in the pressure stabilisation preconditioner block.
        const Switch preconditionerLeastSquaresPressureStabilisation_;

        //- Counter used to throttle momentum stabilisation diagnostics
        mutable label momentumStabDiagnosticCall_;

        //- Report every evaluated PETSc SNES line-search trial
        const Switch reportSnesTrialDiagnostics_;

        //- Return inadmissible PETSc line-search trials as function-domain
        //  errors so a compatible line search can backtrack
        const Switch petscDomainSafeTrials_;

        //- Iteration associated with the current trial counter
        mutable label lastSnesTrialDiagnosticIteration_;

        //- Per-iteration line-search trial counter
        mutable label snesLineSearchTrialCounter_;

        //- Write pressure-constraint diagnostic fields at output times
        const Switch writePressureConstraintDiagnostics_;

        //- Write geometric Jacobian diagnostic fields at output times
        const Switch writeGeometricJacobianDiagnostics_;

        //- Write quasi-static diagnostic fields at output times
        const Switch writeQuasiStaticDiagnostics_;

        //- Report extended incompressibility diagnostics
        const Switch reportExtendedIncompressibilityDiagnostics_;

        //- Threshold used for incompressibility bad-cell masks
        const scalar incompressibilityBadCellThreshold_;

        //- Enable optional Land Problem 3 operator diagnostics
        const Switch landProblem3DiagnosticsEnabled_;

        //- Write exact momentum residual decomposition diagnostics
        const Switch landWriteMomentumDecomposition_;

        //- Write exact pressure residual decomposition diagnostics
        const Switch landWritePressureDecomposition_;

        //- Write cell/face/geometric Jacobian comparison diagnostics
        const Switch landWriteJacobianComparison_;

        //- Write/report fixed-base audit diagnostics
        const Switch landWriteBasalAudit_;

        //- Report follower-pressure force and moment checks
        const Switch landWritePressureTractionCheck_;

        //- Run Land diagnostics at every completed time step
        const Switch landWriteEveryTimeStep_;

        //- Run Land diagnostics at final time
        const Switch landWriteAtFinalTime_;

        //- Basal z reference used by the fixed-base audit
        const scalar landBasalZReference_;

        //- Basal z tolerance used by the fixed-base audit
        const scalar landBasalZTolerance_;

        //- Enable prescribed-field momentum-stabilisation consistency tests
        const Switch momentumStabilisationConsistencyTestEnabled_;

        //- Write prescribed-field momentum-stabilisation test fields
        const Switch momentumStabilisationConsistencyWriteFields_;

        //- Run prescribed-field momentum-stabilisation tests at startup
        const Switch momentumStabilisationConsistencyRunAtStartup_;

        //- Enable prescribed affine kinematics consistency test
        const Switch affineKinematicsConsistencyTestEnabled_;

        //- Write prescribed affine kinematics consistency fields
        const Switch affineKinematicsConsistencyWriteFields_;

        //- Run prescribed affine kinematics test at startup
        const Switch affineKinematicsConsistencyRunAtStartup_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Apply the mixed pressure-displacement stress split
        void applyMixedPressureStressSplit();

        //- Allocate the reciprocal of the bulk modulus
        void makeRKappa() const;

        //- Return const access to the reciprocal of the bulk modulus
        const volScalarField& rKappa() const;

        //- Update solidModel::rAUf() if the mesh or deltaT changed
        //  since it was last computed. The result depends only on
        //  impKf_, rho, mesh and deltaT (not on the values of D or p),
        //  so it can safely be reused across MFFD residual evaluations
        //  within a single Newton solve. Does nothing when
        //  solvePressure() is disabled
        void updateRAUfIfStale();

        //- Refresh Jminus1 from the current Jacobian field
        void updateJacobianDiagnostics();

        //- Report incompressibility diagnostics for the current Jacobian
        void reportIncompressibilityDiagnostics() const;

        //- Update and optionally write/report final-state diagnostics
        void updateFinalStateDiagnostics();

        //- Update unscaled pressure-constraint diagnostic fields
        void updatePressureConstraintDiagnostics();

        //- Update geometric Jacobian diagnostic fields
        void updateGeometricJacobianDiagnostics();

        //- Update quasi-static momentum diagnostic fields
        void updateQuasiStaticDiagnostics();

        //- Calculate deformed cell volumes from a point-displacement field
        void calcDeformedCellVolumes
        (
            const pointVectorField& pointD,
            scalarField& deformedCellVolumes
        ) const;

        //- Write pressure-constraint diagnostics when enabled
        void writePressureConstraintDiagnostics() const;

        //- Write geometric Jacobian diagnostics when enabled
        void writeGeometricJacobianDiagnostics() const;

        //- Write quasi-static diagnostics when enabled
        void writeQuasiStaticDiagnostics() const;

        //- Report pressure-constraint term statistics
        void reportPressureConstraintDiagnostics() const;

        //- Report quasi-static diagnostic statistics
        void reportQuasiStaticDiagnostics() const;

        //- Return global active Cauchy stress for decomposition diagnostics
        tmp<volSymmTensorField> activeCauchyStressDiagnostic() const;

        //- Enforce traction boundary conditions on decomposed face forces
        void enforceDecomposedTractionBoundaries
        (
            surfaceVectorField& passiveForce,
            surfaceVectorField& activeForce,
            surfaceVectorField& pressureForce,
            const volVectorField& D,
            const surfaceVectorField& nCurrent,
            const surfaceScalarField& magSfCurrent
        ) const;

        //- Apply the selected physical-boundary treatment to the isolated
        //  momentum-stabilisation face force
        void applyMomentumStabilisationBoundaryTreatment
        (
            surfaceVectorField& stabilisationFaceForce
        ) const;

        //- True if Land diagnostics should run for this completed time step
        bool runLandProblem3DiagnosticsNow() const;

        //- Update optional Land Problem 3 diagnostics
        void updateLandProblem3Diagnostics();

        //- Update exact momentum decomposition diagnostics
        void updateLandMomentumDecompositionDiagnostics() const;

        //- Update exact pressure decomposition diagnostics
        void updateLandPressureDecompositionDiagnostics() const;

        //- Update Jacobian comparison diagnostics
        void updateLandJacobianComparisonDiagnostics() const;

        //- Report/write fixed-base audit diagnostics
        void updateLandBasalAuditDiagnostics() const;

        //- Report follower-pressure force and moment diagnostics
        void updateLandPressureTractionDiagnostics() const;

        //- Run prescribed-field momentum-stabilisation consistency tests
        void runMomentumStabilisationConsistencyTests() const;

        //- Run one prescribed-field momentum-stabilisation consistency test
        void runMomentumStabilisationConsistencyTest
        (
            const word& testName,
            const volVectorField& prescribedD,
            const volTensorField& prescribedGradD,
            const bool shouldVanish
        ) const;

        //- Run prescribed affine kinematics consistency test
        void runAffineKinematicsConsistencyTest() const;

        //- Report a cell location for a global scalar-field extremum
        void reportScalarExtremumLocation
        (
            const word& name,
            const volScalarField& field,
            const scalar extremumValue
        ) const;

#ifdef USE_PETSC
        //- Evaluate the production mixed residual once at an exact affine
        //  state, without entering a Newton iteration. Runtime opt-in only.
        void evaluateExactAffineResidual();

        //- Evaluate the production mixed residual once using DExact/pExact
        //  fields supplied by an MMS case. Runtime opt-in only.
        void evaluateExactMmsResidual();

        //- Report an evaluated SNES line-search candidate. A null residual
        //  denotes an inadmissible candidate caught before constitutive
        //  evaluation.
        void reportSnesLineSearchTrial
        (
            const Vec x,
            const Vec residual,
            const bool forceReport = false
        );

        //- Copy the PETSc solution vector x into the D field (and p
        //  when solvePressure() is active), refresh the dependent
        //  kinematic fields (gradD, DD, gradDD, F, Finv, J, sigma) and
        //  correct boundary conditions. Applies the pressure unknown
        //  scale (p = pressureUnknownScale_ * pHat). Centralises the
        //  mapping used by formResidual, formJacobian and evolveSnes.
        //  Returns false only for a recoverable inadmissible PETSc
        //  line-search candidate.
        //  Vec is a pointer type, so it is taken by value
        bool unpackSolution(const Vec x, const bool allowDomainError);

        //- Copy the current D field (and, when solvePressure() is
        //  active, p) into the PETSc solution vector x, applying the
        //  inverse pressure unknown scale (pHat = p/pressureUnknownScale_).
        //  Vec is a pointer type, so it is taken by value
        void packSolution(Vec x);
#endif

        //- Enforce traction boundary conditions on the provided traction field
        void enforceTractionBoundaries
        (
            surfaceVectorField& force,             // force at the faces
            const volVectorField& D,               // displacement field
            const surfaceVectorField& nCurrent,    // deformed unit normals
            const surfaceScalarField& magSfCurrent // deformed area magnitudes
        ) const;

        //- Return the intensive cell force from physical stress plus
        //  momentum stabilisation after traction boundary enforcement
        tmp<volVectorField> momentumSurfaceForceDensity
        (
            const volVectorField& D,
            const surfaceVectorField& nCurrent,
            const surfaceScalarField& magSfCurrent
        ) const;

        //- Report/write optional momentum stabilisation diagnostics
        void reportMomentumStabilisationDiagnostics
        (
            const surfaceVectorField& momentumStabFaceForce,
            const vectorField& momentumStabCellForce,
            const volVectorField& physicalSurfaceForceDensity
        ) const;

        //- Evolve the solid solver using the segregated implicit approach
        bool evolveImplicitSegregated();

        //- Evolve the solid solver using the PETSc nonlinear solver (SNES)
        //  A compact approximate Jacobian (same one as used by the segregated
        //  approach) is provided to PETSc. A Jacobian-free Krylov-Newton method
        //  works well with this approach
        bool evolveSnes();

        //- Disallow default bitwise copy construct
        nonLinGeomTotalLagTotalDispSolid
        (
            const nonLinGeomTotalLagTotalDispSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomTotalLagTotalDispSolid&);


protected:

    // Protected member functions

        //- True when the mixed split must retain the complete passive stress
        //  The generic model deliberately returns false to preserve its
        //  established dev(sigmaPassive) - p*I contract.
        virtual bool retainFullPassiveStressInMixedSplit() const
        {
            return false;
        }

        //- Complete a direct constitutive face Cauchy stress.
        //  The generic no-op preserves all existing direct-face laws.
        virtual void completeDirectFaceCauchyStress
        (
            surfaceSymmTensorField& sigmaFace
        ) const;

        //- Return the positive finite-strain constraint g(J)
        //  The pressure residual contains -g(J).
        virtual tmp<volScalarField> mixedVolumetricConstraint
        (
            const volScalarField& J
        ) const;

        //- Scalar form of g(J), used to document and test derivatives
        virtual scalar mixedVolumetricConstraintValue(const scalar J) const;

        //- Exact scalar derivative dg/dJ
        virtual scalar mixedVolumetricConstraintDerivative
        (
            const scalar J
        ) const;

        //- Human-readable g(J) expression for diagnostics
        virtual string mixedVolumetricConstraintDescription() const;

        //- Refresh pressure diagnostics for an opt-in derived model
        void refreshMixedPressureConstraintDiagnostics()
        {
            updatePressureConstraintDiagnostics();
        }

        //- Read-only pressure diagnostic fields for an opt-in derived model
        const volScalarField& pressureCompressibilityDiagnostic() const
        {
            return pressureCompressibilityTerm_;
        }

        const volScalarField& pressureStabilisationDiagnostic() const
        {
            return pressureStabilisationTerm_;
        }

        const volScalarField& volumetricConstraintDiagnostic() const
        {
            return volumetricConstraintTerm_;
        }

        const volScalarField& completePressureResidualDiagnostic() const
        {
            return pressureConstraintResidual_;
        }

        //- Construct with an explicit runtime model name. This overload
        //  permits an isolated derived runtime model to own its coefficients
        //  dictionary while leaving the existing constructor unchanged.
        nonLinGeomTotalLagTotalDispSolid
        (
            const word& modelType,
            Time& runTime,
            const word& region
        );

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Audit-only decomposition of the current momentum residual.  The
    //  decomposition is populated from the same face-force path used by the
    //  production residual; it does not assemble or alter the PETSc matrix.
    struct momentumResidualDecompositionData
    {
        wordList names;
        List<vectorField> cellActions;
        List<vectorField> faceForces;
        symmTensorField passiveFaceStress;
        symmTensorField viscousFaceStress;
        symmTensorField pressureFaceStress;
        symmTensorField directPassiveFaceStress;
        symmTensorField directViscousFaceStress;
        vectorField faceAreaVectors;
    };

    //- Runtime type information
    TypeName("nonLinearGeometryTotalLagrangianTotalDisplacement");

    // Constructors

        //- Construct from components
        nonLinGeomTotalLagTotalDispSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~nonLinGeomTotalLagTotalDispSolid()
        {}


    // Member Functions

        //- Populate the audit-only term-by-term momentum decomposition at the
        //  current trial state.  The returned cell actions are intensive
        //  residual contributions and the returned face forces use the global
        //  polyMesh face numbering.
        void momentumResidualDecomposition
        (
            momentumResidualDecompositionData& decomposition
        ) const;

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

#ifdef USE_PETSC
            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            virtual label initialiseJacobian(Mat& jac);

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x);

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x    // Solution
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            );
#endif // USE_PETSC

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
