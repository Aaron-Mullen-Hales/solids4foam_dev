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

        //- Surface (face) gradient of the total displacement, registered as
        //  "grad(D)f". This is the field mechanicalLaw::updateFf() and
        //  mechanicalLaw::updateF(surfaceSymmTensorField&, ...) look up, and
        //  it is what makes the mechanical law's face deformation gradient
        //  Ff() authoritative rather than a stale identity
        surfaceTensorField gradDf_;

        //- False when the most recent unpackSolution(x) produced a trial
        //  deformation gradient that is not admissible; see
        //  checkTrialGeometry() for the exact condition. Some mechanical
        //  laws have no volumetric barrier term to keep Newton's trial
        //  steps away from this region, so this flag lets formResidual()
        //  and formJacobian() report the failure through PETSc's own
        //  domain-error mechanism instead of the mechanical law evaluating
        //  stress on invalid geometry.
        bool trialGeometryValid_;

        //- Smallest J_ of the most recent trial state. Retained so that a
        //  rejected trial can be reported quantitatively
        scalar trialMinJ_;

        //- Number of trial states rejected by checkTrialGeometry() since the
        //  start of the run. Diagnostic only
        label trialRejectionCount_;

        //- Number of faces whose fvc::fGrad deformation gradient was
        //  inadmissible and replaced by the cell displacement gradient
        label faceSubstitutionCount_;

        //- Number of times the face kinematics were actually rebuilt
        label faceRefreshCount_;

        //- Number of residual/Jacobian probes that skipped the rebuild
        label faceRefreshSkipCount_;

        //- Emit the solidDirectionMixed uncorrected-patch warning only once
        bool warnedUncorrectedPatches_;

        //- True when the most recent nonlinear solve ended with a negative
        //  SNES reason and its result was therefore not accepted. Used to
        //  stop updateTotalFields() advancing the accumulated history as
        //  though the step had converged
        bool solveFailed_;

        //- Constitutive stress retained without a deviatoric projection
        volSymmTensorField sigmaPreserved_;

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

        //- Viscous implicit stiffness surface field, used ONLY as a
        //  preconditioner contribution in formJacobian() (Change 6).
        //  Rebuilt whenever deltaT changes, since it scales with the
        //  time-derivative coefficient; null when the law is
        //  rate-independent or the contribution is switched off
        mutable autoPtr<surfaceScalarField> viscousImpKfPtr_;

        //- deltaT for which viscousImpKfPtr_ was last built
        mutable scalar viscousImpKfDeltaT_;

        //- Enable the viscous P_DD contribution (default on when the law
        //  reports a non-zero viscous stiffness)
        const Switch viscousJacobian_;

        //- How the internal momentum face force is built in formResidual().
        //
        //    interpolatedCell  (DEFAULT -- existing behaviour, unchanged)
        //        Interpolate the CELL stress and the CELL J*Finv^T to faces:
        //          SfCurrent = I_h[J Finv^T] & Sf
        //          force     = |SfCurrent| * (nCurrent & I_h[sigma])
        //
        //    directFaceState   (EXPERIMENTAL, diagnostic)
        //        Build the face state directly from the corrected face
        //        kinematics rather than interpolating cell quantities:
        //          Ff        = I + gradDf^T,  Jf = det(Ff),  Ff^-1
        //          sigma_f   = material face stress - I_h[p]*I
        //          SfCurrent = Jf Ff^-T & Sf
        //          force     = SfCurrent & sigma_f
        //
        //  Only the INTERNAL face force differs. Traction-boundary
        //  enforcement and momentum stabilisation are identical in both, as
        //  are the pressure residual, constitutive law, BCs, time scheme and
        //  the preconditioner.
        const word faceForceTreatment_;

        //- Compare the two face-force constructions each residual evaluation
        //  and report the largest relative difference and where it occurs.
        //  Diagnostic only: both are computed, the selected one is used
        const Switch faceForceDiagnostics_;

        //- Faces for which admissibleDeformation() rejected the direct face
        //  state in the directFaceState path, so the owner-cell state was
        //  substituted. Never silent
        mutable label directFaceFallbackCount_;

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

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

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

        //- Build/refresh the viscous implicit stiffness used by
        //  formJacobian(). It depends only on the law's viscosity and
        //  deltaT, so it is rebuilt only when deltaT changes and is
        //  reused across every Jacobian evaluation in a Newton solve.
        //  Returns nullptr when there is no viscous contribution.
        const surfaceScalarField* viscousImpKfIfAny();

        //- Report the largest and rms difference between the two face-force
        //  constructions, and the owner/neighbour/location of the worst face
        void reportFaceForceDifference
        (
            const surfaceVectorField& forceA,
            const surfaceVectorField& forceB
        ) const;

        //- Refresh the face kinematics: interpolate D to the vertices,
        //  evaluate the face displacement gradient with the established
        //  solids4foam face-gradient operator (fvc::fGrad via
        //  mechanicalModel::grad(D, pointD, gradDf)) and hand it to the
        //  mechanical law so that Ff() = I + gradDf^T is up to date.
        //
        //  The momentum and pressure residuals do not read Ff(); this exists
        //  so that face-based consumers (face viscous fields, the nominal and
        //  material tangents, and the active face stress) never see a stale
        //  identity
        //- Rebuild pointD, grad(D)f and the law's Ff from the current D, make
        //  the traction-boundary normal row kinematic, and guard admissibility
        void updateFaceKinematics();

        //- Decide whether the current F_ is an admissible deformation
        //  gradient, and set J_ = det(F_) while doing so. Returns false for
        //  a state that is mathematically or numerically invalid:
        //
        //    - any component of F_ non-finite;
        //    - J_ non-finite;
        //    - J_ <= 0, i.e. a locally inverted (orientation-reversing)
        //      deformation, for which the constitutive law's J^(-2/3) is
        //      undefined;
        //    - J_ too small for inv(F_) = adj(F_)/det(F_) to be
        //      representable.
        //
        //  No physical deformation limit is imposed: an arbitrarily large
        //  but finite, positive-J_ deformation is admissible. A merely
        //  ill-conditioned F_ with positive J_ is also admissible; the
        //  resulting large residual is handled by the line search.
        //
        //  F_ is screened for finiteness before the determinant is formed,
        //  because det() of a tensor containing Inf can evaluate Inf - Inf and
        //  raise FE_INVALID under FOAM_SIGFPE. The screen uses only
        //  std::isfinite, which performs no arithmetic.
        //
        //  On success J_ and Finv_ are filled, for internal cells AND for every
        //  boundary patch, with the determinant and the genuine 3-D inverse
        //  adj(F)/det(F). Foam::det()/Foam::inv() are not used: Foam::inv()
        //  applies Tensor::safeInv(), which can return a finite but incorrect
        //  surrogate for a valid deformation gradient (e.g. an exact 90-degree
        //  rotation). See admissibleDeformation() in the .C file
        bool checkTrialGeometry();

        //- Human-readable reason the most recent trial geometry was rejected.
        //  Distinguishes "F itself was rejected, so no determinant exists"
        //  from a genuine non-positive/non-finite determinant
        string trialGeometryDescription() const;

#ifdef USE_PETSC
        //- Copy the PETSc solution vector x into the D field (and p
        //  when solvePressure() is active), refresh the dependent
        //  kinematic fields (gradD, DD, gradDD, F, Finv, J, sigma) and
        //  correct boundary conditions. Applies the pressure unknown
        //  scale (p = pressureUnknownScale_ * pHat). Centralises the
        //  mapping used by formResidual, formJacobian and evolveSnes.
        //  Vec is a pointer type, so it is taken by value
        //  refreshFaceState rebuilds the face kinematics (pointD, grad(D)f,
        //  Ff) before the mechanical law is evaluated. It is false for ordinary
        //  MFFD residual and Jacobian probes, which do not consume Ff, and true
        //  for the final accepted unpack, whose face state is what
        //  updateTotalFields() shifts into the viscous face history.
        void unpackSolution(const Vec x, const bool refreshFaceState = false);

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

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

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
            //- Update total accumulated fields. Refuses to advance the
            //  accumulated history after a nonlinear solve whose result was
            //  not accepted
            virtual void updateTotalFields();

            //- True when the most recent evolve() did not converge and its
            //  result was therefore not accepted. The time loop uses this to
            //  avoid accepting the attempted time step
            virtual bool solutionFailed() const;

            //- Report the face-kinematics lifecycle counters
            virtual void end();

            //- Rebuild the face kinematics (pointD, grad(D)f and the law's Ff)
            //  from the current displacement field.
            //
            //  Ordinary MFFD residual and Jacobian probes deliberately do NOT
            //  do this, because the residual is built from cell quantities and
            //  cannot observe Ff. Any routine that evaluates a face
            //  constitutive quantity or a face tangent dP_f/dF_f must therefore
            //  call this first, at the state it is linearising about, or it
            //  will read the face state left by the last accepted step.
            void refreshFaceKinematics();

            //- This model's formJacobian() can report a Jacobian-domain
            //  error, so PETSc's Jacobian-domain check must be enabled
            virtual bool reportsJacobianDomainError() const;

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
