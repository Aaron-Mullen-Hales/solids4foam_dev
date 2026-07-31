# Aróstica monoventricle mechanics: source-architecture audit

Date: 2026-07-27  
Repository: solids4foam  
Branch: solid-pressure-land3Fix  
Commit: ebf998ff25945f611a44e5376656690c62768d32 (Fix electro-mechanical deformation gradient ownership)

The working tree was already dirty when this audit began, including changes in several of the files inspected below. The checked-out working-tree contents, not an older copied implementation, are the authority for the source findings. No source or accepted case file was changed by this audit.

The accepted case was located at
/Volumes/OpenFoam/aaronmullen-hales-v2312/run/CardiacMechanics/ventricle/arostica.
The mesh/fibre case was left untouched.

## Executive result

The existing total-Lagrangian PETSc model supplies the required nonlinear displacement/pressure architecture, inertia, follower-pressure geometry, stabilisation, histories, and domain-safe trial handling. It does not supply the Aróstica constitutive law, viscosity, activation history, or the two support boundary conditions.

The exact mixed stress split requires a minimal derived solid model:

    class ArosticaNonLinGeomTotalLagTotalDispSolid
    :
        public nonLinGeomTotalLagTotalDispSolid
    {
        virtual bool retainFullPassiveStressInMixedSplit() const { return true; }
    };

This is a plan, not an implementation recommendation to copy the base solid model. The current base class already provides the protected runtime-name constructor needed by such a derived runtime model. The derived class must reuse the complete base residual and override only the split contract unless testing exposes an unrelated boundary-hook requirement.

The new passive law should return the complete non-volumetric Cauchy stress:

    sigma_law = sigma_isochoric + sigma_fibre + sigma_sheet + sigma_fs
                + sigma_viscous + sigma_active   [when wrapped]

with no U(J) contribution. The derived solid model then applies sigma = sigma_law - p I. This preserves the anisotropic spherical part and inserts the volumetric stress exactly once.

## Evidence and reference sources

Inspected source files include:

- src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/nonLinGeomTotalLagTotalDispSolid.{H,C}
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.H
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalModel.{H,C}
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.{H,C}
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GuccioneElastic/GuccioneElastic.{H,C}
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/HolzapfelGasserOgdenElastic/HolzapfelGasserOgdenElastic.{H,C}
- src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic/GultekinTwoFibreElastic.{H,C}
- src/solids4FoamModels/solidModels/fvPatchFields/solidTraction/solidTractionFvPatchVectorField.C
- src/solids4FoamModels/solidModels/fvPatchFields/pressureTraction/pressureTractionFvPatchScalarField.C
- src/solids4FoamModels/solidModels/solidModel/solidModel.{H,C}
- src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/README.md

The pinned repository in the accepted case is
reference/cardiac_benchmark_toolkit at commit
e8d47553cfc83eb274eba3e177de33148e7f441c. It is a provenance/data toolkit for meshes, fibres, and benchmark fields; it does not contain the mechanics runtime law or boundary-condition implementation. The useful independent constitutive cross-check is the local reference implementation at:

/Volumes/OpenFoam/aaronmullen-hales-v2312/run/PHDS1/cardiacFoam_backup/tutorials/LandEtAl2015/problem3/problem2/run_Apple_M4_Pro_20250320_154952/paper2/finsberg-cardiac_benchmark-bb849fb/src/cardiac_benchmark/material.py

with activation in its adjacent activation_model.py. The accepted case also contains a standalone src/arosticaBenchmarkModels/ArosticaConstitutiveMath scaffold. That scaffold is not a project runtime law and is incomplete for the sheet compression switch; it must not be treated as the current source implementation.

## Exact Aróstica equations and parameters

The CMAME paper defines

    F = I + Grad(u),       J = det(F),       C = F^T F,
    E = 1/2 (C - I),       S = J F^-1 T F^-T.

The material is stated in second-Piola form as

    S(t) = d Psi_aniso/dE + d Psi_visc/d(dot E) + tau(t) f0 tensor f0.

The invariants are

    I1bar = J^(-2/3) tr(C)
    I4f   = f0 dot C dot f0
    I4s   = s0 dot C dot s0
    I8fs  = f0 dot C dot s0.

The passive energy is

    Psi_aniso = a/(2b)       ( exp(b (I1bar - 3)) - 1 )
              + af/(2bf)     chi(I4f) ( exp(bf  (I4f - 1)^2) - 1 )
              + as/(2bs)     chi(I4s) ( exp(bs  (I4s - 1)^2) - 1 )
              + afs/(2bfs)   ( exp(bfs I8fs^2) - 1 )
              + kappa/4      (J^2 - 1 - 2 ln J).

The exact paper switch is

    chi(x) = x, if x > 1;
             0, otherwise,       x in R+.

It is therefore neither a discontinuous indicator nor merely a positive-part factor: above the threshold it multiplies the exponential by I4, and in compression it is zero. The paper suggests, separately, the logistic approximation 1/(1 + exp(-k (x - 1))). The independent Finsberg reference uses that logistic by default for both the fibre and sheet terms, with k=100; its discontinuous alternative is an indicator-style conditional, not the exact paper chi(x)=x. The accepted-case scaffold currently uses a logistic or a fibre-only I4f > 1 ? I4f : 0 branch, so it is not a complete implementation of the paper's two-term switch. This discrepancy must be made an explicit dictionary choice rather than hidden in a generic heaviside utility.

The volumetric energy is

    U(J) = kappa/4 (J^2 - 1 - 2 ln J).

The viscous energy and stress are

    Psi_visc = eta/2 tr(dot(E)^2),
    S_visc   = d Psi_visc/d(dot E) = eta dot(E).

The active second-Piola stress is

    S_active = tau(t) f0 tensor f0,
    dot(tau) = -|a(t)| tau + sigma0 |a(t)|_+,
    |a|_+ = max(a, 0),
    a(t) = alphaMax f(t) + alphaMin (1 - f(t)),
    f(t) = S+(t - tSys) S-(t - tDias),
    S+(q) = 1/2 (1 + tanh(q/gamma)),
    S-(q) = 1/2 (1 - tanh(q/gamma)).

The pressure waveform in the paper is a separate ODE for the benchmark pressure case:

    dot(p) = -|b(t)| p + sigmaMid |b(t)|_+ + sigmaPre |gPre(t)|_+,
    b(t) = aPre(t) + alphaPre gPre(t) + alphaMid,
    aPre(t) = alphaMax fPre(t) + alphaMin (1 - fPre(t)),
    fPre(t) = S+(t - tSysPre) S-(t - tDiasPre),
    gPre(t) = S-(t - tDiasPre).

The benchmark values are:

| quantity | value |
|---|---:|
| rho | 1000 kg/m3 |
| eta | 100 Pa s |
| kappa | 1e6 Pa |
| compression k | 100 |
| a, af, afs, as | 59.0, 18472, 216, 2481 Pa |
| b, bf, bfs, bs | 8.023, 16.026, 11.436, 11.12 |
| sigma0 | 1.5e5 Pa |
| gamma | 0.005 s |
| alphaMin, alphaMax | -30, 5 |
| tSys, tDias | 0.16, 0.484 s |
| alphaEpi, betaEpi | 1e8 Pa/m, 5e3 Pa s/m |
| alphaTop, betaTop | 1e5 Pa/m, 5e3 Pa s/m |
| pressure alphaPre, alphaMid, sigmaPre, sigmaMid | 5, 1, 7000, 16000 Pa as used by the paper's pressure ODE |
| pressure tSysPre, tDiasPre, gamma | 0.17, 0.484, 0.005 s |

The independent reference solves the active ODE with a stiff Radau integrator from tau(0)=0; it does not advance the state once per constitutive call. The paper reports a maximum active tension of approximately 118817.07 Pa and pressure of approximately 16117.52 Pa for its reference curves.

The strong-form boundary and inertia requirements are:

    rho ddot(u) - Div(J T F^-T) = 0,
    endocardium: P N = p J F^-T N,
    epicardium: (P N) dot N + alphaEpi (u dot N)
                + betaEpi (dot(u) dot N) = 0,
                (P N) cross N = 0,
    base:       P N + alphaTop u + betaTop dot(u) = 0.

Here P=J T F^-T; the support expressions are nominal/reference-area tractions. That measure matters for the custom support boundary conditions.

## Existing solid-model audit

### PETSc, mixed vector, residual, histories, and trial safety

nonLinGeomTotalLagTotalDispSolid registers the total-Lagrangian total-displacement runtime model and has a PETSc SNES path. With solvePressure true, the solution vector contains the displacement components plus one pressure unknown per cell; the case uses the PETSc path and mixed-field scaling.

The production residual is assembled through formResidual after unpackSolution. The latter reconstructs grad(D), sets

    F = I + grad(D)^T,  J = det(F),  Finv = inv(F),

rejects non-finite or non-positive trial Jacobians, calls the mechanical law, then inserts the mixed pressure. The SNES/MFFD paths therefore repeatedly call the same trial evaluation. The current model has a domain-safe line-search option and returns a recoverable error for invalid trial states. A future viscoelastic or activation provider must not write accepted history from any of these paths.

The momentum residual includes the constitutive surface force, body force, -rho*d2dt2(D), optional -rho*dampingCoeff*ddt(D), and momentum stabilisation. The pressure residual is an extensive, scaled form of

    -p/kappa - g(J) + pressureStabilisation = 0.

The pressure row has the corresponding analytic dg/dJ contribution, while the PETSc Jacobian is an approximate/preconditioner structure; SNES still uses the complete nonlinear residual and can use matrix-free derivatives.

The constructor forces the required old-time fields, and the solid model maintains old-time D, DD, gradD, gradDD, sigma, and p structures. The selected case has ddtSchemes { default backward; }; this is OpenFOAM backward finite-difference time integration through fvm/fvc ddt and d2dt2, not a generalized-alpha implementation. It supplies transient inertia, but the benchmark reproduction must document the resulting time discretisation and use the benchmark time step deliberately.

### Volumetric constraint and pressure sign

The source implements

    g(J) = 0.5 (J^2 - 1)/J = 0.5 (J - 1/J),
    g'(J) = 0.5 (1 + 1/J^2).

For the paper's volumetric energy,

    dU/dJ = kappa/4 (2J - 2/J)
           = kappa/2 (J - 1/J)
           = kappa g(J).

Using dJ/dC = J/2 C^-1,

    S_U = 2 (dU/dJ)(dJ/dC)
        = (dU/dJ) J C^-1,
    sigma_U = (1/J) F S_U F^T
             = (dU/dJ) I
             = kappa/2 (J - 1/J) I.

The source pressure row without stabilisation is

    -p/kappa - g(J) = 0,

so

    p = -kappa g(J) = -kappa/2 (J - 1/J),
    -p I = kappa/2 (J - 1/J) I = sigma_U.

The sign is therefore correct: the pressure replacement reproduces the derivative of U(J) exactly, including the sign, when pressure is inserted as -p I.

### Mixed stress split: the architectural decision

The current default path does this for a single passive law:

    sigma <- dev(sigma_law) - p I.

For an electro-mechanical law it preserves the active part separately:

    sigma <- dev(sigma_law - sigma_active) + sigma_active - p I.

That default is appropriate only when the law's stress entering the mixed model has no spherical non-volumetric component, or when the intended model explicitly defines a deviatoric split. Aróstica does not satisfy that premise. I1bar is isochoric, but I4f, I4s, and I8fs are not all made traceless by the J^(-2/3) operation. The pushed-forward dyads from the fibre, sheet, and fibre-sheet derivatives generally have non-zero trace. Applying dev would remove real benchmark stress and cannot be repaired by -p I.

Therefore the exact choice is neither “full passive stress including U(J) plus another -pI” nor the generic deviatoric option. It is:

    sigma_exact = sigma_nonvolumetric_passive
                  + sigma_viscous
                  + sigma_active
                  - p I,

where sigma_nonvolumetric_passive contains the complete traces of the anisotropic derivatives and does not contain sigma_U.

The current protected retainFullPassiveStressInMixedSplit() hook returns false, while applyMixedPressureStressSplit() already has the true branch sigma = sigma - p I. The protected constructor overload taking (modelType, runTime, region) is present specifically to allow an isolated derived runtime model to own its dictionary. This proves that the minimal derived ArosticaNonLinGeomTotalLagTotalDispSolid route is supported by the current architecture. No copy of the approximately 7000-line base source is needed.

## Equation-to-code mapping table

Classification meanings: exact means the current implementation satisfies the benchmark contract; reusable with configuration means the mechanism is present but requires benchmark dictionaries; reusable framework only means the surrounding interface is useful but the required physics is absent; and missing means a new runtime component is required.

| component | exact benchmark requirement | current source implementation | classification | required change | proposed source files | tests required |
|---|---|---|---|---|---|---|
| mixed solid model | total-Lagrangian nonlinear D plus cell pressure unknown | PETSc SNES, mixed vector, trial unpacking, residual and approximate Jacobian are present; protected full-stress hook and constructor overload are present | reusable with configuration | add the minimal derived runtime model returning true from retainFullPassiveStressInMixedSplit() | new solidModels/ArosticaNonLinGeomTotalLagTotalDispSolid/{H,C}, Make/files | runtime selection; mixed affine deformation; exact pressure row; SNES/MFFD repeated residuals |
| stress split | complete non-volumetric passive/viscous/active stress minus pI | default passive path applies dev(sigma); special active path only preserves active stress | missing for exact benchmark | use the derived full-stress branch and ensure the law excludes U(J) | derived solid model; new Aróstica law | non-zero-trace fibre/sheet/shear stress; compare A/B stress decomposition |
| volumetric constraint | U=kappa/4(J²-1-2lnJ) replaced by -pI, with p=-kappa g(J) | g=0.5(J-1/J), residual -p/kappa-g+stab, exact derivative | exact | return benchmark kappa and disable law volumetric insertion in mixed mode | existing model; new law dictionary | scalar J sign/derivative test; p/sigma_U identity; J residual with stabilisation off |
| passive material | exact four-term Holzapfel-Ogden energy, supplied triad fields, exact/specified chi | no Aróstica law; current case names ArosticaHolzapfelOgdenElastic, but no matching project source was found | missing | add independent runtime-selectable law with explicit dictionary and switch convention | new ArosticaHolzapfelOgdenElastic.{H,C}, Make/files | one-cell invariant/energy/stress reference; logistic and exact-switch branches; supplied triad preservation |
| viscosity | S_visc=eta dot(E) from current trial and accepted old-time state | existing visco law patterns use old-time internal fields, but no Aróstica dotE law | reusable framework only | add trial-safe Edot evaluation and push-forward; do not mutate fields in correct() | new law; optional diagnostic fields | backward-step E-Eold; repeated residual/MFFD/rejected-line-search history hash |
| active stress wrapper | F (tau f0 tensor f0) F^T/J in cell and face paths | electroMechanicalLaw wraps a runtime passive law, owns/forwards F/Ff, reads f0/f0f, forms dyads, supports field Ta, and correctly push-forwards cell/face stress | reusable with configuration | use unchanged for Aróstica once Ta is supplied; do not use constant ramp as benchmark activation | existing electroMechanicalLaw unchanged; new provider | cell/face rank-one push-forward; identical F ownership; no registry/history mutation during residual calls |
| time-dependent Ta provider | deterministic activation ODE/history, restartable, current-time consistent, no residual-time advancement | constant/ramp plus object-registry field lookup; no current Land activation ODE provider | reusable framework only | add ArosticaActiveTension with accepted-time update, AUTO_WRITE Ta, old-time state, and deterministic cell/face values; permit a reference-generated tabulated fallback | new activation/coupling model and Make/files; optional table reader | reproduce tau(t) and max; restart equivalence; repeated-call idempotence; update exactly once per accepted step |
| follower pressure | PN=p J F^-T N, time-dependent Case B pressure | solidTraction has pressure series/field and the total-Lagrangian model forms current area vectors J F^-T Sf; case uses useUndeformedArea false | reusable with configuration, pending sign test | configure pressure series and verify patch-normal sign; do not use activation ODE for pressure unless reproducing paper pressure ODE | existing solidTraction/base model unchanged; case-only configuration later | undeformed/deformed area force identity; constant pressure force/moment; patch orientation/sign |
| epicardial support | nominal normal traction -[alpha(u dot N)+beta(dot u dot N)]N, zero tangent | no exact normal spring/dashpot runtime BC found | missing | add small BC that evaluates current D and U without state; enforce nominal/reference measure exactly | new arosticaNormalSpringDashpotTraction.{H,C} | pure normal displacement; pure tangential displacement; velocity damping; force equals reference-area expression |
| basal support | nominal vector traction -alphaTop D-betaTop U | no exact full-vector spring/dashpot runtime BC found | missing | add small stateless vector spring/dashpot BC | new arosticaVectorSpringDashpotTraction.{H,C} | each Cartesian component; zero displacement/velocity; time-step and restart consistency |
| inertia/time integration | rho ddot u, benchmark time step and specified scheme | rho*fvm/fvc*d2dt2(D) and optional damping; accepted old-time fields; case backward ddt | reusable with configuration | set/record the benchmark time step and selected OpenFOAM second-derivative scheme; do not claim generalized-alpha | existing solid model/case system/fvSchemes | manufactured acceleration; restart; time-step convergence |
| JFNK state safety | residual calls must be pure with respect to accepted history | SNES residual, MFFD perturbations, domain-safe Jacobian trials, and line-search invalid-J handling are present | reusable with configuration | make law/provider read-only during trial calls; commit viscous/activation state only after accepted step | new law/provider; existing solid model unchanged | identical residual on repeated calls; rejected line search leaves Eold/Taold unchanged |
| cell stress | second-Piola derivatives pushed to Cauchy consistently | mechanical-law correct(volSymmTensorField); electro wrapper adds active cell stress | reusable framework only | implement exact S derivatives, full non-volumetric stress, and push-forward | new Aróstica law | analytic finite-difference stress; cell affine deformation; trace audit |
| face stress | same stress and active term at faces with accepted face fibre/sheet fields | mechanical model has cell and surface correct; direct constitutive face path exists; electro wrapper has face active stress | reusable framework only | read f0f/s0f/n0f and evaluate face law directly; avoid reconstructing triads or relying only on cell interpolation for validation | new Aróstica law | cell/face consistency on affine field; face traction patch test; supplied face fields |
| tangent/preconditioner support | stable mixed solve; exact tangent optional but useful | mechanicalLaw exposes impK, bulk/shear and material-tangent hooks; existing laws use approximate/numerical infrastructure; PETSc Jacobian is structural | reusable framework only | return positive representative impK/bulk modulus; start with robust approximate tangent, then add constitutive tangent if needed | new law; existing base interfaces unchanged | tangent finite difference; positive impK; SNES convergence and pressure conditioning |

## Passive-law implementation plan

### Milestone 1: purely elastic benchmark law

Proposed name: ArosticaHolzapfelOgdenElastic.

The law should inherit mechanicalLaw, register through the existing nonlinear geometry runtime table, own the deformation-gradient fields through updateF, and implement the established cell and surface correct methods. It should read the accepted fields directly:

    f0, s0, n0, f0f, s0f, n0f.

n0 is useful for validation/orthogonality diagnostics even though the stated energy uses f0 and s0; it must not be regenerated from f0. Suggested dictionary entries are the eight material coefficients, kappa, rho, fibreCompressionSwitch (exact or logistic), kCompressionSwitch, field names for the six cell/face triad fields, includeVolumetricEnergy false in the mixed benchmark, and an explicit bulkModulus/impK policy.

The cell algorithm is:

1. obtain the current trial F from updateF, calculate J and C;
2. calculate I1bar, I4f, I4s, I8fs;
3. differentiate the four non-volumetric energy terms with respect to C;
4. form S=2*dPsi/dC and push it forward as sigma=(F&S&F.T())/J;
5. return the complete non-volumetric Cauchy stress, including its trace, with no U(J) term when used by the mixed model.

The invariant derivatives needed by the implementation are:

    d I1bar/dC = J^(-2/3) [ I - (tr C)/3 C^-1 ],
    d I4f/dC   = f0 tensor f0,
    d I4s/dC   = s0 tensor s0,
    d I8fs/dC  = symm(f0 tensor s0).

For each fibre/sheet term let

    q_i = I4i - 1,
    W_i = ai/(2bi) chi(I4i) exp(bi q_i^2)   [minus constant irrelevant to stress].

For the exact paper switch away from the threshold,

    dW_i/dI4i = ai/(2bi) [chi'(I4i) exp(bi q_i^2)
                           + chi(I4i) exp(bi q_i^2) 2bi q_i],

with chi=I4 and chi'=1 in tension, and both zero in compression. At the threshold the exact law is non-smooth; the implementation must document the chosen one-sided/subgradient convention. For the logistic branch, chi'=k chi(1-chi) and the expression is smooth. The fibre-sheet derivative is

    dW_fs/dC = afs I8fs exp(bfs I8fs^2) symm(f0 tensor s0).

The isotropic term is

    dW1/dC = a/2 exp(b(I1bar-3)) dI1bar/dC.

The face algorithm must repeat this process with Ff and the supplied f0f/s0f/n0f, not with an interpolated or reconstructed basis. The wrapper's active face path can then use the same Ff ownership.

Return bulkModulus() as the finite positive benchmark kappa used by the mixed constraint. Return a positive representative shearModulus() and impK() suitable for the existing stabilisation/preconditioner; a conservative benchmark-scale value is acceptable for the first milestone if documented and validated. Do not include U(J) in correct() when the mixed solid model owns pressure. If a standalone compressible mode is offered, make its volumetric stress a separate explicit mode so it cannot be double-inserted.

Material tangent work should initially reuse the base/numerical tangent interface with a finite-difference-safe, read-only constitutive evaluation. An analytic tangent can follow once the affine and one-cell tests pass. Write diagnostic fields for J, I1bar, I4f, I4s, I8fs, chi_f, chi_s, sigma_nonVol, trace(sigma_nonVol), and the pressure decomposition.

### Milestone 2: add viscosity

Use the current trial field and accepted old-time field:

    E_trial     = 1/2 (F_trial^T F_trial - I),
    E_old      = E at the accepted old time,
    dotE       = (E_trial - E_old)/deltaT,
    S_visc     = eta dotE,
    sigma_visc = F_trial S_visc F_trial^T/J_trial.

F.oldTime() is not the safest source for this material because deformation gradient ownership is inside the mechanical law and its old-time semantics must be established for every cell and face. E.oldTime() is safe only if the new law explicitly owns and writes an accepted E field. D.oldTime() plus recomputing the old F through the same gradient operator is safer than assuming a stored F.oldTime() is available, but it can differ from the accepted gradient if boundary/face reconstruction changes. The preferred implementation is an AUTO_WRITE accepted E field (and a face Ef field if face viscosity is evaluated directly), initialized and committed once per accepted time step. If that is not added, use the existing accepted gradD.oldTime() and reconstruct Fold=I+gradD.oldTime().T() using exactly the same gradient policy as the solid model.

The constitutive call must calculate E_trial and dotE into temporaries only. It must not call storeOldTime, overwrite Eold, or update any Maxwell/ODE state from correct(). Accepted-state advancement belongs to the time-step commit hook after SNES convergence. This is required for repeated residuals, MFFD perturbations, failed line searches, and rejected time steps.

## electroMechanicalLaw and activation plan

The current wrapper is reusable unchanged for the active stress. It:

- selects a nested runtime passive law;
- obtains the same cell/face F owned by that law;
- reads and normalizes f0 and f0f and forms reference dyads;
- computes F & (Ta*f0f0) & F.T()/J for cells and faces;
- supports a registry field named Ta, with face values interpolated from the cell field;
- leaves repeated stress evaluations computational and does not itself mutate accepted state.

The constant activeTension/rampTime path is not the Aróstica ODE. The current source has no production Land Problem 3 activation/coupling model; the landProblem3 code in the selected solid model is diagnostic and residual-decomposition support, not an activation-state integrator. Existing Land cases use constant/ramped active tension and do not solve this benchmark activation history.

The cleanest architecture is a new small ArosticaActiveTension runtime provider, not code inside electroMechanicalLaw. It should register a restartable volScalarField Ta (AUTO_WRITE), provide deterministic values at the current physical time, and expose an explicit accepted-time update(). That update must be called once per accepted time step, with the old accepted Ta as its initial state, and never from correct, formResidual, or a material tangent call. The provider should also maintain Ta.oldTime() and correct boundary values so cell and face calls are consistent.

For the first activation validation, a table provider generated from the reference ODE is a safe fallback. It removes ODE-state mutation from the nonlinear solve, but the table must be pinned with its parameters, time interpolation, initial value, and restart semantics. The eventual provider should integrate the scalar ODE with an accepted-step method and compare its curve against the table/reference values.

## Template audits

### GuccioneElastic

Reusable structure:

- mechanicalLaw inheritance and runtime registration;
- dictionary/dimension validation and updateF ownership;
- cell and surface correct entry points;
- second-Piola storage and push-forward conventions;
- pressure-displacement/volumetric separation as an architectural pattern;
- rho, bulkModulus, shearModulus, impK, and tangent hooks;
- restart fields, validation, and diagnostic conventions.

Not reusable as constitutive code: the current Guccione law constructs s0 and n0 from f0 with a seed-axis basis. It reads f0 and f0f but does not use the accepted Aróstica s0/s0f/n0/n0f fields as the constitutive triad. A new Aróstica law must read all six fields directly. Reconstructing a sheet vector from an arbitrary seed changes I4s; reconstructing the normal can change the orientation/orthogonality used by diagnostics and any sheet/fibre coupling. It therefore changes the sheet stress and I8fs stress, even if f0 is correct. Do not copy Guccione's basis construction or its own pressure insertion.

### HolzapfelGasserOgdenElastic

Reusable boilerplate includes runtime selection, field/dimension handling, updateF, cell/face push-forward structure, bulk/shear/impK, and tangent interfaces. Its constitutive implementation is not the benchmark law. The current HGO law uses cylindrical Ec/Ea/Er fields, two symmetric fibre families, I4/I6, and its own pressure insertion. It has no supplied sheet invariant and no I8fs shear term. It must remain unchanged.

### GultekinTwoFibreElastic

This is the closest project-side pattern for a runtime anisotropic law with supplied fibre fields, finite bulk modulus, cell/face constitutive evaluation, validation, and stabilisation-scale reporting. Reuse its architecture and testing style only. Do not reuse its two-fibre constitutive formula, pressure split, or field names for Aróstica.

## Boundary-condition audit and traction measure

solidTraction is runtime-selectable and supports prescribed vector traction, constant pressure, pressure fields, and time interpolation. In the current total-Lagrangian residual, the deformed area vector is

    Sf_current = interpolate(J Finv^T) & Sf,

and the physical face force is current traction times current area. The tractionBoundarySnGrad path uses the current normal, current Cauchy stress, Finv, and implicit stiffness. Thus a prescribed pressure with useUndeformedArea false is geometrically capable of producing the follower force p J F^-T N; the sign must be checked against the accepted case's endocardial patch-normal orientation because solidTraction forms traction - n*pressure.

The Aróstica epicardial and basal equations are nominal tractions PN, not generic current Cauchy traction vectors. A support BC must either use the existing reference-area force path deliberately (useUndeformedArea true) or convert

    q = PN,
    t_current = q / |J F^-T N|,

before calling the current-traction gradient interface. It must not silently multiply a nominal spring traction by current area. The force and tangent handling should be tested together.

No exact normal spring/dashpot or basal vector spring/dashpot BC was found. Contact penalty models and bulk viscoelastic laws do not express these conditions. The accepted case currently names arosticaNormalSpringDashpotTraction and arosticaVectorSpringDashpotTraction, but those runtime types are absent from the current project source; the case dictionary is a forward-looking scaffold, not evidence that they are available.

## Case A and Case B order

Case A means passive material + viscosity + inertia + both supports + active stress, with zero chamber pressure. Case B means passive material + viscosity + inertia + both supports, with active stress zero and time-dependent chamber pressure.

Complete Case B first. It validates the new passive law, the mixed stress split, viscosity history, inertia/time discretisation, support tractions, and follower pressure without activation-state coupling. Then activate the unchanged wrapper with a prescribed/table Ta, and finally validate the restartable active provider. This order isolates constitutive and boundary errors before adding the activation ODE. A one-cell active-wrapper test can run earlier in parallel, but the first full accepted mesh run should be Case B.

## Proposed implementation order

1. Add unit/reference tests for invariants, chi, stress derivatives, and the U(J)/pressure identity without changing existing laws.
2. Add ArosticaHolzapfelOgdenElastic with supplied cell/face triads and the full non-volumetric stress; validate affine and one-cell cases.
3. Add the minimal derived mixed solid model and verify mixed pressure, non-zero stress traces, SNES, MFFD, and domain-safe trials.
4. Add the stateless nominal epicardial and basal support BCs and verify force and traction measures.
5. Add viscosity using accepted old-time E/gradD data, with rejection and repeated-residual tests.
6. Configure and complete Case B, including pressure curve, time step, follower sign, support parameters, and restart.
7. Reuse electroMechanicalLaw unchanged with a pinned Ta table for a one-cell and then full-mesh activation smoke test.
8. Add ArosticaActiveTension, accepted-time commit/update semantics, restart tests, and the full Case A run.

## Final audit summary

### Existing components reusable unchanged

- mechanicalLaw runtime-selection, deformation-gradient, tangent, and modulus interfaces.
- mechanicalModel law selection, cell/face correction, and submesh support.
- electroMechanicalLaw active stress push-forward and field-based Ta path.
- PETSc/SNES, mixed vector packing, residual/JFNK trial handling, pressure stabilisation, and current-area follower-pressure geometry in the base solid model.
- solidTraction for prescribed pressure after an explicit patch-orientation and force-measure validation.
- Existing Guccione/HGO/Gultekin classes as unchanged structural references.

### Minimal derived solid model required: yes

Exact reason: the default mixed split removes dev from the passive stress before adding -pI, but Aróstica's I4f, I4s, and I8fs derivatives can have non-zero Cauchy traces. The benchmark requires those traces plus the independent pressure stress. The current protected hook already provides the correct sigma_law - pI branch, and the protected runtime-name constructor allows a small derived runtime class. Copying the base solid-model source is unnecessary and would risk divergence.

### New source classes proposed

- ArosticaHolzapfelOgdenElastic (Milestone 1 passive law).
- Optional ArosticaHolzapfelOgdenViscoelastic only if keeping viscosity in a separately selectable class is preferred; otherwise add viscosity to the same law in Milestone 2.
- ArosticaNonLinGeomTotalLagTotalDispSolid.
- ArosticaActiveTension plus its accepted-time/restart provider interface.
- arosticaNormalSpringDashpotTraction.
- arosticaVectorSpringDashpotTraction.

### Existing classes that must remain unchanged

- nonLinGeomTotalLagTotalDispSolid base implementation.
- electroMechanicalLaw.
- GuccioneElastic.
- HolzapfelGasserOgdenElastic.
- GultekinTwoFibreElastic.
- mechanicalLaw and mechanicalModel.
- solidTraction and existing spring/contact classes.

### Blockers and uncertainties

- The pinned toolkit is not a mechanics reference implementation; the local Finsberg code is the constitutive/activation cross-check and should be pinned as an explicit external reference in implementation tests.
- The paper's exact chi(x)=x switch is non-smooth at I4=1, while the published/reference code defaults to a logistic approximation and the accepted case currently requests logistic. The implementation must choose and report which benchmark convention is being reproduced.
- The current accepted case dictionaries already name absent Aróstica runtime laws/BCs and an ArosticaJRelation setting. They are scaffolding and must not be mistaken for compiled source.
- The exact nominal-to-current traction conversion for the support BCs must be tested against the base model's tractionBoundarySnGrad and patch normal orientation. This is the main remaining boundary-interface design issue.
- The source has old-time infrastructure, but the safest viscosity design still needs a deliberate accepted E/face-E ownership decision before implementation.
- The current case uses backward OpenFOAM time differencing, whereas the paper benchmark table reports several teams' time-integrator choices. The selected scheme and time step must be declared as part of the reproduction result.

