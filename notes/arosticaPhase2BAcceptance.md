# Aróstica Phase 2B acceptance campaign

## Decision

**PHASE 2B SOURCE ACCEPTED; CASE/SOLVER INTEGRATION NOT YET ACCEPTED.**

The critical constitutive, trial-state, history, MFFD, zero-load ventricular
topology, and process-restart checks pass. Acceptance is withheld because a
controlled non-zero loading run on the scratch monoventricle mesh reaches
the production PETSc solve but currently returns
`DIVERGED_LINEAR_SOLVE`. This campaign does not change the accepted benchmark
case and does not begin activation, pressure-history generation, or Case B.

## Recorded implementation

The repository was audited on branch `solid-pressure-land3Fix`, commit
`ebf998ff25945f611a44e5376656690c62768d32` before the campaign. The core
Phase 2B law is
`ArosticaHolzapfelOgdenViscoelastic`, derived from
`ArosticaHolzapfelOgdenElastic`; its runtime name and the elastic runtime
name are unchanged.

The history lifecycle is:

| Operation | Function | State effect |
|---|---|---|
| Construct/read | `ArosticaHolzapfelOgdenViscoelastic::ArosticaHolzapfelOgdenViscoelastic` | Constructs registered cell/face `E`, `EOld`, `EOldOld`; reads them if present. |
| Current trial | `updateViscousFields()` | Recomputes current cell/face `E`, `Edot`, `Sviscous`, and `sigmaViscous`; does not shift accepted fields. |
| Stress | `correct(volSymmTensorField&)`, `correct(surfaceSymmTensorField&)` | Adds viscous stress to the shared passive stress. |
| First/restart seed | `initialiseHistoryIfRequired()` | Seeds missing accepted history once from the current accepted strain. |
| Accepted update | `updateTotalFields()` | The only function allowed to shift accepted history. |
| Restart write policy | `setRestart()` | Sets history fields to `AUTO_WRITE`; it does not advance them. |
| Temporal coefficients | `temporalCoefficients()` | Reads the configured scheme and OpenFOAM v2312 `deltaT`/`deltaT0`. |

`mechanicalModel::updateTotalFields()` calls each law’s update, and
`solidModel::updateTotalFields()` is called by the solver after a completed
physical timestep. `updateTotalFields()` stores the last OpenFOAM
`timeIndex`; a second call in the same timestep is a no-op. Thus the measured
number of history shifts is exactly one per accepted timestep.

## History evidence

The one-cell production smoke recorded history fields before/after trial
calls, after a rejected C trial, after one commit, and after a second same-time
commit. Maximum differences were:

| Check | Cell | Direct face |
|---|---:|---:|
| First shift (`EOldOld <- EOld`, `EOld <- E`) | 0 | 0 |
| Second update at same time index | 0 | 0 |
| Any trial/rejected-trial history change | 0 | 0 |
| First-step Edot error | 0 | 0 |
| Five unequal-step backward Edot error | 0 | 0 |

The production PETSc residual harness also observed maximum accepted-history
change `0` after every residual and MFFD perturbation.

## PETSc trial and MFFD evidence

The production harness used mixed displacement-pressure unknowns,
`solvePressure true`, non-zero `eta`, logistic compression, and both
Aróstica spring-dashpot variants. The first residual after changing the PETSc
trial vector changed the patch displacement, velocity, traction, and residual.
Repeated-state maximum differences were `0`.

The MFFD-style sweep was:

| Patch | epsilons | reference norm | max difference from `1e-7` |
|---|---|---:|---:|
| Normal | `1e-7 ... 1e-2` | `7.20872e6` | `2.13393e-9` |
| Vector | `1e-7 ... 1e-2` | `2.27423e4` | `7.49989e-12` |

The six finite-difference norms were constant to the displayed precision in
both cases. Trial evaluation order did not affect the result, and histories
were unchanged.

## Restart evidence

The ventricular scratch case was copied from the monoventricle mesh topology
to `/tmp`; the accepted benchmark directory was not modified. An
uninterrupted process ran through `0.005`. A second process ran through
`0.002`, wrote, terminated, and restarted from `latestTime` through `0.005`.
At each of `0.003`, `0.004`, and `0.005`, the following were byte-identical:

`D`, `p`, `F`, direct-face `Ff`, `EOld`, `EOldOld`, face `E`, face `EOld`,
face `EOldOld`, `Edot`, `Sviscous`, `sigmaViscous`, and total stress. The
solver reported convergence at every zero-load timestep and `min(J)=1`.
Therefore the measured L2 and Linf differences are both `0` for every
compared field. The write/read path includes all explicit history fields;
they are not reconstructed from only `D`.

The restart run also exposed a real boundary-condition serialization defect.
The common spring-dashpot writer now emits a dimensioned coefficient name
token, which is required by the v2312 dimensioned-scalar reader. The boundary
runtime regression was rerun after this correction.

## Two-patch ventricular topology

The scratch case selected:

- `arosticaNormalSpringDashpotTraction` on `epicardium`;
- `arosticaVectorSpringDashpotTraction` on `base`;
- `ArosticaHolzapfelOgdenViscoelastic`, `eta=100`, logistic switch;
- `arosticaNonLinearGeometryTotalLagrangianTotalDisplacement`;
- `solvePressure true`; zero active stress and zero endocardial pressure.

The zero-load equilibrium run completed three timesteps, then the restart
campaign completed five. Both boundary types retained reference-area force
integration; the normal condition is scalar-normal and the basal condition is
full-vector. Fibre and `tbar` fields were copied unchanged. Written cell and
direct-face diagnostics were finite, and `min(J)=1`.

A controlled non-zero pressure/traction perturbation was also attempted on
the same scratch topology. It entered the PETSc production solve but returned
`DIVERGED_LINEAR_SOLVE` before a loaded timestep could be accepted. This is
the remaining critical blocker; no source workaround was introduced because
the failure is a solver/integration issue requiring a separate diagnosis.

## ddt cost

The source contains one `fvc::ddt(D)` per Aróstica patch update, hence two
full-volume constructions per residual with epicardium and base. The probe
measured:

| Case | Pair time for two constructions | Fraction |
|---|---:|---:|
| One-cell runtime case | below timer resolution | not measurable |
| 17,625-cell ventricular scratch | `0.0005 s` | about `0.64%` of the `0.078 s` residual estimate |

No shared cache was added. The current calculation remains trial-state
dependent and safe under order changes, MFFD perturbations, rejected trials,
and restart.

## Regression summary

Passing during this campaign:

- dependency-free Aróstica elastic, boundary, and viscoelastic tests;
- elastic runtime smoke;
- viscoelastic runtime smoke;
- PETSc boundary runtime smoke for normal and vector patches;
- mixed-solid runtime smoke;
- compiled Gultekin runtime regression;
- zero-load two-patch ventricular scratch run;
- five-step two-process ventricular restart;
- `wmake libso` for OpenFOAM v2312.

The NumPy-only Gultekin script remains environment-blocked because NumPy is
not available; it is not classified as pass or fail. The loaded two-patch
case is blocked by `DIVERGED_LINEAR_SOLVE` as described above.

## Remaining work before Case B

Diagnose and repair the loaded two-patch PETSc linear solve, then rerun the
non-zero loading and loaded restart/MFFD checks. After that, repeat the full
regression matrix and only then begin Case B integration. Activation,
activation ODEs, chamber-pressure histories, fibre generation, and changes
to the accepted benchmark remain out of scope.

## Loaded linear-solve diagnosis

The exact loaded-case investigation is recorded in
`notes/arosticaLoadedLinearSolveDiagnosis.md`. The first failing step is
`t=0.001`, `deltaT=0.001`, with endocardial traction `(0.1 0 0)`. The
residual is finite at `J=1`, `D=0`, and `p=0`. The exact PETSc hierarchy is
SNES `DIVERGED_LINEAR_SOLVE`, caused by outer LGMRES `DIVERGED_ITS` at 1000
iterations; there is no PC factorisation failure. The matrix-free case has
`PC type none` because its options contain `-snes_mf` without
`-snes_mf_operator`.

An assembled MUMPS LU test reports a 70,500-by-70,500 matrix with zero null
pivots and one-step linear convergence. It nevertheless gives a poor
nonlinear direction; backtracking remains admissible but ends in
`DIVERGED_LINE_SEARCH`. Adding the assembled fieldsplit/hypre preconditioner
causes its subsolves to converge while the outer LGMRES still ends in
`DIVERGED_ITS`. Elastic, eta=0, eta=100, either support alone, both supports,
pressure/traction, and either face-stress route all fail at the first
non-zero load. The zero-load case passes.

The approximate displacement Jacobian contains inertia and stabilisation but
omits passive, viscous, spring, and dashpot tangents. The resulting scale
separation explains the generic solver limitation. No production source or
accepted benchmark correction was made. The loaded multi-timestep and
loaded restart tests therefore remain outstanding before Case B.

## Phase 2C preconditioner development

The intended `-snes_mf_operator` arrangement was verified directly with
PETSc views (`A=MFFD`, `P=MPIAIJ`). Exact LU of the compact P and exact
fieldsplit block LU did not converge the loaded outer Krylov solve. The
production action audit measured approximately 0.98 relative `Jv/Pv`
errors for displacement and coupled modes, with a near-zero compact
pressure-to-displacement response for the pressure-only mode. These results
are documented in `notes/arosticaJfnkPreconditionerDevelopment.md`.

Phase 2C is therefore not complete. No physical residual or constitutive
source correction was made.
