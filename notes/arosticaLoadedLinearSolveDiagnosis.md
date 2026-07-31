# Aróstica loaded two-patch linear-solve diagnosis

## Decision

**PHASE 2B SOURCE ACCEPTED; CASE/SOLVER INTEGRATION NOT YET ACCEPTED.**

The constitutive law, accepted-history lifecycle, trial-state boundary
conditions, and restart path are not the cause of the loaded failure. The
scratch loaded case is still not an accepted production integration case.
The accepted benchmark directory was not modified.

Phase 2C diagnostics are recorded in
`notes/arosticaJfnkPreconditionerDevelopment.md`. They confirm that the
original `-snes_mf` run had no useful preconditioner, while the corrected
`-snes_mf_operator` arrangement still has a large action mismatch between
the complete MFFD operator and the compact assembled matrix.

## Exact failing case

The original failure was preserved at:

```text
/tmp/ArosticaTwoPatchLoad.MBTzKK/log.solids4Foam.traction
```

The case is a serial, one-MPI-process copy of the monoventricle `mesh1`
topology with 17,625 tetrahedral cells and patches endocardium (2,500
faces), epicardium (3,474 faces), and base (256 faces). The command was run
from the case directory:

```text
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc
solids4Foam -case . > log.solids4Foam.traction 2>&1
```

The first non-zero step is `Time = 0.001`, with `deltaT = 0.001`. The
endocardial traction is `(0.1 0 0)` and endocardial pressure is zero. The
epicardial support is normal spring/dashpot with `alpha = 1e8 Pa/m`,
`beta = 5e3 Pa s/m`; the basal support is full-vector with `alpha = 1e5
Pa/m`, `beta = 5e3 Pa s/m`. The material is
`ArosticaHolzapfelOgdenViscoelastic`, `eta = 100 Pa s`, bulk modulus `1e6
Pa`, logistic compression switch `K = 100`, and passive representative
`impK = 275 Pa`. `solvePressure true`, legacy momentum stabilisation, and
the default `faceStressTreatment interpolatedCell` were used.

The failing `petscOptions` contained `-snes_mf` but did not contain
`-snes_mf_operator`; it selected `newtonls` with basic line search, LGMRES,
1000 maximum KSP iterations, and a fieldsplit Schur configuration. The
exact options and original logs are retained in the scratch case.

Repository state at the start of diagnosis:

```text
branch: solid-pressure-land3Fix
commit: ebf998ff25945f611a44e5376656690c62768d32
```

The worktree already contained the Phase 2A/2B and unrelated pre-existing
changes. No accepted benchmark files were changed by this diagnosis.

## PETSc reason hierarchy

With diagnostic PETSc options in
`/tmp/ArosticaLoadedDiagnosis.H3aS4i/petscOptions.diag`, the first residual
norm was `3.248312922653e-05`. The exact hierarchy was:

```text
SNES: DIVERGED_LINEAR_SOLVE, nonlinear iteration 0
KSP : DIVERGED_ITS, 1000 iterations
PC  : no PC failure; PC type none
```

The final true KSP residual was `1.217666807445e-06`, or
`3.748613007552e-02` of the initial norm. The KSP was LGMRES with a
70,500-by-70,500 matrix-free operator. This is not
`DIVERGED_NANORINF`, breakdown, zero-pivot, or a failed AMG/factorisation
setup.

The `-snes_mf_operator` diagnostic variant built the intended fieldsplit
preconditioner. Both `fieldsplit_0` and `fieldsplit_1` reported one-step
`CONVERGED_ITS` applications, while the outer LGMRES again terminated in
`DIVERGED_ITS` after 1000 iterations. Its final true residual ratio was
approximately `3.82295e-02`. Successful sub-block applications therefore
do not mean that the preconditioner approximates the matrix-free physical
operator well enough.

The fresh `-snes_view` output removes an important ambiguity. In the
original failing run, PETSc printed:

```text
SNES: type: newtonls
Jacobian is applied matrix-free with differencing, no explicit Jacobian
PC: type: none
linear system matrix = precond matrix:
  type: mffd
```

Thus the original `-snes_mf` run was configuration 2: `A=MFFD`, with no
useful `P`. It was not configuration 1 (`A=P=assembled approximate`). The
helper supplies the same `A_.m` as both `jac` and `B` to `SNESSetJacobian`;
PETSc's `-snes_mf` option replaces the applied matrix path with MFFD. With
`-snes_mf_operator`, the view instead printed `linear system matrix followed
by preconditioner matrix`, with `type: mffd` followed by `type: mpiaij`.
That is the intended configuration 3: `A=MFFD`, `P=assembled approximate`.

Configuration 3 was tested separately and still reached outer
`DIVERGED_ITS`; therefore the missing option is a real scratch-case
configuration defect, but it is not by itself the explanation for the
intended JFNK arrangement's failure.

## Residual admissibility

The failing solve is the first Newton linear solve from the accepted zero
state. At that state:

```text
min(J) = max(J) = 1
F = I
D = 0
p = 0
E = 0 and Edot = 0
sigmaViscous = 0
spring traction = 0
dashpot traction = 0
applied endocardial traction magnitude = 0.1
```

The pressure equation is zero at `J=1, p=0`; the reported combined PETSc
residual norm is finite. Repeated residual evaluations of the same state
were deterministic. There was no exponential overflow, NaN/Inf field, or
non-positive `J` before the failed KSP solve. In the direct assembled basic
Newton experiment, a later unrestricted step did reach `J = -5.41919578363`;
that is a consequence of the poor nonlinear direction, not the state at the
original KSP failure. A safeguarded backtracking experiment rejected the
step and ended in `DIVERGED_LINE_SEARCH` without an invalid constitutive
state.

The current production log reports the combined SNES norm rather than
separate displacement/pressure block norms immediately before the KSP
call. The source audit of `formResidual()` shows that the pressure block is
zero at the initial state and the finite combined norm is therefore a
finite momentum-load residual.

## Direct solve discrimination

The assembled test was run from `/tmp/ArosticaLoadedDirectSolve` with
`-snes_mf` removed, `-ksp_type preonly`, and `-pc_type lu`. PETSc/MUMPS
reported a 70,500-by-70,500 assembled matrix, zero null pivots, and a
one-iteration converged linear solve. The full-load direct solve therefore
does not report singularity or factorisation failure. It subsequently takes
a non-admissible unrestricted Newton step. With backtracking, the residual
decreases only from `3.248312922653e-05` to approximately
`3.248238079145e-05` and the line search terminates in
`DIVERGED_LINE_SEARCH`.

## Component isolation and load sweep

Valid isolation runs from the same scratch configuration produced:

| Variant | Result |
|---|---|
| Elastic law, both supports, non-zero traction | `DIVERGED_LINEAR_SOLVE` |
| Viscoelastic law, `eta=0`, both supports | `DIVERGED_LINEAR_SOLVE` |
| Viscoelastic law, `eta=100`, both supports | `DIVERGED_LINEAR_SOLVE` |
| Epicardial support only | `DIVERGED_LINEAR_SOLVE` |
| Basal support only | `DIVERGED_LINEAR_SOLVE` |
| Both supports | `DIVERGED_LINEAR_SOLVE` |
| Small non-zero pressure | `DIVERGED_LINEAR_SOLVE` |
| Original non-zero traction | `DIVERGED_LINEAR_SOLVE` |
| `interpolatedCell` face stress | `DIVERGED_LINEAR_SOLVE` |
| `directConstitutive` face stress | `DIVERGED_LINEAR_SOLVE` |

`solvePressure false` is rejected by the Aróstica solid model and is not a
valid isolation. A case with neither support would introduce physical
rigid-body freedom and was not used as evidence against the law.

The smallest valid trigger is a non-zero load in the mixed PETSc
configuration, not eta, either support, pressure versus traction, or face
stress treatment. The original matrix-free path fails even at load fraction
`1e-4` (traction magnitude `1e-5`); the zero-load case converges.

## Approximate-Jacobian audit

The production `formJacobian()` constructs the displacement block from

```text
momentumStabilisation().vectorJacobian(D, &impKf_)
- rho()*fvm::d2dt2(D)
- dampingCoeff()*rho()*fvm::ddt(D)   [when enabled]
```

The pressure block contains compressibility, pressure stabilisation, and
approximate grad/div coupling. The approximate displacement block does not
contain the current-trial derivative of passive material stress, viscous
material stress, epicardial spring/dashpot traction, or basal
spring/dashpot traction. Those terms are present in the physical residual
and matrix-free operator.

For the failing case the characteristic scales are:

| Contribution | Scale | Ratio to `impK=275` |
|---|---:|---:|
| passive representative `impK` | `2.75e2` | `1` |
| inertia `rho/deltaT^2` | `1.00e9` | `3.64e6` |
| viscosity `eta/deltaT` | `1.00e5` | `3.64e2` |
| epicardial spring | `1.00e8` | `3.64e5` |
| epicardial dashpot/deltaT | `5.00e6` | `1.82e4` |
| basal spring | `1.00e5` | `3.64e2` |
| basal dashpot/deltaT | `5.00e6` | `1.82e4` |

`pressureEqnScale=550`, the pressure unknown scale is `550`, and the bulk
modulus is `1e6`. `rAUf` is formed with the documented negative inverse of
the assembled momentum diagonal, giving a positive reciprocal scale for
pressure stabilisation. MUMPS reports no null pivots. Finite inertia and the
finite `1/kappa` pressure term also prevent exact rigid-translation and
constant-pressure null modes at this timestep; the issue is approximation
quality rather than an exact nullspace.

## Controlled solver experiments

| Experiment | Physical residual changed? | Result |
|---|---|---|
| assembled `preonly + LU` | no | linear solve succeeds; unrestricted Newton reaches negative `J` |
| assembled `preonly + LU` with backtracking | no | admissible, but `DIVERGED_LINE_SEARCH` |
| matrix-free operator plus assembled fieldsplit/hypre P | no | subsolves converge; outer `DIVERGED_ITS` at 1000 |
| matrix-free operator without `-snes_mf_operator` | no | `PC none`; outer `DIVERGED_ITS` at 1000 |
| FGMRES with assembled LU P | no | Krylov solves converge, but SNES `DIVERGED_DTOL` |

No random production PETSc option change was committed. Adding material,
viscous, and boundary tangents to the approximate displacement block is a
possible future generic preconditioning improvement, but it needs a broader
regression campaign.

## Root cause and remaining work

The evidence supports category **B: generic mixed-solver/preconditioner
limitation**, with a contributing scratch configuration issue: the original
options used `-snes_mf` without `-snes_mf_operator`, which intentionally
leaves the matrix-free KSP without the assembled preconditioner. Adding the
assembled preconditioner exposes, but does not remove, the larger mismatch
between the physical material/boundary tangent and the segregated
approximate Jacobian.

No production source correction was required. Loaded multi-timestep,
loaded process-restart, and loaded MFFD campaigns remain blocked because no
justified solver configuration produced a converged loaded first timestep.
The zero-load multi-timestep/restart and all constitutive, trial-state, and
MFFD source tests remain passing.

## Phase 2C pressure-to-momentum correction

The original compact pressure-to-momentum block was assembled with the
least-squares `fvm::grad(p)` insertion only.  The production residual uses
the total-Lagrangian nominal pressure stress

```text
Pnominal,p = -J*p*F^-T
```

and replaces the pressure traction on explicitly prescribed traction patches.
Consequently, the compact matrix must contain both the cell/face gradient
contribution and the corresponding pressure-face correction on those patches.
The corrected assembly uses the Gauss finite-volume pressure-gradient
insertion with the configured pressure-unknown scale, followed by the
reference/current face-area pressure traction correction for solid-traction
patches.  This changes only `P`; the physical residual and constitutive
equations are unchanged.

Before the correction the pressure-only production audit gave
`||Jv_D|| ~= 2.94e-1` and `||Pv_D|| ~= 4.86e-7`.  After the correction the
compact `P_Dp` block has:

| Quantity | Corrected value |
|---|---:|
| nonzero count | `242269` |
| Frobenius norm | `7.25793879554e-1` |
| maximum absolute coefficient | `8.59725447703e-3` |
| minimum non-zero coefficient | `4.23027752967e-9` |

Constant, smooth, localised, checkerboard, and one-cell pressure modes match
the MFFD displacement action to approximately `1e-13` relative error, with
unit sign correlation.  The physical residual reference norm remains
`3.248312922653e-05` before and after the compact-matrix correction.

## Phase 2C controlled displacement-block experiments

Two optional compact-preconditioner switches were added, both defaulting to
`false` so existing production behaviour is preserved:

```text
preconditionerMaterialTangent  false;
preconditionerBoundaryTangent  false;
```

With corrected `P_Dp`, exact LU applied to the assembled compact `P` still
reached KSP iteration 100 and diverged in the outer MFFD solve.  The true
residual ratios at iteration 100 were:

| Compact-P variant | KSP result | `||r_100||/||r_0||` |
|---|---|---:|
| corrected `P_Dp` only | `DIVERGED_ITS` | `0.56763` |
| plus passive material tangent | `DIVERGED_ITS` | `0.44712` |
| plus boundary spring/dashpot tangent | `DIVERGED_ITS` | `0.48667` |
| plus passive and boundary tangents | `DIVERGED_ITS` | `0.34902` |

The combined audit confirms that the pressure coupling is corrected: smooth
displacement and coupled-mode action errors are approximately `4.3e-16` and
`5.1e-16`, respectively.  The remaining displacement mismatch is dominated
by geometric/material terms: approximate rotation error `0.9403` and
localised displacement error `0.8036`.  Thus the next required work is a
broader, validated `P_DD` approximation, not another pressure-block fix.

The current classification is therefore:

```text
PRESSURE COUPLING CORRECTED; CONTINUE P_DD DEVELOPMENT
```

## Phase 2C Stage 1 displacement finite-difference reference

`tests/ArosticaViscoelastic/preconditionerAudit/ArosticaJvPvAudit.C` now
constructs central displacement-only finite differences at fixed pressure and
accepted history. It reports separate `J_DD^FD` and `J_pD^FD` actions,
compares them with the compact P action, and tracks the four explicit
accepted-history fields for every perturbation.

The utility was rebuilt with OpenFOAM v2312 and run on separate scratch
copies for both `interpolatedCell` and `directConstitutive` face stress. The
epsilon sweep was `1e-4, 1e-5, 1e-6, 1e-7, 1e-8`; the stable range for the
localised mode is `1e-6` through `1e-8`. Both runs completed successfully.

Common results were:

```text
baseline repeated-residual difference = 0
maximum accepted-history checksum change = 0
P_Dp nonzero count = 242269
P_Dp Frobenius norm = 0.725793879554
```

At stable epsilon, the interpolated-cell path has relative J/P errors of
approximately `0.0278` for rotation, `0.0249` for volumetric deformation,
`0.0323` for fibre deformation, `0.0257` for sheet deformation, and
`0.0211` for shear. The localised displacement remains the largest error at
`0.8036` with cosine `0.7573`.

For direct constitutive face stress, the corresponding errors are
approximately `0.0260`, `0.0268`, `0.0561`, `0.0670`, and `0.0364`, with
cosines above `0.997`. The direct localised mode remains the largest at
`0.8735`, although its cosine is `0.9779`; this indicates magnitude and
stencil under-representation rather than an opposite-sign action.

The corrected pressure-only action remains at approximately `1e-13` relative
error. No physical residual, constitutive equation, or accepted-history
behavior changed. Audit-only mode normalisation and smooth representative
fibre/sheet directions were required to keep the largest finite-difference
perturbations admissible on the direct-face path.

Stage 1 confirms that the next implementation target is term-by-term
`J_DD` decomposition, with particular attention to nominal stress/geometric
mapping and localised direct-face constitutive response.

## Phase 2C Stage 2 decomposition result

The Stage 2 audit was completed on copied scratch cases only. The production
residual, corrected `P_Dp`, PETSc A/P arrangement and accepted-history
lifecycle were unchanged. The audit decomposition closed at relative errors
`3.5881e-16` (`interpolatedCell`) and `3.0292e-16`
(`directConstitutive`). Repeated residual and accepted-history differences
were zero.

At `epsilon = 1e-6`, the localised interpolated-cell action is dominated by
the current viscoelastic face response (`1712.11`) and an oppositely directed
face interpolation/reconstruction action (`1208.05`, cosine `-0.5421`). The
passive action is only `58.25`; inertia is `12.985` and stabilisation is
`0.955`. The corresponding direct-face decomposition has the same passive,
viscous, inertia and stabilisation actions, but the interpolation term is
zero to `2.8e-17`. This explains the direct/interpolated Stage 1 cosines of
approximately `0.9779` and `0.7573`: the paths differ by a face stencil and
the interpolated route introduces cancellation, not merely a uniform scale.

The zero-load accepted state has negligible passive stress, so the audit
separates the localised passive response as constitutive rather than
geometric: determinant/inverse-transpose contributions are below `8.1e-15`
and the passive material/geometric cross remainder is `1.12e-4` in the
direct-face force audit. A non-zero prestressed loaded state is still needed
before making a production geometric-tangent decision.

The current optional passive tangent is useful but incomplete: it changes the
localised interpolated `P_DD` action from norm `29.8606` (switches disabled) to
`229.576` (passive enabled), but leaves relative error `0.8036`. The boundary
tangent is not active for the interior localised perturbation. No optional
contribution was promoted to a production default in Stage 2.

## Viscous `P_DD` follow-up

The bounded direct-face viscous preconditioner experiment used the copied
case `/private/tmp/ArosticaTwoPatchLoadViscousDirect`, with
`faceStressTreatment directConstitutive`, non-zero `eta=100`, corrected
`P_Dp`, and `preconditionerViscousTangent true`. The first non-zero step was
`t=0.001`, `deltaT=0.001`; the residual remained finite and the direct tangent
reported 32,135 internal faces.

With the production fieldsplit/Hypre configuration, the outer Krylov solve
still ended as `DIVERGED_ITS` at 1,000 iterations. This is an improvement
targeted at `P_DD`, not a change to the physical residual. The audit-only
direct localised action improved from relative error `0.873543888734` without
the viscous term to `0.170972333645` with viscous only; cosine improved from
`0.977877910708` to `0.990991433726`. Combining the pre-existing passive and
boundary diagnostic blocks gave error `0.250961909622`, cosine
`0.991802171161`, because those compact additions are not exact for this
interior localised mode.

The exact diagnostic arrangement was verified by PETSc view output as
`linear system matrix: mffd` followed by `preconditioner matrix: mpiaij`.
`-ksp_type preonly -pc_type lu` selected MUMPS on the 70,500-unknown serial
case. The first three SNES linear solves converged in one iteration and
MUMPS reported zero null pivots, but the repeated factorisation crashed with
an illegal-instruction/memory-corruption signal at SNES iteration 3. The
exact-LU result is therefore diagnostic/backend-limited, not a successful
loaded solve.

The opt-in term is explicitly rejected for `interpolatedCell`. Its separate
interpolation/reconstruction mismatch remains unaddressed. Physical
residual checks, repeated residual checks, and accepted-history checks stayed
unchanged; the source default remains `preconditionerViscousTangent false`.

## Passive nominal-force follow-up

`preconditionerPassiveNominalTangent` is a new default-off direct-face-only
compact-preconditioner switch. The law returns the full nine-component
passive nominal-force derivative of `Ppassive = F&S`; the solid model adds it
in `formJacobian()` using the existing owner-neighbour stencil. No physical
residual or accepted history is changed, and the older scalar passive switch
cannot be combined with it.

The genuinely prestressed elastic smoke reference passed with nominal-force
magnitude `11960.7` and directional FD error `7.30483e-11`. The copied PETSc
audit remained deterministic and closed its term decomposition to
`1.75081968491e-16`. Its direct localised stable-epsilon action was
`||J_DD v||=1769.3965`, `||P_DD v||=1989.8468`, relative error `0.1854693`,
cosine `0.9916076`. The audit state was kept tiny because larger affine
prestress values were incompatible with the ventricular traction-boundary
scratch field and correctly triggered the existing positive-`J` guard. Thus
this stage does not claim geometric-stiffness or loaded-solver acceptance.

The elastic runtime bus error was a stale executable after the mechanicalLaw
virtual-interface change, not a constitutive failure. Rebuilding the smoke
executable passed runtime selection and the prestressed nominal-force test;
the Alltest script now rebuilds first.

The bounded loaded follow-up used a fresh copy with only the corrected
pressure block, the direct-face viscous tangent and the new passive nominal
tangent enabled. The exact diagnostic setup was
`-snes_mf_operator`, `-ksp_type preonly`, `-pc_type lu`; PETSc displayed
`mffd` as the linear-system matrix and `mpiaij` as the preconditioning matrix.
MUMPS factored `P` with zero null pivots and the first three linear solves
returned `CONVERGED_ITS` in one iteration, but the run aborted at the next
factorisation with signal 4 (illegal instruction/memory corruption). The
production fieldsplit copy stopped at the first loaded time `t=0.001` with
`SNES DIVERGED_LINEAR_SOLVE`; no loaded timestep was accepted. These results
are diagnostic/backend and coupled-preconditioner failures, not evidence of
a changed physical residual.

The runtime smoke harnesses were repaired independently: each now rebuilds
after the `mechanicalLaw` virtual interface change, and the viscoelastic and
boundary scripts use portable `sed -i.bak` commands on macOS. Elastic,
viscoelastic, boundary, mixed-solid and Gultekin regression checks pass after
the repair.

## Coupled-block follow-up

The later bounded mixed-block investigation used the retained mesh1
synthetic tiny-load case in scratch copies only. The physical stabilisation
settings were preserved as `scaleFactor = 10` and `scaleFactorJacobian = 10`;
no production case used a zero scale. Temporary zero-scale copies were
decomposition controls only and are not candidate configurations.

The current assembled blocks are:

```text
P_DD: 737055 nonzeros, Frobenius norm 3.88677198695e5
P_Dp: 242269 nonzeros, Frobenius norm 7.25793879554e-1
P_pD: 242269 nonzeros, Frobenius norm 7.18883915619e-1
P_pp:  81895 nonzeros, Frobenius norm 4.58485262431e-3
```

`P_Dp` remains corrected and agrees with the complete MFFD
pressure-to-momentum action. The constant pressure mode is also correctly
scaled (`J_pp/P_pp` error approximately `1.1e-12`). Nonconstant pressure
modes expose a structural mismatch: smooth-gradient error `0.9240` with
cosine `0.4085`, checkerboard error `0.9303` with cosine `0.3948`, and
one-cell error `0.9476` with cosine `0.9341`.

The reason is that the Rhie-Chow residual differentiates a least-squares
gradient inside a face high-pass flux, whereas the compact
`diffStencilLaplacianStab::scalarJacobian()` supplies only an
`fvm::laplacian` neighbour stencil. The exact derivative is wider than the
current compact preconditioner stencil. This is not fixed by multiplying the
existing Jacobian by zero or by blindly changing its scalar multiplier.

Schur factorisation experiments did not converge the first loaded step at
the diagnostic iteration cap. A properly configured exact serial LU of the
assembled compact P had no zero pivot and retained `A = MFFD`, `P = MPIAIJ`,
but the true residual ratio was still approximately `1.25e-1` at 100 outer
iterations. Full/selfp exact block LU reached approximately `3.09e-1` at
100. FGMRES with the current fieldsplit reached approximately `5.29e-1` at
100; disabling PETSc mixed-field scaling worsened the ratio to approximately
`1.49`. These are diagnostic results, not loaded acceptance.

The current classification is therefore:

```text
PHASE 2B SOURCE ACCEPTED;
COUPLED COMPACT-P/SCHUR DEVELOPMENT STILL REQUIRED
```

### Follow-up pressure-block implementation

An opt-in `preconditionerLeastSquaresPressureCoupling` assembly now matches
the production least-squares gradient used by the pressure equation. Its
scratch localized `P_pD` finite-difference action error is below `1e-9`, with
cosine `1`, while the physical residual and accepted histories remain
unchanged. The default remains false pending a complete coupled-preconditioner
acceptance result.

An additional opt-in least-squares `P_pp` correction was tested with the
existing nonzero compact Laplacian retained. It did not close the nonconstant
pressure action mismatch and did not produce a loaded acceptance result, so it
remains diagnostic only. No test configuration used `scaleFactorJacobian 0` as
a fix: with nonzero physical stabilization that setting would indeed zero the
corresponding compact preconditioner contribution. Scratch loaded runs kept
`scaleFactor = 10` and `scaleFactorJacobian = 10`.

The next source-side study must either provide a wider-stencil or otherwise
consistent preconditioner for the least-squares Rhie-Chow derivative, or
redesign the mixed Schur approximation. The physical residual and accepted
history lifecycle remain unchanged.

## Full wide pressure-block source correction

The source-side correction now assembles the actual production
least-squares/Rhie-Chow pressure-stabilisation derivative behind the existing
opt-in switch `preconditionerLeastSquaresPressureStabilisation`. It is not a
replacement pressure residual and does not alter the accepted history path.
The wider stencil covers the least-squares gradient graph, face interpolation,
corrected `snGrad`, cell divergence, boundary terms, and processor-interface
columns. The finite bulk term, `pressureUnknownScale`, `pressureEqnScale`,
cell-volume weighting, and `scaleFactorJacobian=10` are applied once.

The pressure-only finite-difference audit passed in serial and in a genuine
two-rank decomposition. Constant pressure had relative error about
`2.3e-12`; smooth, localized, checkerboard, one-cell, random, and
boundary-adjacent modes were at machine precision with cosine 1. The source
library was rebuilt successfully and the audit utility was rebuilt; no
accepted Case B parent files were edited.

The next bounded source-side investigation is the standalone global
explicit-J-versus-P diagnosis described in
`notes/arosticaGlobalJVersusPDiagnosis.md`. It is deliberately separate from
the accepted-history and physical-residual lifecycle and begins on a
six-tetrahedron scratch case before transferring only identified action modes
to a temporary Case B copy.

The zero-load scratch regression passed for ten steps, with `D=0`, `p=0`,
`J=1`, and finite fields. A separate-process restart from `t=0.005` to
`t=0.010` reproduced the uninterrupted written fields byte-for-byte.
The fixed 0.016071182 Pa synthetic load remains blocked at the first
nonzero step for full+selfp, lower+selfp, and full+a11. Their true residuals
remain near the initial residual (approximately 1.00010, 1.00000, and
1.00258 ratios respectively), so the wide `P_pp` action is accepted but the
coupled solver integration is not.
