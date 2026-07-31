# Aróstica Phase 2C JFNK preconditioner development

## Status

**PRECONDITIONER DEVELOPMENT NOT YET COMPLETE.**

The Phase 2B physical residual, constitutive law, accepted-history lifecycle,
restart path, and MFFD trial-state tests remain unchanged and passing. The
loaded two-patch case is not yet a robust production integration case.

The accepted benchmark directory was not modified.

## Frozen physical implementation

The repository state used for this campaign was:

```text
branch: solid-pressure-land3Fix
commit: ebf998ff25945f611a44e5376656690c62768d32
```

The worktree contained the existing Phase 1/2A/2B and unrelated user
changes. The Phase 2C diagnostic test added in this campaign is
`tests/ArosticaViscoelastic/preconditionerAudit/`. The physical residual was
kept unchanged; the compact preconditioner and its optional diagnostics were
updated.

The Python constitutive, boundary-condition, and viscoelastic tests all
passed. No production constitutive or residual source was changed during
this campaign.

## A/P configurations

The helper calls `SNESSetJacobian` with the compact assembled matrix as the
preconditioner matrix. The intended PETSc configuration is:

```text
A = MFFD derivative of the complete production residual
P = compact assembled MPIAIJ matrix
```

This requires `-snes_mf_operator`. The original failing scratch options
used `-snes_mf`, which gave `A=MFFD, P=none`. A separate baseline using
`-snes_mf_operator` confirmed `A=MFFD, P=MPIAIJ`; that arrangement still
failed, so the missing option was a real case defect but not the complete
cause.

The loaded first-step state is the serial 17,625-cell monoventricle copy
with both Aróstica spring-dashpot patches, `solvePressure true`, backward
time discretisation, `eta=100 Pa s`, `rho=1000 kg/m3`, bulk modulus `1e6 Pa`,
epicardial `alpha=1e8`, `beta=5e3`, basal `alpha=1e5`, `beta=5e3`, and
endocardial traction `(0.1 0 0)`.

Baseline results at the first non-zero load step:

| Configuration | Result |
|---|---|
| `-snes_mf`, no useful PC | KSP `DIVERGED_ITS` at 1000; true residual ratio `3.7486e-2` |
| `-snes_mf_operator`, current fieldsplit/hypre | KSP `DIVERGED_ITS` at 1000; true residual ratio `3.8230e-2` |
| compact matrix as both A and P | KSP `DIVERGED_ITS` at 1000; true residual ratio `4.5845e-1` |
| MFFD A + exact LU(P) | KSP `DIVERGED_ITS` at 100; true residual ratio `9.9920e-1` |
| MFFD A + exact block LU, lower/a11 | KSP `DIVERGED_ITS`; at 100 iterations true ratio about `5.25e-1` |
| MFFD A + exact block LU, upper | KSP `DIVERGED_ITS`; at 20 iterations true ratio `7.42e-1` |
| MFFD A + exact block LU, full | KSP `DIVERGED_ITS`; at 20 iterations true ratio `9.32e-1` |
| MFFD A + exact block LU, diagonal | KSP `DIVERGED_ITS`; at 20 iterations true ratio `7.30e-1` |
| MFFD A + unlagged current fieldsplit | KSP `DIVERGED_ITS`; at 100 iterations true ratio `5.26e-1` |

The exact LU(P) view showed an MFFD operator followed by an MPIAIJ
preconditioner matrix and zero MUMPS null pivots. Thus factorisation failure,
preconditioner lagging, and fieldsplit subsolver singularity are not the
dominant causes.

## Jv/Pv action audit

`ArosticaJvPvAudit` forms the production residual difference quotient and
compares it with the compact matrix action at the same trial state. At
`epsilon=1e-7` the loaded first-step results were:

| Mode | `||Jv||2` | `||Pv||2` | `||Jv-Pv||/||Jv||` | cosine |
|---|---:|---:|---:|---:|
| smooth displacement | `2.5771e4` | `2.0648e3` | `9.8379e-1` | `2.4071e-1` |
| pressure only | `2.9350e-1` | `4.8625e-7` | `1.0000` | `1.6567e-6` |
| coupled | `1.9590e4` | `1.8607e3` | `9.8265e-1` | `2.2855e-1` |
| rigid translation | `8.9533e3` | `1.6074e3` | `9.7254e-1` | `2.4062e-1` |
| approximate rotation | `3.2438e5` | `4.2474e3` | `9.9665e-1` | `2.6199e-1` |
| localised displacement | `8.7034e2` | `3.2711e1` | `9.6979e-1` | `8.1040e-1` |

The pressure-only result is blockwise diagnostic: the compact matrix
contains the pressure-row response (`4.8625e-7`) but essentially no
pressure-to-displacement response, while the complete residual derivative
contains displacement components approximately `(0.0841, 0.1989, 0.1988)`.
The displacement and coupled modes show the same order-of-magnitude
under-representation. This is a physical-Jacobian approximation deficiency,
not a history mutation or MFFD ordering problem.

## Current compact matrix audit

The displacement block currently contains momentum stabilisation and
inertia/damping terms. The mixed pressure blocks contain compressibility,
pressure stabilisation, and approximate grad/div coupling. The compact
matrix does not include the current-trial tangent of passive material stress,
viscous stress, or spring-dashpot boundary forces. The complete MFFD
operator contains all of those terms.

The characteristic omitted scales for the failing case are approximately:

```text
passive representative impK       2.75e2
inertia rho/deltaT^2              1.00e9
viscosity eta/deltaT              1.00e5
epicardial spring                 1.00e8
epicardial dashpot/deltaT         5.00e6
basal spring                      1.00e5
basal dashpot/deltaT              5.00e6
```

Changing only `impKcoeff` to 10, 30, or 100 in scratch cases did not change
the measured outer Krylov history through 300 iterations. A scalar
representative stiffness is therefore not a sufficient preconditioner
improvement.

## Interpretation

The loaded failure is classified as a **generic mixed-JFNK/preconditioner
deficiency in the current compact P, aggravated in the original scratch
case by omission of `-snes_mf_operator`**. It is not evidence of a Phase 2B
physical residual defect. The next justified development stages are a
block-consistent pressure coupling audit and controlled addition of passive,
viscous, and reference-area spring/dashpot tangents to P, with each
contribution independently switchable and verified against the Jv/Pv audit.

No such production change has yet been committed because the current
action audit does not identify a safe single scalar replacement for the
missing block physics.

## Pressure-to-momentum correction

The PETSc solution is interlaced by cell as `[D_x,D_y,D_z,pHat]`, with block
size 4 and pressure column offset 3. The physical pressure is
`p = pressureUnknownScale*pHat`; the momentum row is not pressure-equation
scaled. The compact matrix is therefore

```text
P = [ P_DD  P_Dp ]
    [ P_pD  P_pp ]
```

The pressure part of the total-Lagrangian nominal stress is

```text
delta Pnominal_p = -J*delta(p)*F^-T
```

The production finite-volume momentum residual uses the current face force
`-p*SfCurrent`, where `SfCurrent = (J*F^-T) & Sf`. The old compact insertion
used the least-squares gradient helper and omitted the traction-boundary face
term. Consequently, a pressure-only perturbation had
`||Jv_D|| = 2.93504809324e-1` but `||Pv_D||` only
`4.86247941867e-7`.

The correction in
`nonLinGeomTotalLagTotalDispSolid::formJacobian()` uses the Gauss finite-volume
pressure-gradient stencil with `pressureUnknownScale_`, then adds the
opposite `pressureUnknownScale_*SfCurrent` coefficient on each
`solidTraction` boundary face. This cancels the pressure face force that the
production residual replaces with prescribed/trial traction. It changes P
only; `formResidual()` and the physical residual value are untouched.

The corrected P_Dp statistics on the 17,625-cell scratch mesh are:

```text
nonzero count       242269
Frobenius norm      0.725793879554
maximum coefficient 0.00859725447703
minimum nonzero     4.23027752967e-09
```

Constant, smooth-gradient, localised, checkerboard, and one-cell pressure
modes now have displacement actions matching MFFD to approximately
`1e-13` relative error and cosine 1. The physical residual initial norm is
unchanged at `3.248312922653e-05`.

## Controlled P_DD stages

Two optional compact-P-only switches were added for controlled development:

```text
preconditionerMaterialTangent false;
preconditionerBoundaryTangent false;
```

The first reuses the Aróstica passive `mat66` tangent through the new
`mechanicalLaw::passiveMaterialTangentField()` interface and forms a positive
row-sum scalar face bound for the existing compact Laplacian. It is disabled
by default and does not alter stress or residual evaluation.

The second adds the reference-area spring/dashpot displacement tangent as a
block-local P_DD contribution. For the normal patch it uses
`-A*(alpha + beta*c0)*N outer N`; for the vector patch it uses
`-A*(alpha + beta*c0)*I`. The coefficient `c0` is `1/deltaT` for Euler and
the current coefficient of OpenFOAM backward ddt, including the first-step
fallback and variable-deltaT BDF2 coefficient. Unsupported ddt schemes fail
clearly. The scratch first-step coefficient was `c0=1000 s^-1`, with tangent
coefficient range `5.1e6` to `1.05e8` and 3,730 patch faces.

These stages were measured independently with matrix-free A and exact LU(P):

| P additions | KSP result at 100 iterations | true residual ratio |
|---|---|---:|
| corrected P_Dp only | `DIVERGED_ITS` | `5.6763e-1` |
| passive material tangent | `DIVERGED_ITS` | `4.4712e-1` |
| boundary spring/dashpot tangent | `DIVERGED_ITS` | `4.8667e-1` |
| passive + boundary tangents | `DIVERGED_ITS` | `3.4902e-1` |

The combined Jv/Pv audit remains exact for pressure-only, smooth
displacement, coupled, and rigid-translation modes. The remaining relative
errors are approximately `0.9403` for the deterministic rotation mode and
`0.8036` for the localised displacement mode. Thus these contributions are
useful, justified diagnostics but are not yet a robust loaded preconditioner.

No loaded multi-timestep or loaded process-restart acceptance is claimed
after these P changes. The next required stage is a block-consistent
material/geometric P_DD approximation, followed by the loaded campaign.

## Stage 1 reduced finite-difference J_DD reference

The preconditioner audit was extended without changing the production
residual or adding another production `P_DD` contribution. The audit now
constructs central differences

```text
J_DD^FD(v) = [R_D(D + epsilon*v,p) - R_D(D - epsilon*v,p)]/(2*epsilon)
```

at fixed pressure and accepted material history. Deterministic modes include
smooth, localised, rigid-translation, coordinate-rotation, volumetric,
fibre, sheet, and fibre-sheet shear displacement modes, in addition to the
existing pressure and coupled modes. The epsilon sweep is
`1e-4, 1e-5, 1e-6, 1e-7, 1e-8`. The audit reports displacement and pressure
residual components separately and tracks `EOld`, `EOldOld`, `EfOld`, and
`EfOldOld` for every perturbation.

The utility was rebuilt successfully with OpenFOAM v2312 and run on separate
copies of the 17,625-cell two-patch scratch case:

```text
/tmp/ArosticaPhase2CAuditCombined
/tmp/ArosticaPhase2CAuditDirect
```

Both runs completed with `PASS`. The corrected pressure block remained at
242269 nonzeros and Frobenius norm `0.725793879554`. The repeated baseline
residual difference was `0`, and the maximum accepted-history checksum
change was `0` for both face-stress treatments.

The stable central-difference range was `1e-6` to `1e-8` for the localised
mode; the larger perturbations show finite-amplitude variation. At stable
epsilon, representative relative J/P errors and cosines were:

| Mode | Interpolated-cell | Direct-constitutive |
|---|---:|---:|
| smooth displacement | `~4e-16 / 1` | `~4e-16 / 1` |
| coupled displacement-pressure | `~5e-16 / 1` | `1.20e-5 / 0.99999999993` |
| rigid translation | `~8e-16 / 1` | `~8e-16 / 1` |
| rotation | `0.02780 / 0.999614` | `0.02597 / 0.999663` |
| localised displacement | `0.80361 / 0.75732` | `0.87354 / 0.97788` |
| volumetric displacement | `0.02488 / 0.999787` | `0.02684 / 0.999736` |
| fibre deformation | `0.03234 / 0.999615` | `0.05611 / 0.998467` |
| sheet deformation | `0.02572 / 0.999681` | `0.06701 / 0.997772` |
| fibre-sheet shear | `0.02105 / 0.999782` | `0.03640 / 0.999337` |

The pressure-only displacement action continues to agree with corrected
`P_Dp` to approximately `1e-13`. No material, geometric, viscous, or
boundary tangent was added by this stage. The next stage is term-by-term
decomposition of `J_DD`, beginning with nominal-force/geometric mapping and
the localised direct-face constitutive response.

## Stage 2 term-by-term `J_DD` decomposition

Stage 2 added an audit-only
`nonLinGeomTotalLagTotalDispSolid::momentumResidualDecomposition()` path and
the corresponding diagnostics in
`tests/ArosticaViscoelastic/preconditionerAudit/ArosticaJvPvAudit.C`. It
reuses the current production force path, including the selected face-stress
route, traction-boundary replacement, momentum stabilisation, inertia and
damping. It does not assemble or modify `P`, and it does not change the
physical residual or accepted histories.

The audit used `epsilon = 1e-6`, inside the Stage 1 stable range, on copied
17,625-cell cases:

```text
/tmp/ArosticaPhase2CAuditCombined   (interpolatedCell)
/tmp/ArosticaPhase2CAuditDirect     (directConstitutive)
```

The term sum reconstructs the complete displacement residual finite
difference to `3.5881e-16` relative error for the interpolated path and
`3.0292e-16` for the direct path. Repeated residual difference and
accepted-history checksum change remain exactly zero.

| Mode/path | Term | Action norm | Fraction of total | Cosine with total | P representation / mismatch |
|---|---|---:|---:|---:|---|
| localised/interpolated | passive | 58.25 | 0.0733 | 0.6165 | scalar diagnostic tangent; combined `P_DD` error `0.8036`, cosine `0.7573` |
| localised/interpolated | viscous | 1712.11 | 2.1553 | 0.8189 | omitted from current compact P |
| localised/interpolated | face interpolation | 1208.05 | 1.5208 | -0.5421 | omitted; cancels part of direct material action |
| localised/interpolated | stabilisation | 0.955 | 0.00120 | 0.5755 | existing compact path |
| localised/interpolated | inertia | 12.985 | 0.01635 | 0.8228 | existing compact path |
| localised/direct | passive | 58.25 | 0.0329 | 0.7759 | scalar diagnostic tangent; combined `P_DD` error `0.8735`, cosine `0.9779` |
| localised/direct | viscous | 1712.11 | 0.9676 | 0.9998 | omitted from current compact P |
| localised/direct | face interpolation | `2.8e-17` | `1.6e-20` | — | correctly zero for direct route |
| localised/direct | stabilisation | 0.955 | 0.000540 | 0.9139 | existing compact path |
| localised/direct | inertia | 12.985 | 0.00734 | 0.8961 | existing compact path |

The smooth mode is boundary-dominated: epicardial spring `141.526` and
dashpot `7.0763`, with exact compact action. Rotation, volumetric and shear
remain close in direction; their direct-path relative errors are `0.02597`,
`0.02298`, and `0.03640`. The current scratch state has zero accepted
stress and zero pressure, so the separated determinant and
inverse-transpose geometric terms are below `8.1e-15` for the localised
mode; the passive stress derivative is therefore the measurable internal
response in this state.

The direct/interpolated difference is structural rather than a scalar scale
factor. Direct constitutive face evaluation has no interpolation action and
the localised best-fit scalar for current `P_DD` is `7.53674`, leaving
`0.20918` relative error. The interpolated route has an oppositely correlated
face-path action (`-0.5421`), giving best-fit scalar `2.62041` but still
`0.65304` relative error. The largest localised `Jv-Pv` action is at cell 0;
the largest decomposed face action is internal face 2, owner cell 0 and
neighbour cell 833.

Isolated copied cases classify the optional P contributions. With both
switches disabled, localised interpolated `P_DD` has norm `29.8606`, relative
error `0.96979`, and cosine `0.81043`. The passive switch alone restores
`P_DD` norm `229.576`, error `0.80361`, and cosine `0.75732`. The boundary
switch alone has the no-tangent result for this interior perturbation because
its support is on traction patches. The scalar passive tangent is therefore
structurally useful but incomplete; the boundary tangent is irrelevant to
this mode and remains diagnostic-only.

The next single production contribution should be a current-trial viscous
`P_DD` approximation using the existing viscoelastic tangent and exact
backward coefficient. Its source scope should be limited to the optional
compact assembly in `nonLinGeomTotalLagTotalDispSolid::formJacobian()` and
the existing material-tangent infrastructure; no residual or constitutive
equation should be changed.

## Phase 2C viscous `P_DD` stage

The current-trial viscous nominal-force tangent was implemented as an opt-in
compact-preconditioner contribution. `ArosticaHolzapfelOgdenViscoelastic`
provides the direct-face derivative of `Pviscous = F*Sviscous`; the solid
model inserts its full tensor coupling on the internal owner/neighbour
stencil during `formJacobian()`. The temporal coefficient is obtained from
the law's existing `temporalCoefficients()` implementation, so startup Euler
fallback, backward/BDF2 and variable `deltaT` use the same `c0` as the
physical law. No history field is written or advanced by this path.

The switch is:

```text
preconditionerViscousTangent false;  // default
```

Requesting it with `faceStressTreatment interpolatedCell` fails clearly; no
direct-face tangent is silently applied to the interpolation/reconstruction
residual. `eta=0` produces a zero viscous matrix contribution.

| Direct localised audit | No viscous tangent | Viscous only | Viscous + passive + boundary |
|---|---:|---:|---:|
| `||J_DD v||` | 1769.3964 | 1769.3964 | 1769.3964 |
| `||P_DD v||` | 229.5759 | 1941.5092 | 2137.0676 |
| relative action error | 0.873544 | 0.170972 | 0.250962 |
| cosine | 0.977878 | 0.990991 | 0.991802 |

The isolated physical viscous action is `1712.11198693`, cosine
`0.999753654407` with the complete localised action. Thus the new term
removes the dominant missing direction and magnitude, but its compact
owner/neighbour approximation still differs from the direct `fvc::fGrad`
face reconstruction. The remaining deficiency is displacement-block
stencil/geometric accuracy, not pressure coupling.

The direct loaded two-patch run with the switch enabled still fails at the
first non-zero step (`t=0.001`, `deltaT=0.001`) with outer LGMRES
`DIVERGED_ITS` at 1000 iterations. Exact LU of the assembled `P` used
`A=MFFD`, `P=MPIAIJ`, and MUMPS; its first three linear solves converged in
one iteration with no null pivots, but the repeated nonlinear solve aborted
in the MUMPS backend with an illegal-instruction/memory-corruption signal.
This is not evidence that the physical residual is invalid, and it does not
constitute loaded acceptance. The switch remains default-off pending a
better direct-face stencil/geometric `P_DD` approximation and a robust
fieldsplit solve.

## Prestressed passive nominal-force stage

The passive contribution is bounded behind the default-off switch
`preconditionerPassiveNominalTangent false`. The Aróstica law exposes
`passiveNominalTangentField(List<tensor>&)` and evaluates the full nine-
component direct-face derivative of `Ppassive = F&S` through the shared
constitutive kernel. The solid model inserts it only for
`faceStressTreatment directConstitutive`, on the existing internal
owner-neighbour stencil, and rejects simultaneous use of the older scalar
`preconditionerMaterialTangent` switch. The physical residual is unchanged.

The independent prestressed elastic runtime-smoke reference uses a controlled
non-identity face deformation and non-zero passive nominal force. Its
directional finite-difference error is `7.30483e-11`. The copied PETSc audit
with the new switch enabled remained deterministic, with zero accepted-history
change and Stage 2 term-sum reconstruction error `1.75081968491e-16`.
At the tiny audit state, the direct localised action was
`||J_DD v||=1769.3965`, `||P_DD v||=1989.8468`, relative error `0.1854693`,
cosine `0.9916076`. This is action-level validation only; the larger affine
prestress attempts on the ventricular scratch boundary field produced invalid
`J`, so geometric stiffness and loaded convergence remain unresolved.

The elastic smoke bus error was classified as a stale executable/library ABI
mismatch after the mechanicalLaw virtual interface changed. Rebuilding the
executable made the smoke pass, and its Alltest script now rebuilds before
invocation.

The passive-nominal loaded gate was also rerun on a copied two-patch case with
`preconditionerMaterialTangent false`, `preconditionerBoundaryTangent false`,
`preconditionerViscousTangent true`, and
`preconditionerPassiveNominalTangent true`. With `-snes_mf_operator`, PETSc
reported an MFFD operator followed by an MPIAIJ preconditioner. The exact
`P` factorisation used MUMPS, converged the first three one-iteration linear
solves, and then aborted with PETSc signal 4 (illegal instruction/memory
corruption) during the next factorisation; it is not a loaded nonlinear
acceptance. The production fieldsplit configuration failed at `t=0.001` with
`SNES DIVERGED_LINEAR_SOLVE` before a converged loaded step. Therefore this
bounded contribution is validated at action level but does not close the
coupled loaded-solver gate.

The elastic, viscoelastic and boundary runtime smoke scripts now rebuild
their executable before running, and use `sed -i.bak` for macOS/Linux
portability. This prevents stale `mechanicalLaw` virtual-table layouts and
macOS `sed -i` syntax from being misreported as source failures.

## Coupled-block and Schur investigation

The retained Case B mesh1 synthetic tiny-load failure was reproduced from a
scratch copy only. The physical configuration remains unchanged: `deltaT =
0.001`, pressure amplitude `0.016071182 Pa`, direct-constitutive face stress,
the two optional direct-face tangents enabled, `scaleMixedPetScFields true`,
and physical/stabilisation scales `scaleFactor = 10` and
`scaleFactorJacobian = 10`.

The PETSc block ordering is interlaced per cell as `[D_x,D_y,D_z,pHat]`, with
physical `p = 550.000000036*pHat` and pressure equation row scale
`550.000000036`. The assembled matrix is therefore

```text
    [ P_DD  P_Dp ]
P = [            ]
    [ P_pD  P_pp ]
```

On the 17,625-cell serial scratch mesh, the current assembled block norms
were:

| block | nonzeros | Frobenius norm |
|---|---:|---:|
| `P_DD` | 737055 | `3.88677198695e5` |
| `P_Dp` | 242269 | `7.25793879554e-1` |
| `P_pD` | 242269 | `7.18883915619e-1` |
| `P_pp` | 81895 | `4.58485262431e-3` |

The corrected `P_Dp` pressure-only momentum action remains consistent with
the MFFD action. The `P_pD` action is zero for smooth and rigid-translation
modes, and is approximate for localised/high-frequency modes: the localised
relative error is `0.9232` with cosine `0.4438`; volumetric error is `0.1190`
with cosine `0.9939`. This is consistent with the compact neighbour stencil
used by `InsertFvmDivUIntoPETScMatrix` versus the wider least-squares-gradient
operator used by the residual.

The `P_pp` audit is more decisive. The constant-pressure mode matches the
finite-difference derivative to `1.1e-12`, confirming the sign and scaling of
the `-p/kappa` term. Nonconstant modes do not:

| mode | relative `J_pp`/`P_pp` action error | cosine |
|---|---:|---:|
| smooth gradient | `0.9240` | `0.4085` |
| localised pressure | `1.1205` | `0.9228` |
| checkerboard | `0.9303` | `0.3948` |
| one-cell basis | `0.9476` | `0.9341` |

An audit-only decomposition confirmed the source. The physical Rhie-Chow
residual is

```text
div(magSf*scaleFactor*(snGrad(p) - n & interpolate(grad(p))))
```

with `grad(p)` supplied by the configured `leastSquaresS4f` scheme. The
current `diffStencilLaplacianStab::scalarJacobian()` instead returns a
compact `scaleFactorJacobian*fvm::laplacian(...)` matrix. The exact derivative
of the residual contains the derivative of the interpolated least-squares
gradient and therefore has a wider stencil than the current compact P. The
observed nonconstant `P_pp` mismatch is structural, not evidence that a
different positive scalar multiplier should be used.

For isolation only, scratch copies temporarily used `scaleFactor = 0` and
`scaleFactorJacobian = 0` to remove the physical stabilisation and verify the
compressibility block. Those copies are not production configurations and
were not used for solver acceptance. The retained and production settings
remain `scaleFactor = 10` and `scaleFactorJacobian = 10`; no physical or
preconditioner block is intentionally multiplied by zero.

### Schur and scaling experiments

With `A = MFFD` and `P = MPIAIJ`, lower/upper/full/diagonal Schur variants,
`a11`, `self`, `selfp`, and additive fieldsplit variants all reached their
diagnostic 100-iteration cap. No PC zero pivot or NaN was reported. An exact
serial LU factorisation of the complete compact P was made through `PC ASM +
sub-PC LU`; factorisation succeeded with no zero pivot and PETSc displayed
`A = mffd`, `P = mpiaij`. At 100 outer iterations its true residual ratio was
approximately `1.25e-1`, so exact inversion of the compact P does not make it
an adequate complete Jacobian.

Exact block LU with a full/selfp fieldsplit reached approximately `3.09e-1`
at 100 iterations. FGMRES with the production fieldsplit was better than the
short LGMRES comparison (`5.29e-1` versus approximately `1.04` at 100
iterations), but the extended run had not converged by the audit window and
is not an acceptance result. Disabling PETSc mixed-field scaling was worse
(`1.49` at 100 iterations), so the current scaling is retained.

The coupled investigation therefore identifies a compact-stencil pressure
stabilisation/Jacobian mismatch and a broader P-versus-MFFD displacement
deficiency, but does not justify a production change to the physical
residual, to either pressure scale, or to `scaleFactor`. A wider-stencil
stabilisation Jacobian or a deliberately redesigned mixed preconditioner is
still required before the loaded case can be accepted.

### Least-squares pressure-block correction audit

The opt-in `preconditionerLeastSquaresPressureCoupling` path was rebuilt and
tested in a scratch audit. It uses the same least-squares displacement-gradient
stencil and the current `J*F^-T` determinant derivative as the pressure
residual. The localised `P_pD` action changed from a relative error of about
`0.9232` and cosine `0.4438` to relative error below `1e-9` and cosine `1`.
The physical residual and accepted histories were unchanged.

The separate opt-in
`preconditionerLeastSquaresPressureStabilisation` experiment retained the
nonzero compact `fvm::laplacian` block and added an internal-face derivative of
the interpolated least-squares pressure gradient. It did not improve the
nonconstant `P_pp` audit sufficiently and is not enabled by default. It is
diagnostic only; it must not be confused with setting `scaleFactorJacobian 0`.

All loaded scratch configurations retained `scaleFactor = 10` and
`scaleFactorJacobian = 10`. A zero Jacobian scale would remove that
preconditioner contribution and is not a valid production fix.

### Opt-in full wide production `P_pp`

The canonical switch remains
`preconditionerLeastSquaresPressureStabilisation`. It is default-disabled.
When enabled, the assembled pressure block retains the finite-bulk-modulus
and compact direct-`snGrad` terms and adds the derivative of the production
Rhie-Chow/least-squares chain:

```text
p cell -> least-squares grad(p) -> interpolated face grad(p)
        -> corrected face stabilisation flux -> cell divergence
```

The implementation is in
`nonLinGeomTotalLagTotalDispSolid.C::formJacobian`. It includes internal
faces, fixed-value boundary contributions, coupled processor-face gradient
stencils, global pressure-column indices, and `ADD_VALUES` assembly. The
physical residual, pressure values, time discretisation, and accepted
histories are not changed. `pressureUnknownScale` and `pressureEqnScale` are
applied once, while the preconditioner uses the configured
`scaleFactorJacobian` (10); the physical residual continues to use
`scaleFactor` (10).

The serial mesh1 scratch audit produced 206487 stored nonzeros and
Frobenius norm `0.00909876488258`. Its constant, smooth, localized,
checkerboard, one-cell, deterministic-random, and boundary-adjacent pressure
modes matched `J_pp^FD` with cosine 1 and relative errors from approximately
`2.45e-12` (constant mode) down to machine precision. A two-rank decomposed
audit exercised processor-face exchange and completed successfully; its local
P_pp storage was 105037 entries per rank and the same pressure-mode action
errors were obtained. Repeated assembly was deterministic and history
checksums remained unchanged.

The wide block therefore passes the pressure-only serial/parallel action
acceptance. It does not by itself solve the coupled loaded case. With the
frozen wide-P_pp runtime and the 0.016071182 Pa synthetic load, full+selfp,
lower+selfp, and full+a11 all failed at `t=0.001` by KSP iteration limit;
the true-residual ratios were approximately 1.0000966, 1.0000000, and
1.0025815 respectively. The compact Schur/factorisation choice remains
insufficient for coupled convergence.

The requested constitutive, boundary, viscoelastic, ddt-cost, mixed-solid,
and Python regression suites passed. The compiled Gultekin bulk-modulus
runtime smoke reached law construction but aborted with signal 6 during its
comparison sequence; this remains a separate unresolved runtime regression
outside the pressure-block source path.

## Global explicit J-versus-P diagnosis

The bounded follow-up investigation is documented in
`notes/arosticaGlobalJVersusPDiagnosis.md`. It uses a standalone executable
and a six-tetrahedron mixed scratch case to construct complete finite-difference
`J_FD`, compare it with PETSc MFFD, replace assembled blocks with explicit FD
blocks, and compare exact and approximate Schur operators. Production source
behavior and the accepted Case B parent remain unchanged.
