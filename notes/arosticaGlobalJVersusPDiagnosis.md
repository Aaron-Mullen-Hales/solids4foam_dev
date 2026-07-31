# Aróstica global explicit J-versus-P diagnosis

## Scope

This investigation is standalone and diagnostic only. It does not modify the
production solid model, `foamPetscSnesHelper`, physical residual, pressure
field, accepted histories, or the accepted Case B parent case.

The executable is
`tests/ArosticaViscoelastic/globalJVersusPDiagnosis/ArosticaGlobalJVersusPDiagnosis`.
It uses the existing public `initialiseSolution`, `initialiseJacobian`,
`formResidual`, and `formJacobian` APIs. PETSc MFFD is constructed locally with
`formResidual` as its callback; no production accessor is required.

## Reduced case

The fixture is a six-tetrahedron unit cube with base, endocardium, and
epicardium patches. It retains mixed displacement/pressure unknowns, finite
bulk modulus, anisotropic fibre/sheet/normal fields, viscoelasticity, inertia,
follower pressure, and finite basal/epicardial spring-dashpot support.

The mixed ordering is cell-major:

```text
(D_x, D_y, D_z, pHat)
```

There are six cells and therefore 24 unknowns. The expected blocks are
18×18, 18×6, 6×18, and 6×6.

## Explicit J construction

For every owned column, the diagnostic evaluates the complete production
residual at positive and negative perturbations and inserts

```text
J_FD e_i = (R(x + epsilon e_i) - R(x - epsilon e_i))/(2 epsilon)
```

The default epsilon pair is `1e-5` and `1e-6`; their relative matrix change is
reported. Histories and old-time state are not advanced by the diagnostic.
The assembled matrix is obtained from `formJacobian` and is named `P` in the
reports.

Serial runs write Matrix Market files and metadata under the temporary case's
`globalJVersusPExport` directory. Parallel runs exercise PETSc global column
ownership and report ownership ranges; dense spectral analysis is serial-only.

## Diagnostics

The executable supports action, block, Schur, spectral, and displacement
inverse studies through the switches documented by `-help`. Block studies
replace one selected block of `P` with the corresponding explicit FD block
and run a GMRES/LU solve against the explicit `J_FD`. Schur studies compute

```text
S_J = J_pp - J_pD J_DD^-1 J_Dp
S_P = P_pp - P_pD P_DD^-1 P_Dp
```

The test runner executes both the action-exact wide pressure block and the
compact pressure block in separate temporary cases. No case output is written
to the accepted Aróstica parent.

The reduced dense spectral path uses Eigen bundled with the repository. It
reports eigenvalues, singular values, condition estimates, and departure from
normality for `J`, `P`, and `P^-1 J`. This avoids any dependency on petsc4py,
NumPy, or SciPy.

## Runtime guard

Explicit export is guarded above 64 cells. The exact Case B mesh1 has 70,500
mixed unknowns; full explicit export would require at least one residual pair
per column and is not part of the normal diagnosis. The executable can instead
run action probes on a supplied scratch copy with `-probeFullCaseModes`.

## Regression boundary

The wide/compact reduced diagnostics are separate from the existing action,
history, restart, zero-load, boundary, and material-law tests. The compiled
Gultekin smoke abort remains a separate shared-library regression and is not
interpreted as Aróstica solver evidence.
