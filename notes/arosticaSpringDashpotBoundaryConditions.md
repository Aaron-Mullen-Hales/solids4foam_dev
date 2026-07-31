# Aróstica spring-dashpot boundary conditions: source equations

This note records the boundary-condition convention used by Phase 2A. The
accepted case is deliberately not changed. Its `alpha`/`beta` entries are
accepted as compatibility aliases by the new source classes.

## Benchmark convention

The Aróstica benchmark applies support tractions per undeformed area. With
`N` the reference patch normal, `D` the current trial displacement, and
`Ddot` the trial velocity, the two support laws are

\[
 t_{0,epi}=-\alpha_{epi}(D\cdot N)N
             -\beta_{epi}(Ddot\cdot N)N,
\]

and

\[
 t_{0,base}=-\alpha_{base}D-\beta_{base}Ddot.
\]

The epicardial law therefore has exactly zero tangential spring and dashpot
traction. The basal law acts independently on all three displacement
components. Both laws have zero pressure traction. The force inserted into
the total-Lagrangian residual is `t0*magSfReference`, not `t0*magSfCurrent`.

The reference implementation and accepted setup use `alphaEpi=1e8 Pa/m`,
`betaEpi=5e3 Pa s/m`, and the monoventricle basal values
`alphaBase=1e5 Pa/m`, `betaBase=5e3 Pa s/m`. The later biventricular setup
uses the same dashpot value and a basal spring coefficient of `1e6 Pa/m`.

## Source cross-check

The benchmark PDF describes the epicardial support as a normal spring and
dashpot and the basal support as a vector spring and dashpot. The pinned
Finsberg/reference implementation supplies these as nominal/reference
tractions. The accepted case's `0/D` uses the names `alpha` and `beta` and
does not explicitly write `useUndeformedArea`; this is compatible with the
new conditions because they force the reference-area route internally.

No contradiction was found in the signs or traction measure. The signs are
restoring: positive outward normal displacement produces traction in the
negative reference-normal direction. Replacing `N` by `-N` changes both
normal projections and the final normal vector, so the epicardial traction is
unchanged.

## Interface implication

`nonLinGeomTotalLagTotalDispSolid::enforceTractionBoundaries` already
recognises `solidTractionFvPatchVectorField` polymorphically. It multiplies
the returned traction by the reference patch area when
`useUndeformedArea()` is true. Phase 2A makes this accessor virtual and the
new common Aróstica base overrides it to return true. No Aróstica-specific
branch is added to the solid model.

The boundary condition obtains `D` from the object registry and evaluates
`fvc::ddt(D)`. Consequently it uses the active OpenFOAM `ddt(D)` scheme and
the current trial field together with accepted old-time fields. It does not
advance time, call `updateTotalFields()`, cache trial values, or modify old
time data.

