# GultekinTwoFibreElastic stress-interface audit

## Scope and baseline

This audit was completed before creating or editing any constitutive-law source
file. It covers the current working-tree versions of:

- `mechanicalLaw.H` and `mechanicalLaw.C`;
- `mechanicalModel.C`;
- `nonLinGeomTotalLagTotalDispSolid.H` and
  `nonLinGeomTotalLagTotalDispSolid.C`;
- `HolzapfelGasserOgdenElastic.H` and
  `HolzapfelGasserOgdenElastic.C`;
- `GuccioneElastic.H` and `GuccioneElastic.C`;
- `electroMechanicalLaw.H` and `electroMechanicalLaw.C`;
- `neoHookeanElastic.H` and `neoHookeanElastic.C`;
- `Make/files.openfoam`, `Make/files.foamextend`, and `Make/options`.

The audit was made on branch `solid-pressure-land3Fix` at starting commit
`ebf998ff25945f611a44e5376656690c62768d32`. The worktree already contained
uncommitted changes to the mixed solver, Guccione and electro-mechanical laws,
the two build lists, and unrelated files. Those pre-existing changes are user
work and are not part of this implementation.

The starting SHA-256 hashes of the source templates are:

```text
420b55c0121c64e86befe08d97470b6b7ddd4b466faa1adf27016c15453d659b  HolzapfelGasserOgdenElastic.H
807c453b0e969aefad4cdbb1b71266cede6f9d63f02ea7d6a0c5c9918c6c301b  HolzapfelGasserOgdenElastic.C
```

## Stress measure returned by `mechanicalLaw::correct`

`mechanicalLaw::correct` does not impose one stress measure at the abstract
base-class level; the comment on `sigma0` explicitly says that nonlinear laws
may choose a law-specific measure. The actual consumer fixes the contract. For
`nonLinGeomTotalLagTotalDispSolid`, and for the nonlinear laws audited here,
the returned `volSymmTensorField` is Cauchy stress.

The evidence is:

1. `neoHookeanElastic` explicitly evaluates and documents Cauchy stress.
2. `GuccioneElastic` pushes second-Piola stress forward as
   `symm(F & S & F.T())/J` before assigning `sigma`.
3. `electroMechanicalLaw` pushes active second-Piola stress forward with the
   same `1/J` conversion and adds the resulting active Cauchy stress.
4. `HolzapfelGasserOgdenElastic` labels its result as Cauchy stress.
5. The solid model contracts `sigma` with the current normal and current area,
   which is the spatial Cauchy-traction operation.

The new law must therefore return Cauchy stress, not Kirchhoff stress and not
second-Piola stress.

## Entry into the total-Lagrangian momentum residual

The solver computes

```text
SfCurrent = interpolate(J F^-T) & SfReference
nCurrent  = SfCurrent/|SfCurrent|
traction  = nCurrent & interpolate(sigma)
force     = |SfCurrent| traction
residual  = div_reference(force) + body/inertia/damping/stabilisation terms
```

This is the discrete Nanson transformation. In continuum notation,

```text
da n = J F^-T dA N,
P = J sigma F^-T,
Div_X(P) = 0.
```

There is no additional constitutive stress conversion after
`mechanical().correct(sigma())`. The solid model performs the equivalent
Cauchy-to-first-Piola force transformation through the deformed face area.

## Mixed pressure insertion and sign

For `solvePressure true`, the current `nonLinGeomTotalLagTotalDispSolid` calls
`applyMixedPressureStressSplit()` after the material update. For an ordinary
passive law it applies

```text
sigma = dev(sigma) - p I.
```

The electro-mechanical special case preserves the complete active Cauchy
stress while deviatorically projecting only the passive part.

Consequently:

- positive `p` is compressive pressure;
- hydrostatic Cauchy stress is `-p I`;
- the solid model inserts the final mixed-pressure contribution once;
- the passive material law should return passive Cauchy stress only.

Several mixed-capable laws also return `-p I` when the pressure field is
available, so that their direct `correct()` result is meaningful to other
consumers. That pattern is redundant for this new target-specific law because
the solver owns pressure explicitly. `GultekinTwoFibreElastic` therefore does
not inspect either `p` or `pf`.

### Cell and face pressure ownership trace

The target cell path is:

```text
mechanical().correct(sigma)       -> sigmaPassive
applyMixedPressureStressSplit()   -> dev(sigmaPassive) - p I
fvc::interpolate(sigma)           -> current-face total Cauchy stress
```

The PETSc residual path performs both operations in `unpackSolution()`. The
predictor also applies the split after predicting `p`. The law must not insert
cell pressure before either call.

The pressure row uses

```text
-p rKappa + pressureStabilisation - (J^2 - 1)/(2J) = 0,
rKappa = 1/bulkModulus.
```

This law requires a finite dictionary-controlled bulk modulus. The resulting
`-p/K` term defines the slightly-compressible pressure constraint and controls
the constant pressure mode. It remains a mixed pressure equation parameter,
not a material volumetric penalty stress. At a finite discrete constraint
residual, pressure stabilisation and the chosen off-constraint equations
select a numerical continuation of the gauge; they do not restore the
spherical unsplit passive term inside the law.

The target momentum path does not call
`mechanical().correct(surfaceSymmTensorField&)`. It interpolates the completed
cell stress instead. Consequently `pf` handling inside the law would be stale
for this solver. The face overload still returns the same passive
constitutive stress for consumers that explicitly request face stress; any
solver-level face pressure remains the caller's responsibility.

## Effect of `dev()` on anisotropic stress

The solid model does apply `dev()` to the passive law output when
`solvePressure true`. It does not do so when mixed pressure is disabled.

A non-deviatoric passive anisotropic stress can cross the `mechanicalLaw`
interface, but its spherical component cannot survive the current mixed-solver
post-processing as a separately identifiable constitutive contribution. The
new law must nevertheless calculate and return the complete anisotropic
Cauchy tensor; it must not pre-emptively apply `dev()` to the unsplit fibre
stress.

For an exactly incompressible material this projection does not discard a
measurable constitutive degree of freedom. If

```text
sigma = sigmaPassive - p I,
```

then the equivalent representation

```text
sigma = dev(sigmaPassive) - pMixed I,
pMixed = p - tr(sigmaPassive)/3
```

produces the same total stress. The mixed pressure variable is therefore a
shifted Lagrange multiplier. In particular,

```text
pMixed = pFull - tr(sigmaPassive)/3,
pFull  = pMixed + tr(sigmaPassive)/3.
```

For the unsplit benchmark this is an admissible and intended incompressible
pressure gauge for displacement and total-stress predictions: the constraint
only fixes stress up to an isotropic Lagrange multiplier, and the two
representations are exactly equivalent at `J=1`. It is not the same reported
pressure as the multiplier associated with the chosen off-constraint
extension of the full unsplit energy. Any comparison of the solver `p` field
with that alternative multiplier must apply the shift above. This equivalence
does not make the original anisotropic stress deviatoric; diagnostics retain
the unprojected passive contribution. At finite `J - 1`, the current mixed
equations are the solver's discrete continuation and should not be described
as the full off-constraint unsplit energy derivative; `J` must therefore be
monitored in the benchmark.

## Required `1/J` conversion

For energy per unit reference volume, the push-forward chain is

```text
S       = 2 dW/dC,
tau     = F S F^T,
sigma   = tau/J.
```

The unsplit fibre expression supplied for this task is Kirchhoff stress:

```text
tauAniso = sum 2 k1 (I4a - 1)
           exp(k2 (I4a - 1)^2) symm(ma ma).
```

The Cauchy stress returned through the audited interface must therefore be
`tauAniso/J`.

The existing HGO implementation omits `1/J` in both matrix and fibre terms and
contains the explicit comment `J=1`. That is intentional for its stated exact
incompressibility assumption: when `J=1`, Kirchhoff and Cauchy stresses are
identical, and isotropic shifts can be absorbed into pressure. It should not be
copied into a numerically near-incompressible mixed implementation, where
`J` is solved only to tolerance. The new law must retain the explicit `1/J`.

## Cell and face equivalence in the existing HGO law

The existing HGO cell and face formulas are algebraically identical. They
differ only in field location and names:

- cells use `F`, `p`, `Ec`, `Ea`, and `Er`;
- faces use `Ff`, `pf`, `Ecf`, `Eaf`, and `Erf`.

Both paths construct the same two reference directions, invariants, matrix
stress, pressure term, and two exponential fibre stresses. Both omit `1/J`.
There is no face-only push-forward or pressure-sign difference.

The new implementation should use one scalar/tensor constitutive helper for
both locations so split selection, tension-only selection, `J`, invariant
evaluation, exponential protection, and stress conversion cannot drift.

## `impK`, effective shear modulus, and the target mixed solver

The existing HGO law:

1. estimates three shear responses by finite perturbation;
2. computes a cellwise `muEff` from perturbed stresses in a principal-stress
   coordinate system;
3. returns `impKcoeff*muEff` from `impK()`;
4. refreshes `muEff` only from `updateTotalFields()`.

`calcInitialShearModulus()` only reports its estimate because the assignment to
`mu_` is commented out. The finite-difference machinery is heuristic, uses a
fixed perturbation, and does not provide a consistent anisotropic material
tangent. It can also become very large under the exponential law.

More importantly, the current `nonLinGeomTotalLagTotalDispSolid` does not call
the law's `impK()` when `solvePressure true`. At construction it uses

```text
impK = 2*mechanical().shearModulus().
```

That field enters the momentum-stabilisation face coefficient, the segregated
approximate displacement Jacobian/preconditioner, the reciprocal diagonal
used by pressure stabilisation, and pressure-row scaling. The original HGO law
does not override `shearModulus()`, so its inherited base implementation is not
suitable for this target path.

For `mu = 10 kPa`, `k1 = 500 kPa`, and `k2 = 2`, copying the HGO machinery is
not justified. The fibre tangent is already much larger than the ground-matrix
modulus near the reference state and grows exponentially in extension. The new
law should expose a documented, controlled constant reference stiffness for
the mixed solver and should not claim that it is the exact large-strain
tangent. This affects conditioning and stabilisation, not the constitutive
stress.

## Model represented by the existing HGO class

Despite its name, the audited `HolzapfelGasserOgdenElastic` is an aligned
two-family exponential law. It constructs two directions from `Ec`, `Ea`,
`Er`, and `fibreAngle`, and uses one `I4` invariant for each family. It has no
dispersion parameter, no generalised structure tensor, and no distributed
orientation integration. It is not a complete dispersed
Holzapfel-Gasser-Ogden implementation.

## Implementation consequences

The independent `GultekinTwoFibreElastic` law will therefore:

- return Cauchy stress and include explicit `1/J` conversion;
- read reference cell and face fibre fields by dictionary-selected names;
- retain the full unsplit anisotropic tensor inside the law and diagnostics;
- derive the split fibre stress directly from `I4bar`, including its
  volumetric derivative, rather than applying `dev()` to an unsplit result;
- return passive Cauchy stress only from both `correct()` overloads and leave
  all `p` and `pf` ownership to the calling solver;
- provide both `shearModulus()` and `impK()` with a documented numerical
  stiffness strategy suitable for the current target solver;
- leave both original HGO source files unchanged.
