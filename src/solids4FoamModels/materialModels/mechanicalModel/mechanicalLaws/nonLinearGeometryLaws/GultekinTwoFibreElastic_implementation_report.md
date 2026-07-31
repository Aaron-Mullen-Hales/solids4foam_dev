# GultekinTwoFibreElastic implementation report

## Result

An independent, runtime-selectable `GultekinTwoFibreElastic` mechanical law
has been added for the Gultekin-Dal-Holzapfel two-fibre artery benchmark. It
uses dictionary-selected reference cell and face fibre fields and does not
construct artery directions from geometry. The affected OpenFOAM library
builds, all 14 independent analytical checks pass, and a compiled one-cell
test calls both `correct()` overloads and compares their stresses with an
independent analytical result. A required finite dictionary-controlled bulk
modulus now supplies `1/K` to the mixed pressure equation without changing
constitutive stress. The original `HolzapfelGasserOgdenElastic` law also
instantiates through its historical dictionary interface.

The artery tube case was not run, and no commit was created.

## Baseline and scope

- Branch: `solid-pressure-land3Fix`
- Starting commit: `ebf998ff25945f611a44e5376656690c62768d32`
- OpenFOAM: `v2312`, build patch `240625`
- `WM_OPTIONS`: `darwin64ClangDPInt32Opt`
- `WM_COMPILER`: `Clang`
- Compiler: `/usr/bin/clang++`, Apple clang version 17.0.0

The worktree was already dirty at the start. It contained user changes to the
mixed solver, Guccione and electro-mechanical laws, both build lists, and
unrelated tutorial and utility files. Those changes were preserved. Apart
from the newly created files, this implementation adds only the
`GultekinTwoFibreElastic.C` entry to each existing solids4foam build list.

## Source-interface audit

The full read-only audit was written before constitutive source editing in
`GultekinTwoFibreElastic_stress_interface_audit.md`. Its findings are:

1. The target nonlinear law interface returns Cauchy stress.
2. The total-Lagrangian solver evaluates current traction with the current
   area vector. By Nanson's formula this is equivalent to using
   `P = J sigma F^-T` in the reference momentum residual.
3. With `solvePressure true`, the current solver inserts pressure after the
   law update through `sigma = dev(sigma) - p I` for a passive law.
4. Positive `p` is compressive, so hydrostatic Cauchy stress is `-p I`.
5. The target mixed solver does apply `dev()` to the passive law output.
6. The law interface itself carries a complete non-deviatoric anisotropic
   tensor. In the mixed solver its spherical part is absorbed into the
   pressure Lagrange multiplier rather than remaining separately identifiable.
7. The existing HGO law's missing explicit `1/J` is intentional under its
   documented `J=1` assumption. The new near-incompressible mixed law must
   explicitly convert Kirchhoff to Cauchy stress with `1/J`.
8. The existing HGO cell and face formulas are algebraically identical.
9. The HGO finite-perturbation `muEff` is a heuristic, not a consistent
   tangent. More importantly, the current mixed solver uses
   `2*mechanical().shearModulus()` rather than `impK()`.
10. The existing HGO class is an aligned two-family exponential law, not a
    complete dispersed Holzapfel-Gasser-Ogden implementation.

A follow-up ownership trace confirmed that the target mixed solver expects
passive Cauchy stress from the law. Its cell path calls `correct(sigma())` and
then applies `sigma = dev(sigma) - p I`. Its momentum path obtains face
traction by interpolating that completed cell stress; it does not call the
surface `correct()` overload. The surface overload remains available for
other callers and likewise returns passive stress. Consequently neither
overload reads `p` or `pf`.

## Bulk-modulus source-path audit

The follow-up audit in `GultekinTwoFibreElastic_bulk_modulus_audit.md` was
also written before editing the law. It confirmed:

1. `makeRKappa()` lazily caches
   `1/mechanical().bulkModulus()` once.
2. The physical pressure residual contains
   `-p/K + pressureStabilisation - (J^2 - 1)/(2J)` before volume and row
   scaling.
3. Finite `1/K` gives a nonzero response to a constant pressure perturbation
   and therefore defines that mode.
4. The same `1/K` appears in the pressure-pressure approximate Jacobian and
   any PETSc preconditioner or external Schur split using that assembled
   block.
5. `K` does not set pressure-row scaling, pressure-unknown scaling,
   `rAUf`, pressure stabilisation, or pressure/traction boundary conditions.
6. The zero volumetric argument to both `updateF()` calls must remain zero so
   an enforced-linear fallback also remains passive and mixed-pressure owned.
7. No target path adds a `K`-based volumetric constitutive stress.

## Why a new law was created

`HolzapfelGasserOgdenElastic.H` and `.C` were used only as structural source
templates. Their runtime name, artery-coordinate field interface,
constitutive equations, and implicit-stiffness machinery were not changed or
shared. A separate class avoids changing behaviour for existing cases and
gives the benchmark an explicit direct-reference-fibre interface.

The new class is declared and registered as:

```cpp
class GultekinTwoFibreElastic
TypeName("GultekinTwoFibreElastic")
```

It has its own `addToRunTimeSelectionTable` entry and derives directly from
`mechanicalLaw`.

## Mathematical formulation

The passive strain energy per unit reference volume is

```text
W = mu/2 (I1bar - 3)
  + sum_a k1/(2 k2) [exp(k2 Ea^2) - 1],

I1bar = J^(-2/3) tr(C),
C = F.T() & F,
ma = F & Ma.
```

For the default unsplit benchmark formulation,

```text
Ea = I4a - 1,
I4a = Ma.C.Ma.
```

For `anisotropicSplit true`,

```text
Ea = I4bar_a - 1,
I4bar_a = J^(-2/3) I4a.
```

The implemented Cauchy stresses are

```text
sigmaIso = mu/J dev(J^(-2/3) b),

sigmaFibre_a = 2 k1 Ea exp(k2 Ea^2)/J (ma ma)
               [unsplit],

sigmaFibre_a = 2 k1 Ea exp(k2 Ea^2)/J
               J^(-2/3) dev(ma ma)
               [split].
```

The split tensor follows by differentiating `I4bar` and is not produced by
blindly projecting an unsplit stress. No `dev()` is applied to the complete
unsplit anisotropic stress.

When `fibresTensionOnly true`, the energy strain and stress coefficient use
the positive part `max(Ea, 0)`. Its default is `false`, so compressed fibres
retain the benchmark's squared-strain response.

### Stress transformation chain and pressure

For energy per reference volume,

```text
S = 2 dW/dC,
tau = F & S & F.T(),
sigma = tau/J.
```

The new law performs the final `1/J` conversion explicitly and returns only
`sigmaPassive` from both overloads. The target solid model owns pressure and
forms `dev(sigmaPassive) - pMixed I`, so pressure is present exactly once.

In the exactly incompressible `K -> infinity` limit, this projection can be
represented as a pressure-gauge choice rather than a change to the unsplit
material response. The full-stress representation

```text
sigmaTotal = sigmaPassive - pFull I
```

is identical to the mixed-solver representation when

```text
pMixed = pFull - tr(sigmaPassive)/3,
pFull  = pMixed + tr(sigmaPassive)/3.
```

For an exactly incompressible formulation the Lagrange multiplier is defined
only up to this constitutive spherical shift, so the solver split produces the
same displacement and total stress. A reported `pMixed` must be gauge-shifted
before comparison with a source that reports `pFull`. This equivalence does
not make the unsplit fibre stress deviatoric.

For the required finite `K`, however, `-p/K` fixes the pressure constant mode.
The shifted `pFull` above would generally not satisfy the same finite-`K`
pressure equation. Therefore `dev(sigmaPassive) - pMixed I` is the specified
slightly-compressible mixed model, not a freely shiftable pressure gauge away
from the incompressible limit.

At exactly `J=1`, split and unsplit energies and deviatoric stresses agree.
Their full stresses differ by a spherical tensor because the consistent split
derivative contains a volumetric term; this difference is absorbed by mixed
pressure. At `J != 1`, the two options are physically distinct.

## Fibre-field interface and validation

The complete benchmark dictionary interface is:

```text
type                 GultekinTwoFibreElastic;

rho                  rho [1 -3 0 0 0 0 0] 1000;
mu                   mu  [1 -1 -2 0 0 0 0] 10000;
k1                   k1  [1 -1 -2 0 0 0 0] 500000;
k2                   k2  [0 0 0 0 0 0 0] 2;

bulkModulus          K
                     [1 -1 -2 0 0 0 0] 1.0e6;

fibreField1          f0;
fibreField2          f1;
faceFibreField1      f0f;
faceFibreField2      f1f;

anisotropicSplit     false;
fibresTensionOnly    false;
useSecondFibreFamily true;

implicitShearModulus implicitShearModulus
                     [1 -1 -2 0 0 0 0] 2.01e6;
impKcoeff            1.0;

fibreUnitTolerance   1e-6;
clipExponent         false;
exponentLimit        650;
writeDiagnostics     true;
```

The cell fields are dimensionless `volVectorField`s and the face fields are
dimensionless `surfaceVectorField`s. All are immutable reference directions.
The constructor checks existence, dimensions, finite values, nonzero
magnitude, minimum and maximum magnitude, and unit length within
`fibreUnitTolerance` (default `1e-6`). It fails on non-unit data rather than
normalising it silently. It does not assume the families are orthogonal.

When `useSecondFibreFamily false`, family-2 names and fields may be omitted;
family-2 stress is explicitly zero. A zero vector is never accepted as a
family-disable signal.

Cell updates use the selected cell fields, face updates use the selected face
fields, and explicit face fibres are never replaced by interpolation. Both
paths call the same pointwise constitutive helper.

## Robustness

The required `bulkModulus` has no fallback. A missing entry produces a clear
`FatalIOError`; wrong dimensions, zero, negative, NaN, and Inf values all
fail. No arbitrary upper limit is imposed. The implementation also terminates
on other invalid material dimensions or values, negative `mu` or `k1`,
non-positive `k2`, invalid fibre data, non-finite or non-positive `J`,
non-finite invariants, and non-finite stress.

Exponential overflow protection is controlled by `clipExponent` and
`exponentLimit`. The default clips arguments above 650 and emits a warning
containing the clipped cell/face count, maximum unclipped argument, and limit.
With clipping disabled, exceeding the configured limit is fatal. Ordinary
states are not silently clipped.

## Finite mixed bulk modulus

The law stores the required dimensioned dictionary value as
`bulkModulus_`. Its `bulkModulus()` API returns a uniform field containing
that finite value. Construction reports `K`, `1/K`, and `K/mu`; for the
recommended values it reports `K/mu = 100`.

The target solver consumes this API only through its lazily cached
`rKappa = 1/K` pressure field. Both constitutive overloads continue to pass
zero volumetric modulus to `updateF()` and return only the isochoric matrix
plus fibre Cauchy stress. No volumetric penalty stress was added.

## Implicit stiffness strategy

The existing HGO perturbation-based `muEff` machinery was not copied.
`implicitShearModulus` is an explicitly reported numerical parameter. The
recommended benchmark dictionary supplies it dimensionally; if omitted, the
documented fallback is

```text
mu + 2 Nfamilies k1.
```

`impKcoeff` scales the explicit value or fallback. Construction reports the
source, unscaled modulus, coefficient, and effective value. Both `impK()` and
`shearModulus()` return the effective value. For `mu = 10 kPa`, `k1 = 500
kPa`, and two families the fallback is 2.01 MPa. This scalar is a numerical
stiffness for stabilisation and the segregated approximate Jacobian, not the
exact large-strain exponential tangent. It may require case-specific tuning
under severe extension and torsion.

`bulkModulus` and `implicitShearModulus` are independent. `impKcoeff` scales
only the latter. The compiled test changes each parameter separately and
verifies that the other API and passive constitutive stress do not change.

## Optional diagnostics

When `writeDiagnostics true`, the law allocates and writes cell fields
`GTF_J`, `GTF_I4`, `GTF_I6`, `GTF_lambda4`, `GTF_lambda6`, `GTF_sigmaIso`,
`GTF_sigmaFibre1`, `GTF_sigmaFibre2`, and `GTF_sigmaPassive`. They are not
created when diagnostics are disabled. The fields are filled by the same
constitutive helper used for stress; `GTF_sigmaPassive` is the complete
unprojected passive Cauchy tensor and excludes mixed pressure.

## Files and build-system changes

New implementation and documentation files:

- `GultekinTwoFibreElastic/GultekinTwoFibreElastic.H`
- `GultekinTwoFibreElastic/GultekinTwoFibreElastic.C`
- `GultekinTwoFibreElastic/README.md`
- `GultekinTwoFibreElastic_stress_interface_audit.md`
- `GultekinTwoFibreElastic_bulk_modulus_audit.md`
- `GultekinTwoFibreElastic_implementation_report.md`
- `GultekinTwoFibreElastic_build.log`

New independent analytical test:

- `GultekinTwoFibreElastic/tests/test_gultekin_two_fibre.py`

New one-cell runtime-smoke utility and fixture:

- `GultekinTwoFibreElastic/tests/runtimeSmoke/GultekinTwoFibreLawSmoke.C`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/Alltest`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/Make/files`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/Make/options`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/D`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/pointD`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Ea`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Eaf`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Ec`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Ecf`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Er`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/Erf`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/f0`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/f0f`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/f1`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/0/f1f`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/dynamicMeshDict`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/g`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtf`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfInfBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfMissingBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfNaNBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfNegativeBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfOneFamily`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfStiff`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfWrongBulkDimensions`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.gtfZeroBulk`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/lawProperties.hgo`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/physicsProperties`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/constant/solidProperties`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/system/blockMeshDict`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/system/controlDict`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/system/fvSchemes`
- `GultekinTwoFibreElastic/tests/runtimeSmoke/case/system/fvSolution`

One source entry was added to each existing build list:

- `src/solids4FoamModels/Make/files.openfoam`
- `src/solids4FoamModels/Make/files.foamextend`

`Make/options` did not require a change. The OpenFOAM-v2312 build was tested;
the foam-extend build-list entry was inspected but not compiled in this
environment.

For the finite-bulk-modulus follow-up, the implementation/header, README,
Python test, compiled smoke source, smoke driver, GTF test dictionaries,
bulk-modulus audit, implementation report, and combined build/test log were
updated. No solid-model, artery-case, HGO, or build-list source was changed by
that follow-up.

## Material-level verification

The independent Python implementation does not call the C++ constitutive
helper. It evaluates the energy and stresses analytically and uses explicit
absolute, relative, and finite-difference tolerances.

| Check | Definition | Result |
| --- | --- | --- |
| A | `F = I`: `J=1`, `I4=I6=1`, zero passive stress | Pass |
| B | Rigid rotation and rotated-state objectivity | Pass |
| C | Isochoric uniaxial extension | Pass |
| D | Extension along family 1 | Pass |
| E | Extension transverse to both families | Pass |
| F | Simple shear | Pass |
| G | Fibre-family interchange invariance | Pass |
| H | Explicitly disabled family 2 gives exactly zero family-2 stress | Pass |
| I | Split/unsplit pressure-equivalence at `J=1`; difference at `J!=1` | Pass |
| J | Tension-only changes compressed response only | Pass |
| K | Finite-difference energy-stress consistency, split and unsplit | Pass |
| L | Homogeneous cell/face formula consistency | Pass |
| M | Implicit-stiffness sensitivity; constitutive stress unchanged | Pass |
| N | Finite `-p/K`, invalid-K rejection, and constitutive separation | Pass |

Automated result:

```text
PASS: 14 material-level checks; atol=2.0e-10, rtol=2.0e-09, energy-stress rtol=3.0e-06
```

Exact material-test command, from `nonLinearGeometryLaws`:

```bash
python3 GultekinTwoFibreElastic/tests/test_gultekin_two_fibre.py
```

## Build and runtime-selection verification

The exact library build command, from `src/solids4FoamModels`, was:

```bash
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc && wmake libso
```

Result: pass.

```text
Library: /Users/aaronmullen-hales/OpenFOAM/aaronmullen-hales-v2312/platforms/darwin64ClangDPInt32Opt/lib/libsolids4FoamModels.dylib
SHA-256: 0b09a5ab890bcdd5255da7079997325a24b5f84d75cb249e86c4f984ffdd6de2
```

The exact runtime-smoke utility build command, from
`GultekinTwoFibreElastic/tests/runtimeSmoke`, was:

```bash
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc && wmake
```

The complete compiled material-test command, from `nonLinearGeometryLaws`,
was:

```bash
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc &&
sh GultekinTwoFibreElastic/tests/runtimeSmoke/Alltest
```

It called both `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)` for homogeneous simple shear. The cell and
face results each agreed with the independent C++ analytical stress to
`9.05733e-10 Pa`. Deliberately nonzero registered `p` and `pf` fields were
ignored. Doubling `implicitShearModulus` from 2.01 MPa to 4.02 MPa doubled
both reported numerical-stiffness APIs while changing constitutive stress by
exactly zero and leaving `bulkModulus()` at 1 MPa. Doubling `bulkModulus` from
1 MPa to 2 MPa changed every returned bulk-modulus value by exactly two while
leaving `shearModulus()`, `impK()`, cell stress, and face stress unchanged.
The explicit one-family dictionary also passed both overload checks with
maximum error `9.05148e-10 Pa`.

Missing, dimensionally wrong, zero, negative, NaN, and Inf bulk-modulus
dictionaries each produced the expected fatal construction failure. The
primary `K = 1 MPa`, `mu = 10 kPa` construction log reported `1/K = 1e-6
1/Pa` and `K/mu = 100`.

The exact original-HGO smoke-test command was:

```bash
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc &&
GultekinTwoFibreLawSmoke \
    -case /tmp/GTF_bulkModulusStress_20260717 \
    -dict lawProperties.hgo
```

Result:

```text
PASS: instantiated runtime mechanical law HolzapfelGasserOgdenElastic from lawProperties.hgo
```

The original-law smoke establishes historical dictionary and runtime-table
compatibility without changing the HGO implementation.

The concise command-and-result record is in
`GultekinTwoFibreElastic_build.log`.

## Original HGO regression safety

Before and after SHA-256 values are identical:

```text
HolzapfelGasserOgdenElastic.H
before: 420b55c0121c64e86befe08d97470b6b7ddd4b466faa1adf27016c15453d659b
after:  420b55c0121c64e86befe08d97470b6b7ddd4b466faa1adf27016c15453d659b

HolzapfelGasserOgdenElastic.C
before: 807c453b0e969aefad4cdbb1b71266cede6f9d63f02ea7d6a0c5c9918c6c301b
after:  807c453b0e969aefad4cdbb1b71266cede6f9d63f02ea7d6a0c5c9918c6c301b
```

Both original files are byte-for-byte unchanged. Their git diff is empty,
and both runtime types instantiate from the rebuilt library.

## Recommended artery-benchmark settings

Use `nonLinGeomTotalLagTotalDispSolid`, `solvePressure true`, direct reference
fields `f0`, `f1`, `f0f`, and `f1f`, and use:

```text
type                  GultekinTwoFibreElastic;
rho                   rho [1 -3 0 0 0 0 0] 1000;
mu                    mu  [1 -1 -2 0 0 0 0] 10000;
k1                    k1  [1 -1 -2 0 0 0 0] 500000;
k2                    k2  [0 0 0 0 0 0 0] 2;
bulkModulus           K
                      [1 -1 -2 0 0 0 0] 1.0e6;
fibreField1           f0;
fibreField2           f1;
faceFibreField1       f0f;
faceFibreField2       f1f;
anisotropicSplit      false;
fibresTensionOnly     false;
useSecondFibreFamily  true;
implicitShearModulus  implicitShearModulus
                      [1 -1 -2 0 0 0 0] 2.01e6;
impKcoeff             1.0;
fibreUnitTolerance    1e-6;
clipExponent          false;
exponentLimit         650;
writeDiagnostics      true;
```

Enable diagnostics during initial one-step or reduced-load checks to inspect
`J`, fibre invariants, and the separated stress terms. The recommended entry
keeps exponent clipping disabled so an excessive exponent is fatal. If
clipping is enabled for diagnosis, treat every clipping warning as evidence
that the state, load increment, or numerical stiffness needs investigation.

## Unresolved concerns and limits

- In the `K -> infinity` incompressible limit, the mixed solver's projection
  of passive spherical stress can be represented by a pressure-gauge shift.
  With finite `K`, `-p/K` defines the pressure mode, so the solver's
  `dev(sigmaPassive) - p I` split is the specified slightly-compressible model,
  not a freely shiftable pressure gauge away from that limit.
- At finite `J - 1`, pressure stabilisation and the mixed pressure equation
  define the solver's off-constraint continuation; it is not the complete
  unsplit energy derivative away from incompressibility. Benchmark runs
  should monitor `GTF_J`.
- The constant implicit stiffness is deliberately conservative and
  controllable, but is not a consistent tangent for exponential large-strain
  extension or torsion. Benchmark convergence may require tuning
  `impKcoeff` or supplying `implicitShearModulus`.
- The compiled C++ test exercises both constitutive overloads on a one-cell
  mesh, but it is not a complete solid solve. The full artery case was
  intentionally not run.
- The foam-extend build entry was added, but only OpenFOAM-v2312 was compiled
  in this session.
- Fibre fields are required to be unit length; upstream field generation must
  satisfy the configured tolerance because the material performs no silent
  normalisation.

## Reproducibility commands

```bash
# Affected library, from src/solids4FoamModels
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc && wmake libso

# Analytical material tests, from nonLinearGeometryLaws
python3 GultekinTwoFibreElastic/tests/test_gultekin_two_fibre.py

# Compiled cell/face, bulk-API, invalid-input, and original-HGO tests
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc &&
sh GultekinTwoFibreElastic/tests/runtimeSmoke/Alltest

# Original-law smoke test, after preparing the supplied one-cell fixture
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc &&
GultekinTwoFibreLawSmoke \
    -case /tmp/GTF_bulkModulusStress_20260717 \
    -dict lawProperties.hgo
```

## Final worktree snapshot

The exact final `git status --short` output is:

```text
 M src/solids4FoamModels/Make/files.foamextend
 M src/solids4FoamModels/Make/files.openfoam
 M src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GuccioneElastic/GuccioneElastic.C
 M src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GuccioneElastic/GuccioneElastic.H
 M src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.C
 M src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.H
 M src/solids4FoamModels/numerics/stabilisationModels/stabilisationModel/stabilisationModel.H
 M src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/nonLinGeomTotalLagTotalDispSolid.C
 M src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/nonLinGeomTotalLagTotalDispSolid.H
 M tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/0/T
 M tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/1/T
 M tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/2/T
 M tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/3/T
 M tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/4/T
?? applications/utilities/stabilisationModelAudit/
?? src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic/
?? src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_build.log
?? src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_bulk_modulus_audit.md
?? src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_implementation_report.md
?? src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_stress_interface_audit.md
?? src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowAffineExactStab/
?? src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowHighPassStab/
?? src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowZeroForceStab/
?? tutorials/solids/hyperelasticity/idealisedVentricle/0/
?? tutorials/solids/hyperelasticity/idealisedVentricle/constant/
?? tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.hypre
?? tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.ilu
?? tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.lu
?? tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.seg.hypre
?? tutorials/solids/hyperelasticity/idealisedVentricle/system/
```

The exact final `git diff --stat` output is:

```text
 src/solids4FoamModels/Make/files.foamextend        |    4 +
 src/solids4FoamModels/Make/files.openfoam          |    4 +
 .../GuccioneElastic/GuccioneElastic.C              | 3657 ++++++++++++--
 .../GuccioneElastic/GuccioneElastic.H              |  103 +-
 .../electroMechanicalLaw/electroMechanicalLaw.C    | 1117 ++++-
 .../electroMechanicalLaw/electroMechanicalLaw.H    |   36 +
 .../stabilisationModel/stabilisationModel.H        |    8 +
 .../nonLinGeomTotalLagTotalDispSolid.C             | 5260 +++++++++++++++++++-
 .../nonLinGeomTotalLagTotalDispSolid.H             |  254 +
 .../hotCylinderTemperatureField/0/T                |    4 +-
 .../hotCylinderTemperatureField/1/T                |    4 +-
 .../hotCylinderTemperatureField/2/T                |    4 +-
 .../hotCylinderTemperatureField/3/T                |    4 +-
 .../hotCylinderTemperatureField/4/T                |    4 +-
 14 files changed, 9752 insertions(+), 711 deletions(-)
```

This stat is dominated by pre-existing user work and, by Git design, excludes
all untracked new `GultekinTwoFibreElastic` files. The only tracked-file edits
made for this implementation are one build-list line in each `Make/files`
variant.
