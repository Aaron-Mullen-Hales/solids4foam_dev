# GultekinTwoFibreElastic

`GultekinTwoFibreElastic` is an aligned two-family exponential anisotropic law
created for the Gultekin-Dal-Holzapfel artery benchmark. It is not a dispersed
Holzapfel-Gasser-Ogden model.

## Dictionary

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

The four field names are dictionary-selected. `f0` and `f1` must be
dimensionless `volVectorField` reference directions. `f0f` and `f1f` must be
dimensionless `surfaceVectorField` reference directions. Cell fields are used
only by the cell constitutive path and face fields only by the face path. The
law never reconstructs or updates these fields from deformed geometry.

All supplied directions must be finite and unit length within
`fibreUnitTolerance`, which defaults to `1e-6`. Zero directions and silent
normalisation are not allowed. The two families need not be orthogonal. When
`useSecondFibreFamily false`, the family-2 dictionary names and fields may be
omitted; the second stress is explicitly zero.

## Energy and stress

The matrix energy is

```text
Wiso = mu/2 (I1bar - 3),
I1bar = J^(-2/3) tr(C).
```

For each enabled family,

```text
Wa = k1/(2 k2) [exp(k2 Ea^2) - 1].
```

For the default unsplit benchmark form,

```text
Ea  = I4a - 1,
I4a = Ma.C.Ma.
```

For `anisotropicSplit true`,

```text
Ea      = I4bar_a - 1,
I4bar_a = J^(-2/3) I4a.
```

The split stress includes the full derivative of `J^(-2/3) I4a`; it is not
formed by blindly applying `dev()` to a previously calculated unsplit stress.
The law returns Cauchy stress and explicitly divides the Kirchhoff terms by
`J`. The passive contributions are

```text
sigmaIso = mu/J dev(J^(-2/3) b),

sigmaFibre_a = 2 k1 Ea exp(k2 Ea^2)/J (ma ma)          [unsplit],

sigmaFibre_a = 2 k1 Ea exp(k2 Ea^2)/J
               J^(-2/3) dev(ma ma)                    [split].
```

No `dev()` is applied to the complete unsplit anisotropic stress.

When `fibresTensionOnly true`, `Ea` in the energy and stress coefficient is
replaced by its positive part. The default is `false`, preserving the
benchmark's compressed-fibre response.

## Mixed pressure

Both `correct()` overloads return passive Cauchy stress only. The law does not
look up `p` or `pf`. With `nonLinGeomTotalLagTotalDispSolid` and
`solvePressure true`, the solid model owns the pressure operation and forms

```text
sigmaTotal = dev(sigmaPassive) - p I.
```

Its momentum faces interpolate this completed cell stress; the solver does not
call the face-law overload or use a law-level `pf` path. Pressure is therefore
inserted exactly once.

In the exactly incompressible `K -> infinity` limit, this is an admissible
pressure gauge for displacement and total stress. Relative to a representation
`sigmaPassive - pFull I`, the solver reports

```text
p = pFull - tr(sigmaPassive)/3.
```

Thus, in that limit, solver `p` must be shifted by
`tr(sigmaPassive)/3` before comparison with the full-unsplit-energy
multiplier. `GTF_sigmaPassive` retains the complete unprojected constitutive
tensor.

`bulkModulus` is required, finite, and dictionary-controlled. The mixed solver
caches `rKappa = 1/bulkModulus` and includes the physical term

```text
-p/bulkModulus
```

in the pressure residual. A finite value gives a nonzero diagonal response to
a spatially constant pressure and therefore defines that mode in the
slightly-compressible mixed formulation. Exact incompressibility is not
silently selected. At finite `K`, `p` is not freely gauge-shiftable: replacing
it by `p - tr(sigmaPassive)/3` would generally no longer satisfy the same
`-p/K` pressure equation. The solver's `dev(sigmaPassive) - p I` split is the
specified finite-`K` continuation of the benchmark.

The law does not insert `bulkModulus*(J - 1) I`, `bulkModulus*log(J) I`, or
any other volumetric penalty stress. At finite constraint error, the
stabilised mixed equations are the solver's numerical continuation rather
than the full off-constraint unsplit energy derivative; monitor `GTF_J` when
assessing benchmark agreement.

At `J=1`, split and unsplit energies and deviatoric stresses agree, while their
full stresses differ only by a spherical tensor that mixed pressure absorbs.
They are not generally equivalent when `J` differs from one.

## Implicit stiffness

The original HGO finite-perturbation `muEff` machinery is not copied. The
recommended interface supplies the unscaled numerical parameter explicitly:

```text
implicitShearModulus implicitShearModulus [1 -1 -2 0 0 0 0] 2.01e6;
```

If it is omitted, the fallback estimate is

```text
implicitShearModulus = mu + 2 Nfamilies k1.
```

This is a conservative reference-state scalar based on the largest
single-family small-strain fibre modulus. It is not the exact exponentially
growing large-strain tangent. `impKcoeff` scales the specified or fallback
value. Construction reports the source, unscaled value, coefficient, and
effective numerical stiffness separately.

Both `impK()` and `shearModulus()` return this numerical value. This matters
because the current target mixed solver uses `2*shearModulus()` instead of the
law's `impK()` for momentum stabilisation and its approximate displacement
Jacobian. Increasing the value can improve robustness but can also add excess
numerical stabilisation.

`implicitShearModulus` is separate from `bulkModulus`: changing the former
changes the numerical displacement operator and stabilisation stiffness,
whereas changing the latter changes `1/K` in the physical pressure equation.
`impKcoeff` scales only `implicitShearModulus` and never `bulkModulus`.

## Robustness and diagnostics

The law terminates if required `bulkModulus` is missing, has non-pressure
dimensions, or is non-finite or non-positive. It also terminates on other
invalid material dimensions or values, non-positive or non-finite `J`,
invalid invariants, invalid fibre fields, or non-finite stress.
By default exponential arguments above 650 are clipped and a warning reports
the count, maximum unclipped argument, and limit. Set `clipExponent false` to
terminate instead, or choose a lower positive `exponentLimit` up to 650.

When `writeDiagnostics true`, these cell fields are registered and written:

- `GTF_J`;
- `GTF_I4` and `GTF_I6`;
- `GTF_lambda4` and `GTF_lambda6`;
- `GTF_sigmaIso`;
- `GTF_sigmaFibre1` and `GTF_sigmaFibre2`;
- `GTF_sigmaPassive`.

`I4`, `I6`, and their square-root stretches are the invariants actually used
by the selected split or unsplit formulation. `GTF_sigmaPassive` excludes
mixed pressure.

## Material test

Run the independent analytical regression suite with:

```bash
python3 GultekinTwoFibreElastic/tests/test_gultekin_two_fibre.py
```

It checks undeformed and rotated states, objectivity, extensions, shear,
family interchange and disabling, split/unsplit behaviour, tension-only
response, finite-difference energy-stress consistency, and homogeneous
cell/face equivalence. It also confirms that numerical stiffness and finite
bulk modulus remain separate from constitutive stress.

The compiled one-cell test uses a registered minimal `solidModel`, calls both
`correct()` overloads for homogeneous simple shear, checks against an
independent analytical stress, registers nonzero `p` and `pf` to prove they
are ignored, doubles `implicitShearModulus`, doubles `bulkModulus`, checks all
three modulus APIs, exercises a disabled second family, verifies invalid bulk
moduli fail, and instantiates the historical HGO law:

```bash
source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc
sh GultekinTwoFibreElastic/tests/runtimeSmoke/Alltest
```
