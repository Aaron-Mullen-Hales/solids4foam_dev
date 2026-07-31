# Gültekin mixed finite-volume solid implementation report

Date: 2026-07-18

Source workspace:

    /Volumes/OpenFoam/aaronmullen-hales-v2312/s4f-development

This implementation is an isolated, runtime-selectable solids4foam
finite-volume mixed pressure-displacement model. It targets the continuum
equations used by the original Gültekin–Dal–Holzapfel Q1P0+WAS benchmark. It
does not implement Q1 or P0 finite elements, Hu–Washizu element assembly,
finite-element shape functions, element stiffness matrices, or FE quadrature.

No artery case was modified or run. No Git command was used.

## Implementation architecture

The existing `nonLinGeomTotalLagTotalDispSolid` stores its stress split,
deformation Jacobian, pressure residual, pressure diagnostics, and several
support methods as private members. A direct subclass could therefore not
override the continuum contract safely. Copying its approximately 6,000-line
implementation would create a fragile solver fork.

The implementation uses a small derived runtime class and four protected
virtual continuum hooks in the existing model:

1. retain the complete passive stress in the mixed split;
2. evaluate the positive constraint function `g(J)`;
3. evaluate its scalar derivative `dg/dJ`;
4. provide a diagnostic expression for `g(J)`.

The generic implementations exactly preserve the established behavior:

    retainFullPassiveStressInMixedSplit() = false

    gGeneric(J) = 0.5*(J^2 - 1)/J
                = 0.5*(J - 1/J)

    dgGeneric/dJ = 0.5*(1 + 1/J^2)

    sigmaGeneric = dev(sigmaPassive) - p*I

The isolated derived model overrides those hooks with:

    retainFullPassiveStressInMixedSplit() = true

    gGultekin(J) = (J - 1)/J
                 = 1 - 1/J

    dgGultekin/dJ = 1/J^2

    sigmaGultekin = sigmaPassive - p*I

No existing case selects the new runtime type, so no existing case activates
the new branch. The generic runtime name and its dictionary remain unchanged.

## Runtime model

New runtime model name:

    nonLinGeomTotalLagTotalDispGultekinSolid

The compiled library contains the runtime-table constructor symbol:

    addnonLinGeomTotalLagTotalDispGultekinSolid
        dictionaryConstructorTosolidModelTable_

The model requires:

- `solvePressure true`;
- exactly one `GultekinTwoFibreElastic` law;
- `anisotropicSplit false`;
- `fibresTensionOnly false`;
- `useSecondFibreFamily true`.

The constructor terminates with a clear error if any of those paper-contract
requirements is not satisfied.

At startup it prints once:

    Gultekin mixed continuum mode active

    total Cauchy stress:
        sigma = sigmaPassive - p I

    pressure constraint:
        -p/K -(J - 1)/J + S_p = 0

    volumetric Cauchy stress when S_p = 0:
        -p I = K*(J - 1)/J I

## Exact continuum equations

`GultekinTwoFibreElastic` remains unchanged. In its unsplit branch it returns
passive spatial Cauchy stress:

    sigmaPassive = sigmaIso + sigmaFibre1 + sigmaFibre2

    sigmaIso = (mu/J)*dev(bBar)

    sigmaFibre1 =
        2*k1*(I4 - 1)*exp(k2*(I4 - 1)^2)/J * (m1 tensor m1)

    sigmaFibre2 =
        2*k1*(I6 - 1)*exp(k2*(I6 - 1)^2)/J * (m2 tensor m2)

The new solid model does not apply `dev()` to this sum. It inserts pressure
exactly once:

    sigma = sigmaIso
          + sigmaFibre1
          + sigmaFibre2
          - p*I

The pressure equation is:

    R_p = -p/K -(J - 1)/J + S_p

With `S_p=0` and `R_p=0`:

    p = -K*(J - 1)/J

    -p*I = K*(J - 1)/J*I

This is the volumetric Cauchy stress obtained from:

    U(J) = K*(J - ln(J) - 1)

The pressure sign is the existing solids4foam convention. It was not
reversed.

## Residual and derivative paths

The production pressure residual now calls the runtime constraint hook. For
the new model it assembles:

    pressureEqnScale*V*
    (-p*rKappa + pressureStabilisation.cellScalar(rAUf) -(J - 1)/J)

The following paths are consistent:

- PETSc matrix-free residual: evaluates the complete new nonlinear `g(J)`;
- pressure diagnostics: use `-(J - 1)/J` through the same virtual hook;
- Land pressure decomposition diagnostics: print and consume the selected
  runtime constraint rather than a hard-coded old expression;
- residual row scaling and pressure-unknown scaling: unchanged and applied to
  the completed residual as before;
- pressure diagonal and Schur approximation: the `-rKappa` and pressure
  stabilisation blocks are unchanged;
- D-in-pressure approximate Jacobian: explicitly multiplies its established
  `-V*div(D)` stencil by `dg/dJ` evaluated at `J=1`.

Both constraints have derivative one at the reference state:

    dgGeneric/dJ at J=1 = 1
    dgGultekin/dJ at J=1 = 1

Consequently the existing approximate-Jacobian numbers are unchanged, while
the code now documents why. The exact matrix-free directional derivative uses
`dgGultekin/dJ=1/J^2` through the nonlinear residual.

The derivative sign in the residual is negative:

    d[-gGultekin(J)]/dJ = -1/J^2

## Cell, boundary, and face stress paths

Execution order remains:

1. unpack PETSc displacement;
2. update `gradD`, `F`, `Finv`, and `J`;
3. call `mechanical().correct(sigma)`;
4. unpack the physical pressure and correct its boundary conditions;
5. apply the runtime mixed stress split.

For the new model, step 5 is exactly:

    sigma() = sigma() - p()*I;
    return;

There is no preceding or subsequent `dev()` on this execution path.

The material cell correction and every material boundary value contain the
complete passive stress. Subtracting `p*I` operates on the complete volume
field, including its boundary fields. The PETSc momentum residual then forms
the current physical traction from:

    nCurrent & fvc::interpolate(sigma())

Thus the actual face force contains the complete unsplit fibre spherical
part. The material-law face `correct(surfaceSymmTensorField&)` contract also
remains passive-only, but this PETSc momentum path uses interpolated completed
cell total stress. The law never inserts `p` or `pf`.

Pressure is therefore counted exactly once.

## Boundary-condition compatibility

No boundary-condition class was changed. The new model inherits the existing
total-Lagrangian traction handling:

- fixed and prescribed displacement/rotation remain kinematic conditions;
- `solidTraction` and follower pressure overwrite the physical face force
  with the requested current- or reference-area traction according to the
  existing `useUndeformedArea` setting;
- traction-free boundaries continue to receive zero prescribed traction;
- symmetry/slip zero-shear projection remains unchanged;
- reaction and postprocessing access the final written `sigma` field.

The pressure-load sign convention is unchanged.

## Pressure and momentum stabilisation

The pressure stabilisation remains structurally:

    +S_p

No stabilisation type, coefficient, or discretisation was changed. The source
and optional fields distinguish:

    constitutive constraint = -p/K -(J - 1)/J
    numerical pressure term = S_p
    complete residual       = constitutive constraint + S_p

The momentum stabilisation interface and `diffStencilLaplacian` implementation
are unchanged. No artery-specific scale or numerical choice is embedded in
the new model. Continuum equivalence alone does not imply equivalence while a
large finite-mesh numerical stabilisation term is active.

## Optional diagnostics

The model coefficient:

    writeGultekinMixedDiagnostics false;

is disabled by default. When enabled at a write time it writes:

| field | class | dimensions | definition |
|---|---|---|---|
| `sigmaPassiveFull` | volSymmTensorField | Pa | `sigma + p I` |
| `fibreMeanStress` | volScalarField | Pa | `tr(sigmaPassiveFull)/3` |
| `sigmaVolumetricFromP` | volSymmTensorField | Pa | `-p I` |
| `gultekinConstraintPhysical` | volScalarField | dimensionless | `-p/K -(J-1)/J` |
| `pressureStabilisationContribution` | volScalarField | dimensionless | `S_p` |
| `completePressureResidual` | volScalarField | dimensionless | physical plus stabilisation |

These diagnostics do not change the residual and are not allocated or written
unless explicitly requested.

## Material-law audit

`GultekinTwoFibreElastic.C` and `.H` were not modified.

Source and compiled smoke checks confirm:

- `sigmaIso` is trace-free;
- unsplit `sigmaFibre1` is a complete dyad;
- unsplit `sigmaFibre2` is a complete dyad;
- `sigmaPassive` is their direct sum;
- cell and face `correct()` overloads return passive Cauchy stress;
- nonzero registered `p` and `pf` are ignored by the law;
- no volumetric stress or fibre-mean removal occurs in the law.

`HolzapfelGasserOgdenElastic` was not modified.

## Automated tests

### New formulation tests

Command:

    PYTHONDONTWRITEBYTECODE=1 \
        python3 tests/GultekinMixedSolid/test_gultekin_mixed_solid.py

Result: 10 tests passed.

The suite covers:

- K=5e6 Pa scalar pressure at J=1.0, 1.1, 1.2, 1.4, and 1.6;
- zero-stabilisation pressure-residual closure to machine precision;
- analytic `1/J^2` and centred finite-difference derivatives at all five J;
- retention of a synthetic nonzero-trace fibre tensor;
- generic-model `dev(passive)-pI` regression contract;
- generic pressure constraint unchanged;
- pressure inserted exactly once in the new path;
- final written `sigma` remains spatial Cauchy stress in Pa;
- unchanged material cell and face passive-stress contract;
- runtime registration and absence of FE assembly code.

For the requested synthetic tensors:

    sigmaIso   = diag(-1, 0.5, 0.5)
    sigmaFibre = diag( 0, 6.0, 3.0)
    p          = -2

the test obtains:

    sigmaNew   = diag( 1, 8.5, 5.5)
    wrong dev  = diag(-2, 5.5, 2.5)

### Compiled material-law smoke tests

Command:

    source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc
    cd src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/\
        nonLinearGeometryLaws/GultekinTwoFibreElastic/tests/runtimeSmoke
    sh Alltest

Result: all compiled tests passed on a temporary synthetic mesh.

Notable checks:

- seven cell values and six face values;
- maximum passive-stress error approximately 9.06e-10 Pa;
- both `correct()` overloads passed;
- nonzero p and pf were ignored;
- bulk-modulus and implicit-stiffness separation passed;
- missing, dimensionally wrong, zero, negative, NaN, and Inf K inputs failed as
  expected.

The older optional NumPy material-point script was not run because NumPy is
not installed. No dependency was installed; the compiled OpenFOAM smoke suite
provides the relevant cell/face contract test.

## Build

Exact command:

    source /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc
    wmake libso src/solids4FoamModels

Result: success, exit code 0.

Environment:

    OpenFOAM: v2312
    WM_OPTIONS: darwin64ClangDPInt32Opt

Compiled library:

    /Users/aaronmullen-hales/OpenFOAM/aaronmullen-hales-v2312/
        platforms/darwin64ClangDPInt32Opt/lib/
        libsolids4FoamModels.dylib

Final SHA-256:

    7dc01ef6627628edf03f0a9c91fc4853e98706ba99b3b1a3d602f693e5064f65

Clang emitted no warning for either changed translation unit. `wmkdepend`
reported unresolved legacy compatibility-header names including
`PointPatchFieldMapper.H`, `globalPointPatchFields.H`, `cyclicGgiPolyPatch.H`,
`cyclicGgiFvPatchFields.H`, `BasicSymmetryPointPatchField.H`,
`ComponentMixedPointPatchVectorField.H`, `PrimitivePatchTemplate.H`, and
`PrimitivePatchInterpolationTemplate.H`. These dependency-scanner notices
also arise from the inherited existing include graph; compilation and linking
resolved the active OpenFOAM-v2312 paths successfully.

## Files

Added:

    src/solids4FoamModels/solidModels/
        nonLinGeomTotalLagTotalDispGultekinSolid/
            nonLinGeomTotalLagTotalDispGultekinSolid.H
            nonLinGeomTotalLagTotalDispGultekinSolid.C
            README.md

    tests/GultekinMixedSolid/
        test_gultekin_mixed_solid.py
        README.md

    GultekinMixedSolid_implementation_report.md

Modified:

    src/solids4FoamModels/solidModels/
        nonLinGeomTotalLagTotalDispSolid/
            nonLinGeomTotalLagTotalDispSolid.H
            nonLinGeomTotalLagTotalDispSolid.C

    src/solids4FoamModels/Make/files
    src/solids4FoamModels/Make/files.foamextend

Generated by the normal build:

- `lnInclude` links;
- dependency and object files beneath
  `src/solids4FoamModels/Make/darwin64ClangDPInt32Opt`;
- the compiled `libsolids4FoamModels.dylib` above.

No material law, boundary condition, stabilisation model, PETSc helper, or
artery case file was modified.

## Minimal case handoff

The case chat should change only the runtime solid-model selection and the
coefficient-subdictionary name while retaining its controlled numerical
settings. A minimal `solidProperties` structure is:

    solidModel nonLinGeomTotalLagTotalDispGultekinSolid;

    nonLinGeomTotalLagTotalDispGultekinSolidCoeffs
    {
        solutionAlgorithm PETScSNES;
        solvePressure true;

        // Copy the existing pressure/momentum stabilisation, scaling,
        // time-control and solver entries without silently changing them.

        writeGultekinMixedDiagnostics false;
    }

The mechanical law remains:

    type                  GultekinTwoFibreElastic;
    bulkModulus           K [1 -1 -2 0 0 0 0] 5e6;
    anisotropicSplit      false;
    fibresTensionOnly     false;
    useSecondFibreFamily  true;

The exact full material parameter block should continue to come from the
controlled benchmark configuration. The case should load the compiled library
whose SHA-256 is recorded above and require the startup line:

    Gultekin mixed continuum mode active

No case was changed automatically in this source session.

## Known limitations

- No nonlinear solid simulation or artery test was run by instruction.
- Runtime registration was verified from the compiled constructor-table
  symbol and runtime strings, not by launching a solid case.
- The generic approximate Jacobian remains a reference-state approximation;
  this is intentional and unchanged. PETSc matrix-free differentiation sees
  the exact new residual.
- Pressure and momentum stabilisation remain in the finite-volume equations.
  Their benchmark influence must be assessed separately.
- Matching the original continuum equations does not ensure identical
  coarse-mesh Q1P0 FE and cell-centred FV solutions.
- The new runtime class deliberately rejects multi-material, split-fibre,
  tension-only, and one-family configurations because they are not the target
  Figure 6(b) continuum model.

## Final source audit

The new execution path satisfies from source:

    sigma = sigmaIso
          + sigmaFibre1
          + sigmaFibre2
          - p*I

and:

    R_p = -p/K -(J - 1)/J + S_p

Confirmed:

- generic runtime behavior remains `dev(sigmaPassive)-pI` with its original
  `0.5*(J^2-1)/J` constraint;
- no fibre mean stress is removed in the new path;
- pressure is inserted exactly once;
- internal, boundary, interpolated face force, written stress, traction, and
  postprocessing paths use the completed total Cauchy stress;
- no FE solver or FE assembly code was added;
- no artery-specific stabilisation or load setting was embedded;
- no simulation was run;
- no Git command was run.

