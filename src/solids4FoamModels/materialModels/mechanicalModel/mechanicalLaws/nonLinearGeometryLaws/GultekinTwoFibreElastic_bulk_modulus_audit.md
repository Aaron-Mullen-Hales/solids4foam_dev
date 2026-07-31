# GultekinTwoFibreElastic bulk-modulus source-path audit

## Scope and conclusion

This audit was completed against branch `solid-pressure-land3Fix` at commit
`ebf998ff25945f611a44e5376656690c62768d32` before editing the material law.
It covers `GultekinTwoFibreElastic`, `mechanicalLaw`, `mechanicalModel`,
`nonLinGeomTotalLagTotalDispSolid`, its mixed residual and approximate
Jacobian, the PETSc preconditioning path, pressure stabilisation, pressure
boundary handling, all target-solver uses of `mechanical().bulkModulus()`, all
target-solver uses of `rKappa()`, and `updateF()`.

The checked-out source is consistent with the requested architecture:

```text
GultekinTwoFibreElastic::correct()
    -> passive isochoric matrix plus fibre Cauchy stress only

nonLinGeomTotalLagTotalDispSolid
    -> dev(sigmaPassive) - p I in the momentum stress
    -> 1/K in the mixed pressure equation
```

No source-path conflict was found. A finite dictionary-controlled `K` can be
returned through `bulkModulus()` without inserting a volumetric penalty
stress in either constitutive `correct()` overload.

## 1. Exact pressure residual term containing `1/K`

`nonLinGeomTotalLagTotalDispSolid::makeRKappa()` constructs and caches

```text
rKappa = 1/mechanical().bulkModulus().
```

For `solvePressure true`, `formResidual()` evaluates the cell pressure
residual as

```text
R_p = pressureEqnScale * V
      *[-p*rKappa
        + pressureStabilisation.cellScalar(rAUf)
        - 0.5*(J^2 - 1)/J].
```

Consequently, for a homogeneous material,

```text
-p*rKappa = -p/K.
```

The same unscaled decomposition is used by
`updatePressureConstraintDiagnostics()`.

## 2. Constant pressure mode

The pressure contribution to the momentum equation is a gradient or surface
traction of `-p I`. The pressure-stabilisation operators in this path are
Laplacian-like and do not, by themselves, give a nonzero response to a
spatially constant pressure when no pressure value boundary fixes that mode.

A uniform perturbation `deltaP` changes the physical pressure residual by

```text
deltaR_p = -pressureEqnScale*V*deltaP/K.
```

For every finite positive `K`, this provides a nonzero cell-local pressure
derivative and defines the constant pressure mode. Taking `K` to an effectively
infinite value removes that anchoring term and recovers the incompressible
null-mode limit unless another constraint or boundary condition fixes it.

The finite term also defines the slightly-compressible mixed relation. Ignoring
stabilisation for illustration, the pressure row gives

```text
p/K = -0.5*(J^2 - 1)/J.
```

This is a mixed pressure constraint; it is not a volumetric stress added by
the material law.

## 3. Residual, Jacobian, preconditioner, stabilisation, and boundaries

### Physical residual

`K` enters the physical pressure residual directly and only through `-p/K`.
It does not enter the displacement residual as a constitutive penalty. The
solid model applies the pressure stress separately through
`dev(sigmaPassive) - p I`.

### Residual and unknown scaling

`K` does not set `pressureEqnScale_`. That scale is

```text
pressureScaleFactor
* (pressureScaleByTwoMu ? volumeAverage(2*shearModulus) : 1).
```

Nor does `K` set `pressureUnknownScale_`; that scale is one, a user scalar, or
the same representative `2*shearModulus`. Scaling multiplies the complete
pressure row, including `-p/K`, but is calculated independently of `K`.

### Pressure Jacobian

For the scaled PETSc unknown `pHat`, where
`p = pressureUnknownScale*pHat`, `formJacobian()` assembles

```text
-pressureEqnScale*pressureUnknownScale*fvm::Sp(rKappa, p)
+ pressureEqnScale*pressureUnknownScale
  *pressureStabilisation.scalarJacobian(p, rAUf).
```

Thus finite `1/K` appears explicitly on the pressure-pressure block diagonal.
The displacement-pressure gradient and pressure-displacement divergence blocks
do not use `K`.

### Preconditioner and Schur approximation

`foamPetscSnesHelper` passes the matrix populated by the target model's
`formJacobian()` to PETSc as the assembled Jacobian/preconditioner matrix.
`nonLinGeomTotalLagTotalDispSolid` has no separate custom pressure
preconditioner and no in-source explicit Schur formula. If PETSc field-split
or Schur options are selected externally, they operate on the assembled
blocks, so their pressure block contains the same finite `1/K` diagonal.

### Pressure stabilisation

`K` is not an input to pressure stabilisation. The stabilisation residual and
Jacobian use `p` and `rAUf`. `rAUf` is built from the approximate momentum
diagonal using `impKf`, density, mesh, and time-step information. This keeps
`bulkModulus` distinct from `implicitShearModulus`.

### Boundary-condition pathways

The pressure field is read if present and otherwise created with
`zeroGradient` boundaries. Its boundary conditions are corrected before the
stress split and residual evaluation. `K` is not used to construct or modify
pressure boundary conditions, traction boundary conditions, or fibre fields.
The `-p/K` contribution is cell-local; there is no flux or gradient of `K` in
this pressure path.

## 4. Evaluation frequency

`rKappa()` is lazy. Its first call invokes `makeRKappa()`, which calls
`mechanical().bulkModulus()` and stores the resulting reciprocal field in
`rKappaPtr_`. Later residual, diagnostic, and Jacobian calls reuse the cached
field and do not re-evaluate the law's `bulkModulus()`.

For one material, `mechanicalModel::bulkModulus()` forwards directly to the
single law. For multiple materials, it evaluates each submesh law and maps the
result once when `makeRKappa()` is first called. There is no runtime dictionary
reread or refresh of `rKappaPtr_` in the target solver.

## 5. `updateF()` volumetric argument

`updateF()` always updates the deformation gradient. Its `mu` and `K`
arguments are otherwise used only if the solid model's `enforceLinear` switch
forces the nonlinear law into a linear-elastic fallback stress.

`GultekinTwoFibreElastic::correct()` currently passes a zero volumetric
modulus to both cell and face `updateF()` calls. This should remain unchanged.
It ensures that even an enforced-linear fallback returns only deviatoric
passive stress; the mixed solid model still owns `-p I` and the pressure row
still obtains the physical finite `K` independently through `bulkModulus()`.
Passing the new finite mixed `K` into `updateF()` would create a law-level
volumetric stress in the fallback path and would violate the requested stress
ownership.

## 6. Search for other volumetric constitutive stress

The target law's matrix stress is

```text
mu/J*dev(J^(-2/3)*b),
```

and its fibre terms are unchanged passive Cauchy stresses. Neither
`correct()` overload reads `p` or `pf`. Both assign the shared helper's
`sigmaPassive` result.

After the law update, `nonLinGeomTotalLagTotalDispSolid` applies

```text
sigmaTotal = dev(sigmaPassive) - p I.
```

The only volumetric stress in the target momentum path is therefore the
solver-owned mixed pressure stress. The finite bulk modulus belongs to the
pressure constraint through `1/K`; no `K*(J - 1) I`, `K*log(J) I`, or other
volumetric penalty stress exists in the checked target path.

## Implementation consequence

The safe localized change is to:

- require a finite dimensioned `bulkModulus` entry in the Gultekin law;
- validate and report it;
- return its value from `bulkModulus()`;
- keep both `correct()` overloads and their zero-`K` `updateF()` calls
  unchanged;
- keep pressure and volumetric stress entirely in the mixed solid model.
