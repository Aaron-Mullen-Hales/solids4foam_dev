# `electroMechanicalLaw` deformation-gradient ownership audit

## Conclusion

The bug is confirmed. Before the fix, `electroMechanicalLaw` used its own
inherited `FPtr_` for the active-stress push-forward, while its nested passive
`GuccioneElastic` object updated a different `FPtr_` during every call to
`correct(volSymmTensorField&)`.

In the audited Land Problem 3 L1 case, the outer field was stale rather than
identity: it read an existing `0/F_` and never changed. If no readable `F_`
exists, the same bug leaves the outer field at the identity created by
`mechanicalLaw::makeF()`.

## Diagnostic case and method

The diagnostic used a copied Land Problem 3 L1 case only. It compared, during
ordinary volume-field residual evaluations:

- the inherited outer `electroMechanicalLaw::F()`;
- the deformation gradient owned by the passive `GuccioneElastic` law;
- `I + lookupObject<volTensorField>("grad(D)").T()`; and
- the solid-model registry field `F`.

Norms are volume-weighted RMS and cellwise Frobenius Linf values. The reported
block at each load level is the last instrumented block at the converged state.

## Deformation-gradient evidence

| Time | Field comparison | Linf | RMS |
| ---: | --- | ---: | ---: |
| 0.01 | `F_outer - I` | 0.110524 | 0.0290247 |
| 0.01 | `F_passive - I` | 0.109260 | 0.0286759 |
| 0.01 | `F_from_gradD - I` | 0.109260 | 0.0286759 |
| 0.01 | `F_outer - F_passive` | 0.00397314 | 0.00115384 |
| 0.01 | `F_outer - F_from_gradD` | 0.00397314 | 0.00115384 |
| 0.01 | `F_passive - F_from_gradD` | 0 | 0 |
| 0.5 | `F_outer - I` | 0.110524 | 0.0290247 |
| 0.5 | `F_passive - I` | 1.46705 | 0.483792 |
| 0.5 | `F_from_gradD - I` | 1.46705 | 0.483792 |
| 0.5 | `F_outer - F_passive` | 1.40590 | 0.462954 |
| 0.5 | `F_outer - F_from_gradD` | 1.40590 | 0.462954 |
| 0.5 | `F_passive - F_from_gradD` | 0 | 0 |
| 1.0 | `F_outer - I` | 0.110524 | 0.0290247 |
| 1.0 | `F_passive - I` | 2.02554 | 0.617194 |
| 1.0 | `F_from_gradD - I` | 2.02554 | 0.617194 |
| 1.0 | `F_outer - F_passive` | 1.99728 | 0.601461 |
| 1.0 | `F_outer - F_from_gradD` | 1.99728 | 0.601461 |
| 1.0 | `F_passive - F_from_gradD` | 0 | 0 |

The solid-model `F` agreed with `F_from_gradD` with zero reported Linf and RMS
error at all three load levels.

| Time | `det(F_outer)` min/max | `det(F_passive)` min/max | `det(F_from_gradD)` min/max |
| ---: | --- | --- | --- |
| 0.01 | 0.999372 / 1.00087 | 0.999401 / 1.00077 | 0.999401 / 1.00077 |
| 0.5 | 0.999372 / 1.00087 | 0.975390 / 1.03615 | 0.975390 / 1.03615 |
| 1.0 | 0.999372 / 1.00087 | 0.964426 / 1.04622 | 0.964426 / 1.04622 |

## Object identity and disk output

The addresses were stable throughout the run:

| Object | Address | Name | Registered | Registry-owned | Write option |
| --- | ---: | --- | ---: | ---: | ---: |
| Outer `F` | 932224048 | `F_` | 0 | 0 | 16 (`AUTO_WRITE`) |
| Passive `F` | 933295024 | `F_` | 1 | 0 | 16 (`AUTO_WRITE`) |
| Registry `F_` | 933295024 | `F_` | 1 | 0 | 16 (`AUTO_WRITE`) |
| Solid-model `F` | 662739504 | `F` | 1 | 0 | 16 (`AUTO_WRITE`) |

The passive law is constructed first and successfully registers `F_`. The
outer object subsequently creates another `F_`, but that object is not the
registry entry. Consequently, time-directory writes serialize the current
passive `F_` and hide the stale field used by the active law. At time 1, disk
files `F` and `F_` were byte-for-byte identical after removing only their
different `object` header lines.

## Active-stress impact before the fix

The two diagnostic stresses were

```text
sigma_outer  = F_outer  (Ta f0f0) F_outer.T  / det(F_outer)
sigma_correct = F_gradD (Ta f0f0) F_gradD.T / det(F_gradD)
```

| Time | Stress Linf | Stress RMS | Relative RMS | Residual Linf | Residual RMS | Torque outer | Torque correct | Torque difference |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 2.72830 | 0.637143 | 0.00106740 | 7.57889e3 | 4.74350e2 | -4.53710e-5 | -4.55555e-5 | -1.84506e-7 |
| 0.5 | 4.68868e4 | 1.27301e4 | 0.449183 | 1.32316e8 | 8.12689e6 | -2.14622e-3 | 1.42265e-3 | 3.56887e-3 |
| 1.0 | 1.17314e5 | 3.14482e4 | 0.534458 | 3.26899e8 | 2.07304e7 | -4.16439e-3 | 5.49303e-3 | 9.65742e-3 |

The residual difference is the divergence of the difference between the two
active surface-force fields. Torque is the corresponding global boundary
torque about the ventricular z axis.

### Axial-layer differences

Each row reports stress RMS/Linf, residual RMS/Linf, and z-torque difference.
Layers are equal-width bins over cell-centre z, from apex-side layer 0 to
base-side layer 14.

| Time | Layer | Stress RMS | Stress Linf | Residual RMS | Residual Linf | Torque difference |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 0 | 1.28358 | 2.72830 | 1.95164e3 | 7.57889e3 | 1.59798e-8 |
| 0.01 | 1 | 1.22391 | 2.12520 | 1.11968e3 | 6.58249e3 | -1.75434e-8 |
| 0.01 | 2 | 0.865974 | 1.18768 | 541.060 | 809.328 | 1.50458e-8 |
| 0.01 | 3 | 0.758700 | 1.09117 | 519.459 | 775.705 | 1.87401e-8 |
| 0.01 | 4 | 0.711863 | 1.01245 | 485.505 | 714.171 | 1.84582e-8 |
| 0.01 | 5 | 0.676422 | 0.938091 | 439.850 | 636.940 | 1.31766e-8 |
| 0.01 | 6 | 0.644638 | 0.869041 | 385.616 | 555.794 | 2.28701e-9 |
| 0.01 | 7 | 0.650834 | 0.795638 | 357.323 | 468.660 | -1.25526e-7 |
| 0.01 | 8 | 0.611451 | 0.716864 | 306.844 | 384.714 | -1.41695e-7 |
| 0.01 | 9 | 0.546990 | 0.631615 | 236.450 | 308.647 | -4.64951e-9 |
| 0.01 | 10 | 0.465238 | 0.535730 | 162.052 | 169.398 | 1.84774e-7 |
| 0.01 | 11 | 0.551290 | 0.566524 | 228.640 | 252.244 | -1.88744e-7 |
| 0.01 | 12 | 0.461441 | 0.546235 | 171.562 | 220.339 | -2.66209e-9 |
| 0.01 | 13 | 0.415653 | 0.524864 | 161.330 | 201.598 | -1.60555e-9 |
| 0.01 | 14 | 0.388662 | 0.519399 | 156.043 | 179.866 | -1.54836e-9 |
| 0.5 | 0 | 25164.1 | 46886.8 | 3.14278e7 | 1.32316e8 | -1.45573e-5 |
| 0.5 | 1 | 23479.7 | 34588.1 | 1.95705e7 | 1.11694e8 | 6.26581e-4 |
| 0.5 | 2 | 17682.1 | 23855.7 | 8.88307e6 | 1.28422e7 | 1.89819e-4 |
| 0.5 | 3 | 15246.9 | 22236.4 | 9.10002e6 | 1.25648e7 | -2.64610e-4 |
| 0.5 | 4 | 14030.0 | 20796.9 | 8.74943e6 | 1.19185e7 | -4.29051e-4 |
| 0.5 | 5 | 13092.8 | 19391.7 | 8.07054e6 | 1.10062e7 | -3.55279e-4 |
| 0.5 | 6 | 12336.8 | 17953.2 | 7.13680e6 | 9.84211e6 | -2.62983e-5 |
| 0.5 | 7 | 12482.2 | 16459.8 | 6.62358e6 | 8.43603e6 | 2.73105e-3 |
| 0.5 | 8 | 12035.2 | 14894.6 | 5.54636e6 | 6.85743e6 | 2.45284e-3 |
| 0.5 | 9 | 10785.1 | 13209.9 | 4.21125e6 | 5.29361e6 | 3.57696e-4 |
| 0.5 | 10 | 9590.47 | 11220.6 | 2.71154e6 | 2.87940e6 | -3.25738e-3 |
| 0.5 | 11 | 11351.6 | 11378.0 | 3.92546e6 | 4.10278e6 | 3.50968e-3 |
| 0.5 | 12 | 9962.01 | 11710.9 | 3.01335e6 | 3.72672e6 | 3.79599e-5 |
| 0.5 | 13 | 9091.22 | 11766.4 | 3.10659e6 | 3.80477e6 | -5.97592e-5 |
| 0.5 | 14 | 8525.12 | 11780.8 | 2.95814e6 | 3.24977e6 | -3.92161e-5 |
| 1.0 | 0 | 59112.5 | 117314 | 7.37513e7 | 3.12373e8 | -7.50369e-4 |
| 1.0 | 1 | 61902.6 | 105079 | 5.59068e7 | 3.26899e8 | 1.72560e-3 |
| 1.0 | 2 | 44493.7 | 60590.5 | 2.31399e7 | 3.45801e7 | 8.95654e-4 |
| 1.0 | 3 | 37632.2 | 56605.6 | 2.34981e7 | 3.22598e7 | -1.05715e-3 |
| 1.0 | 4 | 34349.8 | 52972.3 | 2.28943e7 | 3.03814e7 | -2.41383e-3 |
| 1.0 | 5 | 31565.3 | 49281.2 | 2.11682e7 | 2.83041e7 | -2.66593e-3 |
| 1.0 | 6 | 29067.6 | 45416.9 | 1.85569e7 | 2.58083e7 | -1.19241e-3 |
| 1.0 | 7 | 28411.3 | 41504.5 | 1.69456e7 | 2.24573e7 | 7.39932e-3 |
| 1.0 | 8 | 29297.5 | 37750.5 | 1.44997e7 | 1.81113e7 | 6.65131e-3 |
| 1.0 | 9 | 25542.2 | 34020.5 | 1.04310e7 | 1.34280e7 | 1.97047e-3 |
| 1.0 | 10 | 22835.7 | 25866.0 | 6.31659e6 | 6.40395e6 | -7.44278e-3 |
| 1.0 | 11 | 28969.2 | 29814.9 | 9.66765e6 | 9.74724e6 | 8.89500e-3 |
| 1.0 | 12 | 25465.3 | 30282.5 | 7.10435e6 | 8.54906e6 | 3.73382e-4 |
| 1.0 | 13 | 23802.1 | 31265.3 | 7.50863e6 | 9.10614e6 | -1.91235e-4 |
| 1.0 | 14 | 22676.9 | 31920.2 | 7.04452e6 | 7.89663e6 | -2.61545e-4 |

## Source fix

The fix makes the passive law's gradient authoritative:

1. `mechanicalLaw` provides public virtual const accessors for its cell and
   face deformation gradients.
2. `electroMechanicalLaw` overrides both accessors and delegates to
   `passiveMechLawPtr_`.
3. Both active Cauchy-stress paths use those accessors rather than the wrapper's
   protected `F()` or `Ff()`.
4. `electroMechanicalLaw::setRestart()` and `updateTotalFields()` forward to the
   nested passive law.

The wrapper therefore does not lazily create a competing `F_` during active
stress evaluation. Total/updated Lagrangian and incremental/non-incremental
updates remain in the passive law's existing `updateF()` implementation. The
face path delegates in the same way as the volume path, and restart writes are
applied to the authoritative passive fields.

The residual-level invariant instrumentation was temporary and is not retained
in the production patch.

## Validation

### Residual-level invariant and Land Problem 3 L1

The corrected L1 case converged through time 1. At times 0, 0.01, 0.5, and 1:

```text
active/passive same object = true
max |F_active - (I + grad(D).T())| = 0
```

The time-1 step began at audit evaluation 71,472; the check remained active for
all subsequent residual evaluations, and no fail-fast condition was triggered.
The corrected run completed in 135.22 s, compared with 513.47 s for the
instrumented pre-fix run. The final reported maximum von Mises stress changed
from 76,622.6 Pa to 67,012.2 Pa, consistent with the active load being pushed
forward by the evolving deformation rather than a frozen startup field.

### Manufactured push-forward checks

Using `Ta = 60000 Pa` and `f0 = (1, 0, 0)`:

- zero displacement: maximum component error from `Ta*f0f0` was 0;
- affine deformation
  `F = ((1.2, 0.1, 0), (0, 0.9, 0.2), (0.05, 0, 1.1))`,
  `J = 1.189`, produced

  ```text
  ((72666.10597140454, 0, 3027.754415475189),
   (0,                 0, 0),
   (3027.754415475189, 0, 126.156433978133)) Pa
  ```

  with maximum component error `2.84217e-14 Pa`;
- rigid rotation by 0.63 rad about z had `J = 1` and zero objectivity error;
- fibre stretch with `F = diag(1.3, 0.8, 1.1)` had `J = 1.144` and
  `sigma_ff = Ta*1.3^2/J = 88636.36363636363 Pa`, with zero error and zero
  off-dyad components.

The Land 3 runtime fibre check reported zero cell and face Linf/RMS difference
between the active and passive fibre dyads.

### Passive Problems 1 and 2

Both passive cases set `activeTension` to zero, so `hasActiveStress()` prevents
entry into the changed active push-forward path. A copied Problem 1 mesh1 case
converged after the fix. A clean paired numerical comparison was not completed:
the available Problem 2 archived case is not self-contained after cleaning
because its initial `f0` field is absent. This does not affect the branch-level
active-path isolation, but a regenerated Problem 2 mesh1 case should be used
for a formal archived-result comparison.

## Expected Land Problem 3 effect

Land Problem 3 is the only one of the three cases with non-zero active tension.
The fix changes its active Cauchy stress from a push-forward by a frozen startup
gradient to a push-forward by the same current total deformation gradient used
by the passive Guccione response. The effect grows with load: the pre-fix
relative active-stress RMS error was approximately 0.45 at time 0.5 and 0.53 at
time 1. The corrected displacement, stress, residual, and torque response is
therefore expected to move materially, while zero-active-tension cases retain
the passive stress path.
