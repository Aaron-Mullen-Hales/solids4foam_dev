# Aróstica material viscosity (Phase 2B)

## Constitutive definition

`ArosticaHolzapfelOgdenViscoelastic` extends the existing
`ArosticaHolzapfelOgdenElastic` runtime law.  The passive kernel, supplied
fibre/sheet directions, invariant evaluation, Cauchy push-forward, field
validation, representative shear modulus, and numerical tangent helpers are
shared with the elastic implementation.

The benchmark viscous energy and stresses are

```text
E             = 0.5*(C - I),       C = F.T()*F
PsiViscous    = eta/2 * tr(Edot^2)
SViscous      = eta*Edot
sigmaViscous = symm(F*SViscous*F.T())/J
```

The returned material stress is the complete non-volumetric stress

```text
sigmaNonVolumetricFull = sigmaPassiveFull + sigmaViscous
```

No `dev()` operation is applied to the viscous stress.  The mixed Aróstica
solid model adds `-p*I`; pressure and volumetric penalty terms are not part of
this law.

The runtime name is exactly:

```text
ArosticaHolzapfelOgdenViscoelastic
```

`eta` is a dimensioned scalar with dimensions pressure*time:

```text
eta eta [1 -1 -1 0 0 0 0] 100;
```

It is validated as finite and non-negative.

## History and trial lifecycle

The law owns explicit, registered strain fields.  For a law instance named
`arosticaViscoSmoke`, the fields are:

```text
ArosticaE_arosticaViscoSmoke
ArosticaEOld_arosticaViscoSmoke
ArosticaEOldOld_arosticaViscoSmoke
ArosticaEf_arosticaViscoSmoke
ArosticaEfOld_arosticaViscoSmoke
ArosticaEfOldOld_arosticaViscoSmoke
```

The cell fields use the current `F()` field and the face fields use the
supplied direct `Ff()` field.  On each residual/stress call, only the current
`E`/`Ef` fields and derived trial diagnostics are overwritten.  The accepted
`EOld`, `EOldOld`, `EfOld`, and `EfOldOld` fields are changed only by
`updateTotalFields()`, which is called by the mechanical-model timestep
acceptance path.  Thus a rejected PETSc or line-search trial cannot become
history.

At construction, written history fields are read from the current time
directory.  If no explicit old history exists, the current strain field is
used as the initial accepted state.  On the first accepted commit,

```text
EOldOld <- EOld
EOld    <- E
EfOldOld <- EfOld
EfOld    <- Ef
```

The explicit fields are `AUTO_WRITE` and are included by `setRestart()`.  A
restart therefore reads the accepted current/old/old-old cell and face
strains instead of fabricating old-time levels inside a residual evaluation.
`updateTotalFields()` records the OpenFOAM `timeIndex` of its last shift and
returns immediately if called again at that same index.  This prevents a
coupling or solver finalisation path from shifting accepted history twice.

## Time discretisation

The default `ddtSchemes` entry is read from `fvSchemes`.  `backward` and
`Euler` are supported; other schemes fail clearly during construction.

For the first step, and whenever the second accepted history is unavailable,
the startup rule is first-order Euler:

```text
Edot = (E - EOld)/deltaT
```

For later backward steps, with `dt = deltaT` and `dt0 = deltaT0`, the
OpenFOAM v2312 variable-step coefficients are used:

```text
ct   = 1 + dt/(dt + dt0)
ct00 = dt*dt/(dt0*(dt + dt0))
ct0  = ct + ct00

Edot = (ct*E - ct0*EOld + ct00*EOldOld)/dt
```

The same coefficients are applied independently to cell `E` and direct face
`Ef`; face viscous stress is never obtained by interpolating cell stress.

## Cell and face algorithms

For each cell, the law computes `E`, `Edot`, `Sviscous`, and
`sigmaViscous` from the current cell `F()`.  For every internal and boundary
face it computes the corresponding quantities from the supplied `Ff()`.
The complete viscous stress is then added to the passive stress returned by
the shared elastic kernel.

Optional diagnostics are controlled by `writeDiagnostics`:

```text
writeDiagnostics true;
```

They are `Edot`, `Sviscous`, and `sigmaViscous`, plus the direct-face
equivalents `Edotf`, `Sviscousf`, and `sigmaViscousf`.  Diagnostics do not
participate in the stress calculation.

## Tangent and `impK`

`materialTangentField()` uses the existing central finite-difference passive
kernel and adds the current-trial derivative of the viscous stress.  It uses
accepted histories as fixed data and contains no pressure contribution.
The PETSc production path may use this material tangent for constitutive
linearisation; matrix-free residual differences include the same viscous
trial response directly.

The scalar approximate stiffness returned by `impK()` is inherited from the
elastic law and therefore retains the passive representative shear modulus.
No viscous `impK` contribution is enabled in this phase.

## Example

```text
mechanicalLaw
{
    type ArosticaHolzapfelOgdenViscoelastic;
    eta eta [1 -1 -1 0 0 0 0] 100;
    writeDiagnostics false;

    // The remaining entries are the ArosticaHolzapfelOgdenElastic entries.
    // f0, s0, n0 and their face fields are supplied by the case.
}
```

## Tests and limitations

The dependency-free constitutive tests cover zero and constant history,
linear strain history, first-step fallback, variable-step backward
coefficients, rotation objectivity, stress mapping, positive dissipation,
and rejected-trial history immutability.  The OpenFOAM runtime smoke now
also checks the first-step fallback and five unequal timestep values against
independent coefficients; the maximum first-step and variable-deltaT Edot
errors were both `0` in the one-cell production field test.

The PETSc production boundary harness selects this law with `solvePressure
true` and the logistic compression switch.  For both the normal and vector
spring-dashpot patches, `xA -> xB -> xA -> xB` gave repeated residual,
traction, and trial-state differences of `0`, while the residual/traction
sensitivity was `5.80032e6` and `1.8299e4`, respectively.  The MFFD-style
sweep used epsilons `1e-7` through `1e-2`; the maximum difference from the
`1e-7` reference was `2.13393e-9` for the normal patch and `7.49989e-12`
for the vector patch.  Accepted histories remained unchanged (`0`) during
all perturbations.

A genuine two-process zero-load monoventricle restart was also run for five
timesteps.  At times `0.003`, `0.004`, and `0.005`, `D`, `p`, `F`, direct-face
`Ff`, all six history fields, `Edot`, viscous stresses, and total stress were
byte-identical between uninterrupted and restarted runs.  The zero-load
scratch run selected both support boundary types and retained `min(J)=1`.

The constitutive tests use a roundoff-scale tolerance (typically below
`1e-13` in nondimensional tensor quantities).  No active tension, activation
ODE, pressure table, fibre-generator, or accepted benchmark-case changes are
included.  A controlled non-zero pressure/traction perturbation on the
17,625-cell scratch mesh currently reaches PETSc but diverges in the linear
solve; this is recorded as the remaining two-patch integration blocker, not
silently treated as a passing loading test.

The focused loaded-solve audit found the source-level residual finite and
deterministic at the initial state (`J=1`, zero trial strain and pressure).
The failure is the outer Krylov/preconditioner path, not viscous stress or
history mutation. See `notes/arosticaLoadedLinearSolveDiagnosis.md` for the
exact SNES/KSP hierarchy, direct assembled solve, isolation matrix, and
operator scale audit.

## Direct-face viscous `P_DD` contribution

The bounded Phase 2C contribution is selected with
`preconditionerViscousTangent true` in the solid-model dictionary. Its
default is `false`. It is rejected explicitly for `interpolatedCell`; the
direct-face tangent is not reused for that residual path.

The mechanical-law hook returns nine face tensors, one for each current-trial
component of `Ff`. With accepted face histories held fixed,

```text
Sviscous = eta*(c0*E + c1*EOld + c2*EOldOld)
Pviscous = F*Sviscous
dPviscous = dF*Sviscous + F*eta*c0*dE
```

The equality `Pviscous = J*sigmaViscous*F^-T` is exact because
`Sviscous` is symmetric. The compact matrix inserts the full displacement
component coupling using the existing internal-face owner/neighbour force
stencil. It changes `P` only from `formJacobian()` and never changes the
residual or accepted histories. Boundary spring/dashpot blocks remain
separate.

On the direct 17,625-cell audit state, the viscous-only switch combination
(corrected `P_Dp`, viscous tangent enabled, passive and boundary diagnostic
tangents disabled) changed the localised mode to:

```text
||J_DD v|| = 1769.39643664
||P_DD v|| = 1941.50918199
relative error = 0.170972333645
cosine = 0.990991433726
```

The complete decomposition independently measured the viscous action as
`1712.11198693`, with cosine `0.999753654407` against the total. Enabling
the existing passive and boundary diagnostic blocks as well gives relative
error `0.250961909622` and cosine `0.991802171161`; the switch is useful but
is not a complete displacement preconditioner.

At `eta = 0`, the assembled viscous tangent has zero norm on all audited
faces. Repeated residual difference and accepted-history checksum change
remain zero. The loaded fieldsplit case still reaches `DIVERGED_ITS` at 1000
outer iterations, so this contribution is accepted as an action-level
improvement while further `P_DD` work remains necessary.
