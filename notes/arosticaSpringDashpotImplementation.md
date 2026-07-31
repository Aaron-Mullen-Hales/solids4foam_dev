# Phase 2A implementation

## Runtime classes

The independently selectable patch-field runtime names are:

- `arosticaNormalSpringDashpotTraction`
- `arosticaVectorSpringDashpotTraction`

Both concrete classes derive from the shared, non-runtime
`arosticaSpringDashpotTractionFvPatchVectorField`, which derives from
`solidTractionFvPatchVectorField`. Thus the existing total-Lagrangian
`isA<solidTractionFvPatchVectorField>` path consumes them without a model
type check or copied residual implementation.

## Dictionary syntax

The canonical form is dimensioned:

```text
epicardium
{
    type                arosticaNormalSpringDashpotTraction;
    springCoefficient   springCoefficient [1 -2 -2 0 0 0 0] 1e8;
    dashpotCoefficient  dashpotCoefficient [1 -2 -1 0 0 0 0] 5e3;
    useUndeformedArea   true;
    writeDiagnostics    false;
    value               uniform (0 0 0);
}

base
{
    type                arosticaVectorSpringDashpotTraction;
    springCoefficient   springCoefficient [1 -2 -2 0 0 0 0] 1e5;
    dashpotCoefficient  dashpotCoefficient [1 -2 -1 0 0 0 0] 5e3;
    useUndeformedArea   true;
    value               uniform (0 0 0);
}
```

The accepted case's legacy scalar `alpha` and `beta` entries are supported
as aliases, so the case remains untouched. Canonical entries are checked for
`pressure/length` and `pressure*time/length`, respectively. Coefficients
must be finite and non-negative. An explicit `useUndeformedArea false` is a
fatal error; absent or true is forced to true.

The normal condition accepts the historical zero
`tangentialTraction` entry for compatibility, but rejects a non-zero value.
It supplies zero traction and pressure defaults when those generic
`solidTraction` entries are absent.

## Trial-state and velocity behaviour

`updateCoeffs()` computes the traction from the current trial `D` patch field
and `fvc::ddt(D)`. The active case uses `ddtSchemes { default backward; }`, so
the boundary condition follows the same backward scheme as the volume
momentum equation. The implementation is stateless: repeated SNES residuals,
matrix-free perturbations, rejected line searches, and restarts read the
current trial field and accepted old-time fields but never commit history.

The optional `writeDiagnostics true` entry reports maximum displacement,
velocity, and total traction; it does not alter the residual.

## Traction measure and solid-model coupling

The existing `enforceTractionBoundaries` method recognises the common base
class and multiplies `traction()` by the reference `magSf` because the
virtual `useUndeformedArea()` returns true. The pressure field is explicitly
set to zero. Current normals and current areas are not used to define either
support law.

## Scope and limitations

This phase implements only the two support tractions. It does not implement
material viscosity, active tension, pressure histories, Case A/B drivers, or
global Rayleigh damping. The only existing source class modified is the
shared `useUndeformedArea()` accessor, made virtual so this behavior can be
provided polymorphically; the generic `solidTraction` default and all other
existing behavior remain unchanged.

The automated tests cover the scalar formulas, dimensions, normal reversal,
reference-area force, repeated trial evaluation, restart-equivalent history,
mapping/clone invariants at the formula level, and runtime construction. The
OpenFOAM smoke test runs the one-cell Aróstica mixed solid through `evolve()`;
this calls the same production total-Lagrangian traction-enforcement path.

## Phase 2B Stage 0 production-path audit

The OpenFOAM v2312 call sequence is important here. `unpackSolution(x)`
updates the internal trial displacement. The production path constructs the
momentum surface force, and `solidTraction::evaluate()` calls the virtual
`updateCoeffs()` before applying the fixed-gradient correction.
`fixedGradientFvPatchField::evaluate()` then forms the boundary value from
the patch internal field and the updated gradient and resets the patch-field
updated flag.

The old Aróstica implementation read `D.boundaryField()[patchI]` before
calling the base `solidTractionFvPatchVectorField::updateCoeffs()`. That
ordering is not guaranteed to be the current PETSc trial value: after a
previous evaluation, the patch value can still be the value produced by the
previous fixed-gradient evaluation. Repeating one state therefore did not
prove trial safety.

The scoped Phase 2A repair is in the shared Aróstica spring-dashpot helper.
Before reading the displacement, it seeds the patch value from
`patchInternalField()`, which is the current trial predictor available at
this point in the fixed-gradient sequence. It then evaluates `fvc::ddt(D)`
and calls the base traction update. The generic `solidTraction` default
behavior is unchanged.

The Phase 2A `solidTraction` change making `useUndeformedArea()` virtual is
behavior-preserving for existing patches: the base default remains false,
the constructor, mapping, cloning, and runtime type remain unchanged, and
dictionary writing now emits the virtual value. The Aróstica override returns
true, so reference `magSf` is used for the support law.

## Trial-sequence verification

The dependency-free normal and full-vector tests use non-zero spring and
dashpot coefficients and distinct A/B displacement and velocity states. They
evaluate

```text
A
A -> B
B -> A
A -> B -> A
B -> A -> B
```

and compare each occurrence with its state-local reference. The scalar test
maximum difference is `0`. The OpenFOAM runtime smoke now also runs the
production PETSc `formResidual(x)` sequence for both the normal and
full-vector one-cell cases. For both cases, the maximum repeated-state
difference was `0`; the residual/traction sensitivity was `5.80032e6` for
the normal condition and `1.8299e4` for the full-vector condition. A
multi-patch harness remains desirable for the exact ventricular patch layout.

## Restart and velocity audit

The boundary condition remains stateless. Its velocity is constructed by
`fvc::ddt(D)`, so it uses the current trial `D`, the accepted OpenFOAM
old-time fields, and any older fields required by the configured scheme. It
does not write or advance accepted history.

The production PETSc harness tested both normal and full-vector patches with
non-zero spring and dashpot coefficients. The maximum repeated A/B trial
difference was `0`; first-evaluation residual/traction sensitivities were
`5.80032e6` (normal) and `1.8299e4` (vector). A five-step, two-process
zero-load restart on the ventricular scratch mesh completed successfully.
The compared `Ddot`, spring, dashpot, total traction, and complete written
state fields were byte-identical at every post-restart time. The original
campaign exposed and fixed a serialization defect in the common boundary
write method: dimensioned coefficients now retain a name token when written,
so a new process can read `springCoefficient` and `dashpotCoefficient`.

## `fvc::ddt(D)` cost

The source audit finds one full-volume `fvc::ddt(D)` construction per
Aróstica patch update. With one epicardial and one basal spring-dashpot patch,
the static count is therefore two constructions per production residual. The
dedicated probe measured the two-construction pair at below timer resolution
on the one-cell case and `0.0005 s` on the 17,625-cell ventricular scratch
mesh. The five-step zero-load production run took approximately `0.078 s`
per residual, so the pair is about `0.64%` of residual time. Allocation bytes
were not separately measurable. No shared trial-velocity cache was
introduced because this cost is small and the current implementation is
deterministic and does not cache across PETSc trial vectors.
