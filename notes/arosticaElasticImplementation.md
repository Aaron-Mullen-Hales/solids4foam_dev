# Aróstica Phase 1 elastic implementation

## Runtime classes

The Phase 1 runtime classes are:

- arosticaNonLinearGeometryTotalLagrangianTotalDisplacement,
  implemented by nonLinGeomTotalLagTotalDispArosticaSolid.
- ArosticaHolzapfelOgdenElastic.

The derived solid model only overrides
retainFullPassiveStressInMixedSplit() and requires solvePressure true. It
reuses the complete PETSc SNES/JFNK, mixed pressure, inertia, stabilisation,
follower-pressure, residual, Jacobian, and restart implementation in
nonLinGeomTotalLagTotalDispSolid.

Its stress contract is:

    sigma = sigmaNonVolumetricFull - p I.

The generic solid model remains unchanged and retains its established
dev(sigmaPassive) - p I path.

## Elastic energy

The law reads:

    I1bar = J^(-2/3) tr(C)
    I4f   = f0 . C . f0
    I4s   = s0 . C . s0
    I8fs  = f0 . C . s0

and implements:

    Psi_nonVol =
        a/(2b)   [exp(b (I1bar - 3)) - 1]
      + af/(2bf) chi(I4f) [exp(bf (I4f - 1)^2) - 1]
      + as/(2bs) chi(I4s) [exp(bs (I4s - 1)^2) - 1]
      + afs/(2bfs) [exp(bfs I8fs^2) - 1].

The eight material coefficients, finite bulkModulus, and impKcoeff are read
from the law dictionary. U(J) is deliberately excluded from the material
stress. The mixed solid model owns the pressure stress.

The constitutive kernel forms:

    S = 2 d(Psi_nonVol)/dC
    sigma = symm(F S F^T)/J.

The fibre, sheet, and fibre-sheet terms are not deviatorically projected.
Consequently their pushed-forward stress can have a non-zero trace under
dilation. The compression-switch formulas and exact energy derivative are
documented in notes/arosticaCompressionSwitch.md.

## Direction fields

The law requires these dimensionless fields:

    fibreField       f0
    sheetField       s0
    sheetNormalField n0
    faceFibreField       f0f
    faceSheetField       s0f
    faceSheetNormalField n0f

Cell and face fields are read independently. No interpolation, sign change,
cylindrical direction generation, sheet reconstruction, or normal
reconstruction is performed. Magnitudes, pairwise orthogonality, and finite
values are validated. The determinant is reported but its sign is not
constrained.

## Numerical stiffness and tangent

At construction, the law samples three local homogeneous shears in the
supplied f-s-n basis, reports the three moduli, and selects the largest
positive value as muRepresentative. The returned fields are:

    shearModulus = muRepresentative
    impK         = impKcoeff * muRepresentative
    bulkModulus  = bulkModulus

The implicit stiffness affects numerical preconditioning only. It does not
enter the physical stress.

materialTangentField() uses a deterministic centred finite difference of the
shared face constitutive kernel with respect to deformation-gradient
components. It contains no mixed pressure term and has no material state to
mutate. The current PETSc total-displacement model uses its structural
Jacobian/JFNK path rather than requiring this tangent for the nonlinear
residual.

## Example dictionary

    mechanical
    (
        myocardium
        {
            type electroMechanicalLaw;

            passiveMechanicalLaw
            {
                type ArosticaHolzapfelOgdenElastic;

                pressureDisplacement true;

                rho rho [1 -3 0 0 0 0 0] 1000;

                a   a   [1 -1 -2 0 0 0 0] 59.0;
                af  af  [1 -1 -2 0 0 0 0] 18472.0;
                as  as  [1 -1 -2 0 0 0 0] 2481.0;
                afs afs [1 -1 -2 0 0 0 0] 216.0;

                b   b   [0 0 0 0 0 0 0] 8.023;
                bf  bf  [0 0 0 0 0 0 0] 16.026;
                bs  bs  [0 0 0 0 0 0 0] 11.12;
                bfs bfs [0 0 0 0 0 0 0] 11.436;

                bulkModulus kappa [1 -1 -2 0 0 0 0] 1e6;
                impKcoeff 1.0;

                compressionSwitch paperPiecewise;
                compressionSwitchK 100;

                fibreField f0;
                sheetField s0;
                sheetNormalField n0;
                faceFibreField f0f;
                faceSheetField s0f;
                faceSheetNormalField n0f;
            }
        }
    );

This Phase 1 dictionary contains no viscosity, active-tension evolution,
spring-dashpot, or pressure-load entries.

## Tests

tests/ArosticaElastic/test_arostica_elastic.py is a dependency-free
material-point test suite covering identity, frame indifference, isochoric
shear, fibre/sheet extension, fibre-sheet shear, both switch modes,
pure-dilation full-versus-deviatoric stress, second-Piola/Cauchy mapping,
energy differentiation, and cell/face algebra consistency.

The source library build checks runtime compilation and registration of both
new classes. tests/ArosticaElastic/runtimeSmoke/Alltest instantiates the law
from the example dictionary and exercises both correct() paths on a uniform
block mesh, including the supplied face fields.
