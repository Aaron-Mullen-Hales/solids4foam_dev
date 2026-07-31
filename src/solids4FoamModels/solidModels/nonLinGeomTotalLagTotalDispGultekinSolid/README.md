# Gültekin mixed finite-volume solid model

`nonLinGeomTotalLagTotalDispGultekinSolid` is an isolated, runtime-selected
solids4foam finite-volume mixed displacement-pressure model. It targets the
continuum equations used by the original Gültekin–Dal–Holzapfel Q1P0+WAS
benchmark. It does not implement Q1P0 finite elements, Hu–Washizu element
assembly, finite-element shape functions, or finite-element quadrature.

The model requires `solvePressure true` and one
`GultekinTwoFibreElastic` mechanical law. For the paper WAS model that law
must use `anisotropicSplit false`. The law returns passive spatial Cauchy
stress only:

    sigmaPassive = sigmaIso + sigmaFibre1 + sigmaFibre2

The matrix term is isochoric and deviatoric. The complete unsplit fibre dyads
are retained. The final stress is:

    sigma = sigmaPassive - p I

The original volumetric energy and mixed relation are:

    U(J) = K*(J - ln(J) - 1)

    -p/K -(J - 1)/J + S_p = 0

With zero pressure stabilisation this gives:

    -p I = K*(J - 1)/J I

Pressure and momentum stabilisation remain numerical finite-volume terms;
they are not part of the paper's constitutive energy. Matching the continuum
equations does not guarantee identical coarse-mesh FE and FV results.

Select the runtime model in `constant/solidProperties`:

    solidModel nonLinGeomTotalLagTotalDispGultekinSolid;

    nonLinGeomTotalLagTotalDispGultekinSolidCoeffs
    {
        solutionAlgorithm PETScSNES;
        solvePressure true;

        // Existing numerical settings belong here unchanged.

        writeGultekinMixedDiagnostics false;
    }

When enabled, `writeGultekinMixedDiagnostics` writes:

- `sigmaPassiveFull`: volSymmTensorField, Pa;
- `fibreMeanStress`: volScalarField, Pa;
- `sigmaVolumetricFromP`: volSymmTensorField, Pa;
- `gultekinConstraintPhysical`: volScalarField, dimensionless;
- `pressureStabilisationContribution`: volScalarField, dimensionless;
- `completePressureResidual`: volScalarField, dimensionless.

