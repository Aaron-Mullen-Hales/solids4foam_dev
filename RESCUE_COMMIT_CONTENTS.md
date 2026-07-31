# Rescue commit contents

This manifest defines the exact preservation boundary for branch
`rescue/pre-solid-model-cleanup-20260731`. The rescue starts from commit
`ebf998ff25945f611a44e5376656690c62768d32` and preserves the current relevant
working-tree state before the solid-model cleanup.

## Included tracked modifications

- `src/solids4FoamModels/Make/files.foamextend`
- `src/solids4FoamModels/Make/files.openfoam`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.C`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.H`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GuccioneElastic/GuccioneElastic.C`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GuccioneElastic/GuccioneElastic.H`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.C`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.H`
- `src/solids4FoamModels/numerics/foamPetscSnesHelper/foamPetscSnesHelper.C`
- `src/solids4FoamModels/numerics/foamPetscSnesHelper/foamPetscSnesHelper.H`
- `src/solids4FoamModels/numerics/leastSquaresS4fGrad/leastSquaresS4fGrad.C`
- `src/solids4FoamModels/numerics/leastSquaresS4fGrad/leastSquaresS4fGrad.H`
- `src/solids4FoamModels/numerics/leastSquaresS4fGrad/leastSquaresS4fGrads.C`
- `src/solids4FoamModels/numerics/stabilisationModels/stabilisationModel/stabilisationModel.H`
- `src/solids4FoamModels/solidModels/fvPatchFields/solidTraction/solidTractionFvPatchVectorField.C`
- `src/solids4FoamModels/solidModels/fvPatchFields/solidTraction/solidTractionFvPatchVectorField.H`
- `src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/nonLinGeomTotalLagTotalDispSolid.C`
- `src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid/nonLinGeomTotalLagTotalDispSolid.H`

## Included untracked source and headers

Every file recursively beneath these paths is included:

- `applications/utilities/stabilisationModelAudit/`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/ArosticaHolzapfelOgdenElastic/`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/ArosticaHolzapfelOgdenViscoelastic/`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic/`, except the generated build log listed below
- `src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowAffineExactStab/`
- `src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowHighPassStab/`
- `src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowZeroForceStab/`
- `src/solids4FoamModels/solidModels/fvPatchFields/arosticaNormalSpringDashpotTraction/`
- `src/solids4FoamModels/solidModels/fvPatchFields/arosticaSpringDashpotTraction/`
- `src/solids4FoamModels/solidModels/fvPatchFields/arosticaVectorSpringDashpotTraction/`
- `src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispArosticaSolid/`
- `src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispGultekinSolid/`

The local `Make/files`, `Make/options`, and audit dictionaries under those
paths are intentional build configuration or test inputs and are included.

## Included tests and case assets

Every file recursively beneath these paths is included:

- `tests/`
- `tutorials/solids/hyperelasticity/idealisedVentricle/0/`
- `tutorials/solids/hyperelasticity/idealisedVentricle/constant/`
- `tutorials/solids/hyperelasticity/idealisedVentricle/system/`

The following intentional case option files are included:

- `tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.hypre`
- `tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.ilu`
- `tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.lu`
- `tutorials/solids/hyperelasticity/idealisedVentricle/petscOptions.seg.hypre`

## Included documentation and audit records

Every file recursively beneath `notes/` is included, together with:

- `GultekinMixedSolid_implementation_report.md`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_bulk_modulus_audit.md`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_implementation_report.md`
- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_stress_interface_audit.md`
- `RESCUE_COMMIT_CONTENTS.md`

## Explicitly excluded generated files and results

These paths remain present in the original worktree but are not staged or
committed by the rescue operation:

- `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GultekinTwoFibreElastic_build.log` — generated build log
- `tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/0/T` — tracked generated result
- `tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/1/T` — tracked generated result
- `tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/2/T` — tracked generated result
- `tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/3/T` — tracked generated result
- `tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials/hotCylinderTemperatureField/4/T` — tracked generated result

No untracked generated library, executable, processor directory, simulation
time directory, or large solver log was present in the Phase 0 inventory.

## Classification summary

- Source and headers: included.
- Tests and test build configuration: included.
- Production build configuration: included.
- Dictionaries, scripts, and intentional case assets: included.
- Design notes, implementation reports, and audit records: included.
- Generated binaries, executables, logs, and simulation results: excluded.

