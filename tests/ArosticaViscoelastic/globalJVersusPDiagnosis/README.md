# Block-aware MFFD/FD audit

`ArosticaGlobalJVersusPDiagnosis` is an action-only diagnostic for the public
PETSc interfaces of the mixed total-Lagrangian solid model. It deliberately
does not form a complete finite-difference Jacobian, change a production
residual or material law, advance accepted history, or replay a production
solver step.

The diagnostic evaluates standalone-reconstructed search directions at
`1e-4, 1e-5, 1e-6, 1e-7, 1e-8` with central finite differences and compares
them with a reproducible PETSc MFFD action. It reports absolute error, relative
error, cosine, adjacent-epsilon change, local discrepancy cell, and the
assembled production preconditioner action. Complete, momentum-row, and
pressure-row actions are classified independently.

Thresholds are block-specific:

```text
tau_block(eps) = max(
    100*epsmach*Rscale_block/eps,
    100*Aref_block*epsmach,
    1e-12*Aref_block
)
```

`Rscale_block` is the maximum norm over the base and both perturbed residuals.
`Aref_D`, `Aref_p`, and `Aref_full` are separate median plateau references;
the pressure reference is never inferred from a momentum-dominated action.
The acceptance pair is `1e-6` and `1e-7`: both must be suitable, adjacent,
spatially consistent, and unchanged before the `1e-3` relative-error and
`0.999999` cosine gates are applied. Absolute errors and thresholds remain in
the report even when an action is below threshold.

Lifecycle validity is supplied explicitly with `-stateLabel`. Use
`pre-step`, `accepted-incomplete`, `accepted-complete`, or
`accepted-exact`. Pre-step and incomplete reconstructions are diagnostic only
and cannot be the primary accepted-state comparison.

The executable does not claim PETSc LGMRES vectors are Arnoldi basis vectors.
Every retained direction is labelled a standalone-reconstructed search
direction; residuals and solution updates retain those exact labels if a
secondary replay is added. The public model API does not expose exact Arnoldi
vectors, and that limitation does not invalidate the residual/MFFD/FD audit.

## Isolated regression

Run `./Alltest` from this directory after sourcing OpenFOAM v2312. The runner:

- snapshots the dirty and untracked source tree;
- records Git, OpenFOAM, PETSc, `WM_OPTIONS`, compiler environment, hashes,
  timestamps, Make-file target expansion, build logs, `otool` linkage and
  loader output;
- builds the model/dependency libraries, `solids4Foam`, and this diagnostic
  into a private runtime directory;
- proves that Land and Aróstica runtime selection loads that private model
  library and rejects the live platform, archived `1623b2f2...`, and other
  model-library paths;
- audits copied Land pre-step/accepted data, Aróstica pre-step data, and
  Controls A/B without writing their source cases; and
- archives every scratch case, dictionary, state file, epsilon/mode result,
  loader log, binary, checksum, and failure report before cleanup.

The permanent archive defaults to:

```text
/Volumes/OpenFoam/aaronmullen-hales-v2312/run/CardiacMechanics/ventricle/arostica/diagnostics/land_arostica_actual_state_mffd_fd_<UTC timestamp>/
```

Set `S4F_AUDIT_ARCHIVE_BASE`, `S4F_LAND_CASE`, or `S4F_AROSTICA_CASE` when the
local benchmark locations differ. The runner stops before the audit if the
private output paths, linkage, or loader proof cannot be made unambiguous.
