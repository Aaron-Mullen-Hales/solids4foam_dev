# Case-B failing-face lifecycle regression

This regression retains the exact first overflowing face from the archived
Case-B fibre-direction diagnostic. The data checker proves that the
law-owned and independently reconstructed deformation gradients agree and
that `exp(bs*sqr(I4s - 1))` is the first invalid operation. It does not clip,
replace, or evaluate the invalid tensor.

`Alltest` then runs the production diagnostic callback lifecycle on an
immutable copy of the archived 17,625-cell fixture:

1. base and repeated residual;
2. MFFD and repeated MFFD action;
3. central-FD plus;
4. complete accepted-state and boundary restoration; and
5. central-FD minus.

The pre-correction diagnostic traps during step 5. With accepted `sigma` and
physical `p` restored before each trial, every callback is finite and the
primary, deformation, history, boundary, and time fingerprints return to
their exact pre-callback values after each central-FD pair.

Set `S4F_CASEB_REPRO_CASE` only when the immutable evidence archive is stored
outside its default local path.
