# GultekinMixedSolid focused tests

Run from the solids4foam source root:

    python3 tests/GultekinMixedSolid/test_gultekin_mixed_solid.py

The tests cover the scalar finite-K law, zero-stabilisation pressure closure,
the analytic and centred finite-difference derivative, retention of a
synthetic nonzero-trace fibre stress, preservation of the generic mixed-model
contract, single pressure insertion, the written Cauchy-stress dimensions,
runtime registration, and the unchanged unsplit material-law contract.

They are source and algebra tests; they do not run a solids4foam simulation.

