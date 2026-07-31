# Aróstica compression-switch convention

The new law supports two explicit runtime options:

- compressionSwitch paperPiecewise
- compressionSwitch logistic

For either fibre or sheet invariant I4, the energy term is

    W_i = ai/(2 bi) chi(I4) [exp(bi (I4 - 1)^2) - 1].

The paperPiecewise option implements the paper convention exactly:

    chi(I4)  = I4,  if I4 > 1
             = 0,   otherwise

    dchi/dI4 = 1,   if I4 > 1
             = 0,   otherwise

At I4 = 1, the implementation uses the compression-side value zero. The law
is consequently non-smooth at the threshold.

The logistic option follows the pinned Finsberg reference implementation. The
switch is the logistic factor itself; it is not multiplied by I4:

    chi(I4)  = 1/(1 + exp[-compressionSwitchK (I4 - 1)])

    dchi/dI4 = compressionSwitchK chi (1 - chi).

The stress uses the derivative of the selected energy:

    dWi/dI4 = ai/(2 bi)
              [dchi/dI4 (exp(bi (I4 - 1)^2) - 1)
               + chi exp(bi (I4 - 1)^2) 2 bi (I4 - 1)].

compressionSwitchK is required and must be finite and positive for both modes.
Unknown switch names are fatal. The benchmark dictionary must choose the
intended mode explicitly; the law does not hide a default.

