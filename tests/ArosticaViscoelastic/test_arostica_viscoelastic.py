#!/usr/bin/env python3
"""Dependency-free checks for the Aróstica material viscosity contract."""

import math
import unittest


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def transpose(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def scale(a, value):
    return [[value * a[i][j] for j in range(3)] for i in range(3)]


def trace(a):
    return sum(a[i][i] for i in range(3))


def det(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


I = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def strain(f):
    c = matmul(transpose(f), f)
    return scale(add(c, scale(I, -1.0)), 0.5)


def viscous(f, e_old, e_old_old, eta, c0, c1, c2):
    edot = add(scale(strain(f), c0), add(scale(e_old, c1), scale(e_old_old, c2)))
    sigma = scale(matmul(matmul(f, scale(edot, eta)), transpose(f)), 1.0 / det(f))
    return edot, sigma


def max_error(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(3) for j in range(3))


class ArosticaViscoelasticTest(unittest.TestCase):
    def test_zero_and_constant_history(self):
        e0 = [[0.0] * 3 for _ in range(3)]
        f = I
        edot, sigma = viscous(f, e0, e0, 100.0, 50.0, -50.0, 0.0)
        self.assertLess(max_error(edot, e0), 1.0e-15)
        self.assertLess(max_error(sigma, e0), 1.0e-15)

        f = [[1.2, 0.1, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.1]]
        e = strain(f)
        _, sigma = viscous(f, e, e, 100.0, 50.0, -50.0, 0.0)
        self.assertLess(max_error(sigma, [[0.0] * 3 for _ in range(3)]), 1.0e-12)

    def test_backward_coefficients_and_linear_history(self):
        dt, dt0 = 0.02, 0.01
        c_t = 1.0 + dt / (dt + dt0)
        c_00 = dt * dt / (dt0 * (dt + dt0))
        c_0 = c_t + c_00
        c = (c_t / dt, -c_0 / dt, c_00 / dt)
        self.assertAlmostEqual(sum(c), 0.0, places=14)

        # E(t)=E0+q*t is recovered exactly by variable-step BDF2.
        q = scale(I, 0.3)
        e_nm1 = scale(q, -dt0)
        e_n = scale(q, 0.0)
        e_np1 = scale(q, dt)
        edot = add(scale(e_np1, c[0]), add(scale(e_n, c[1]), scale(e_nm1, c[2])))
        self.assertLess(max_error(edot, q), 1.0e-13)

    def test_first_step_fallback(self):
        f = [[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        edot, _ = viscous(f, [[0.0] * 3 for _ in range(3)], [[0.0] * 3 for _ in range(3)],
                           100.0, 1.0 / 0.02, -1.0 / 0.02, 0.0)
        self.assertAlmostEqual(edot[0][0], 5.25, places=14)

    def test_objectivity_and_mapping(self):
        f = [[1.15, 0.2, 0.0], [0.0, 0.95, 0.1], [0.0, 0.0, 1.05]]
        old = strain([[1.1, 0.0, 0.0], [0.0, 0.97, 0.0], [0.0, 0.0, 1.0]])
        old_old = strain(I)
        edot, sigma = viscous(f, old, old_old, 100.0, 4.0, -2.0, 0.5)
        angle = 0.37
        r = [[math.cos(angle), -math.sin(angle), 0.0],
             [math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0]]
        # A common spatial rotation leaves C, E and Edot unchanged; only the
        # pushed-forward Cauchy stress rotates.
        edot_r, sigma_r = viscous(matmul(r, f), old, old_old,
                                   100.0, 4.0, -2.0, 0.5)
        self.assertLess(max_error(edot_r, edot), 1.0e-12)
        self.assertLess(max_error(sigma_r, matmul(matmul(r, sigma), transpose(r))), 1.0e-10)
        self.assertGreaterEqual(sum(edot[i][j] * (100.0 * edot[i][j])
                                     for i in range(3) for j in range(3)), 0.0)

    def test_rejected_trial_does_not_change_history(self):
        e_old = strain([[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        e_old_old = strain([[1.05, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        f_a = [[1.2, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        f_b = [[0.8, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        a1 = viscous(f_a, e_old, e_old_old, 100.0, 2.0, -1.0, 0.0)[1]
        viscous(f_b, e_old, e_old_old, 100.0, 2.0, -1.0, 0.0)
        a2 = viscous(f_a, e_old, e_old_old, 100.0, 2.0, -1.0, 0.0)[1]
        self.assertLess(max_error(a1, a2), 1.0e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
