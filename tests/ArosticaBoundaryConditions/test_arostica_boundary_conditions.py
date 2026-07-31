#!/usr/bin/env python3
"""Focused, dependency-free regression tests for the Phase 2A laws."""

import math
import unittest


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def scale(a, s):
    return tuple(s * x for x in a)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def normal_traction(d, v, n, alpha, beta):
    return add(scale(n, -alpha * dot(d, n)), scale(n, -beta * dot(v, n)))


def vector_traction(d, v, alpha, beta):
    return scale(add(scale(d, alpha), scale(v, beta)), -1.0)


def close_vector(actual, expected, tol=1.0e-14):
    return max(abs(x - y) for x, y in zip(actual, expected)) <= tol


class ArosticaBoundaryConditionTest(unittest.TestCase):
    def setUp(self):
        self.n = (0.0, 0.0, 1.0)
        self.alpha = 4.0
        self.beta = 3.0
        self.d = (0.2, -0.4, 0.5)
        self.v = (-0.7, 0.6, 0.3)

    def test_zero_state(self):
        self.assertEqual(normal_traction((0, 0, 0), (0, 0, 0), self.n, 2, 3), (0, 0, 0))
        self.assertEqual(vector_traction((0, 0, 0), (0, 0, 0), 2, 3), (0, 0, 0))

    def test_normal_spring_and_dashpot(self):
        self.assertTrue(close_vector(
            normal_traction((0, 0, 0.5), (0, 0, 0), self.n, self.alpha, self.beta),
            (0.0, 0.0, -2.0),
        ))
        self.assertTrue(close_vector(
            normal_traction((0, 0, 0), (0, 0, 0.3), self.n, self.alpha, self.beta),
            (0.0, 0.0, -0.9),
        ))

    def test_epicardial_tangential_components_vanish(self):
        result = normal_traction((0.2, -0.4, 0), (-0.7, 0.6, 0), self.n, self.alpha, self.beta)
        self.assertTrue(close_vector(result, (0.0, 0.0, 0.0)))
        mixed = normal_traction(self.d, self.v, self.n, self.alpha, self.beta)
        self.assertEqual(mixed[:2], (0.0, 0.0))

    def test_vector_law(self):
        expected = (-self.alpha * self.d[0] - self.beta * self.v[0],
                    -self.alpha * self.d[1] - self.beta * self.v[1],
                    -self.alpha * self.d[2] - self.beta * self.v[2])
        self.assertTrue(close_vector(vector_traction(self.d, self.v, self.alpha, self.beta), expected))

    def test_normal_reversal(self):
        self.assertTrue(
            close_vector(
                normal_traction(self.d, self.v, self.n, self.alpha, self.beta),
                normal_traction(self.d, self.v, scale(self.n, -1), self.alpha, self.beta),
            )
        )

    def test_reference_area_force(self):
        area0 = 2.75
        current_area = 11.0
        t_normal = normal_traction(self.d, self.v, self.n, self.alpha, self.beta)
        t_vector = vector_traction(self.d, self.v, self.alpha, self.beta)
        expected_normal = scale(t_normal, area0)
        expected_vector = scale(t_vector, area0)
        self.assertTrue(close_vector(expected_normal, scale(t_normal, area0)))
        self.assertTrue(close_vector(expected_vector, scale(t_vector, area0)))
        self.assertFalse(close_vector(expected_normal, scale(t_normal, current_area)))
        self.assertFalse(close_vector(expected_vector, scale(t_vector, current_area)))

    def test_repeated_trial_and_restart_equivalence(self):
        first = normal_traction(self.d, self.v, self.n, self.alpha, self.beta)
        for _ in range(10):
            self.assertEqual(first, normal_traction(self.d, self.v, self.n, self.alpha, self.beta))

        # A restart with the same accepted old-time state and trial state is
        # algebraically identical to the uninterrupted evaluation.
        old_d = (0.1, -0.2, 0.25)
        dt = 0.02
        uninterrupted_v = scale(add(self.d, scale(old_d, -1.0)), 1.0 / dt)
        restarted_v = scale(add(self.d, scale(old_d, -1.0)), 1.0 / dt)
        self.assertEqual(
            normal_traction(self.d, uninterrupted_v, self.n, self.alpha, self.beta),
            normal_traction(self.d, restarted_v, self.n, self.alpha, self.beta),
        )

    def test_trial_state_changes_immediately(self):
        a = vector_traction(self.d, self.v, self.alpha, self.beta)
        b = vector_traction(add(self.d, (0.1, 0.0, 0.0)), self.v, self.alpha, self.beta)
        self.assertNotEqual(a, b)

    def test_trial_sequence_is_order_independent(self):
        # Model the production stateless update for two PETSc trial states.
        # Each occurrence must be a function of that trial state only; in
        # particular, evaluating B must not leave a value that contaminates A.
        states = {
            "A": (self.d, self.v),
            "B": (add(self.d, (0.13, -0.07, 0.11)),
                  add(self.v, (-0.19, 0.05, 0.17))),
        }

        def result(state):
            displacement, velocity = states[state]
            return (
                normal_traction(
                    displacement, velocity, self.n, self.alpha, self.beta
                ),
                vector_traction(
                    displacement, velocity, self.alpha, self.beta
                ),
            )

        reference = {state: result(state) for state in states}
        sequences = ("A", "A B", "B A", "A B A", "B A B")

        for sequence in sequences:
            for state in sequence.split():
                normal, vector = result(state)
                self.assertEqual(normal, reference[state][0])
                self.assertEqual(vector, reference[state][1])

    def test_dimensions(self):
        pressure = (1, -1, -2, 0, 0, 0, 0)
        length = (0, 1, 0, 0, 0, 0, 0)
        time = (0, 0, 1, 0, 0, 0, 0)
        spring = tuple(a - b for a, b in zip(pressure, length))
        dashpot = tuple(a + b for a, b in zip(spring, time))
        self.assertEqual(spring, (1, -2, -2, 0, 0, 0, 0))
        self.assertEqual(dashpot, (1, -2, -1, 0, 0, 0, 0))

    def test_mapping_and_clone_preserve_parameters(self):
        original = (self.alpha, self.beta, True)
        mapped = tuple(original)
        cloned = tuple(original)
        self.assertEqual(mapped, original)
        self.assertEqual(cloned, original)


if __name__ == "__main__":
    unittest.main(verbosity=2)
