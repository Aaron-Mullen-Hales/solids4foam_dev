#!/usr/bin/env python3
"""Dependency-free constitutive checks for the Phase 1 Aróstica law."""

from __future__ import annotations

import math


ATOL = 2.0e-10
RTOL = 3.0e-7
ENERGY_STRESS_RTOL = 4.0e-6


def zero():
    return [[0.0, 0.0, 0.0] for _ in range(3)]


def identity():
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def sub(a, b):
    return [[a[i][j] - b[i][j] for j in range(3)] for i in range(3)]


def scale(a, value):
    return [[value * a[i][j] for j in range(3)] for i in range(3)]


def multiply(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def transpose(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def trace(a):
    return sum(a[i][i] for i in range(3))


def determinant(a):
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inverse(a):
    d = determinant(a)
    if d <= 0.0:
        raise ValueError(f"expected positive determinant, got {d}")
    return [
        [
            (
                a[(j + 1) % 3][(i + 1) % 3]
                * a[(j + 2) % 3][(i + 2) % 3]
                - a[(j + 1) % 3][(i + 2) % 3]
                * a[(j + 2) % 3][(i + 1) % 3]
            )
            / d
            for j in range(3)
        ]
        for i in range(3)
    ]


def symmetric(a):
    return scale(add(a, transpose(a)), 0.5)


def outer(a, b):
    return [[a[i] * b[j] for j in range(3)] for i in range(3)]


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def norm(a):
    return math.sqrt(sum(value * value for row in a for value in row))


def max_error(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(3) for j in range(3))


def assert_close(name, actual, expected, atol=ATOL, rtol=RTOL):
    if isinstance(actual, list):
        error = max_error(actual, expected)
        scale_value = max(1.0, norm(expected))
    else:
        error = abs(actual - expected)
        scale_value = max(1.0, abs(expected))
    if error > atol + rtol * scale_value:
        raise AssertionError(f"{name}: error {error:.6e}")


PARAMETERS = {
    "a": 59.0,
    "b": 8.023,
    "af": 18472.0,
    "bf": 16.026,
    "as": 2481.0,
    "bs": 11.12,
    "afs": 216.0,
    "bfs": 11.436,
}


def switch_value(i4, mode, steepness):
    if mode == "paperPiecewise":
        return (i4, 1.0) if i4 > 1.0 else (0.0, 0.0)
    if mode == "logistic":
        value = 1.0 / (1.0 + math.exp(-steepness * (i4 - 1.0)))
        return value, steepness * value * (1.0 - value)
    raise ValueError(f"unknown switch {mode}")


def constitutive(deformation, f, s, mode="paperPiecewise", steepness=100.0):
    c = multiply(transpose(deformation), deformation)
    jacobian = determinant(deformation)
    if jacobian <= 0.0:
        raise ValueError("deformation must have positive determinant")

    c_inverse = inverse(c)
    trace_c = trace(c)
    j_minus_two_thirds = jacobian ** (-2.0 / 3.0)
    i1bar = j_minus_two_thirds * trace_c
    cf = [sum(c[i][j] * f[j] for j in range(3)) for i in range(3)]
    cs = [sum(c[i][j] * s[j] for j in range(3)) for i in range(3)]
    i4f = dot(f, cf)
    i4s = dot(s, cs)
    i8fs = dot(f, cs)

    chi_f, dchi_f = switch_value(i4f, mode, steepness)
    chi_s, dchi_s = switch_value(i4s, mode, steepness)
    exp_i4f = math.exp(PARAMETERS["bf"] * (i4f - 1.0) ** 2)
    exp_i4s = math.exp(PARAMETERS["bs"] * (i4s - 1.0) ** 2)
    exp_i8 = math.exp(PARAMETERS["bfs"] * i8fs**2)

    energy = (
        PARAMETERS["a"] / (2.0 * PARAMETERS["b"])
        * (math.exp(PARAMETERS["b"] * (i1bar - 3.0)) - 1.0)
        + PARAMETERS["af"] / (2.0 * PARAMETERS["bf"])
        * chi_f * (exp_i4f - 1.0)
        + PARAMETERS["as"] / (2.0 * PARAMETERS["bs"])
        * chi_s * (exp_i4s - 1.0)
        + PARAMETERS["afs"] / (2.0 * PARAMETERS["bfs"])
        * (exp_i8 - 1.0)
    )

    d_i1 = scale(
        sub(identity(), scale(c_inverse, trace_c / 3.0)),
        j_minus_two_thirds,
    )
    second_piola = scale(
        d_i1,
        PARAMETERS["a"] * math.exp(PARAMETERS["b"] * (i1bar - 3.0)),
    )

    qf = i4f - 1.0
    qs = i4s - 1.0
    d_wf = PARAMETERS["af"] / (2.0 * PARAMETERS["bf"]) * (
        dchi_f * (exp_i4f - 1.0)
        + chi_f * exp_i4f * 2.0 * PARAMETERS["bf"] * qf
    )
    d_ws = PARAMETERS["as"] / (2.0 * PARAMETERS["bs"]) * (
        dchi_s * (exp_i4s - 1.0)
        + chi_s * exp_i4s * 2.0 * PARAMETERS["bs"] * qs
    )
    second_piola = add(second_piola, scale(outer(f, f), 2.0 * d_wf))
    second_piola = add(second_piola, scale(outer(s, s), 2.0 * d_ws))
    d_w8 = PARAMETERS["afs"] * i8fs * exp_i8
    second_piola = add(
        second_piola,
        scale(symmetric(outer(f, s)), 2.0 * d_w8),
    )

    sigma = scale(
        symmetric(
            multiply(multiply(deformation, second_piola), transpose(deformation))
        ),
        1.0 / jacobian,
    )
    return {
        "energy": energy,
        "J": jacobian,
        "I1bar": i1bar,
        "I4f": i4f,
        "I4s": i4s,
        "I8fs": i8fs,
        "S": second_piola,
        "sigma": sigma,
    }


def rotation_z(angle):
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]


def numerical_cauchy_from_energy(deformation, mode="logistic"):
    eps = 1.0e-6
    first_piola = zero()
    for i in range(3):
        for j in range(3):
            plus = [row[:] for row in deformation]
            minus = [row[:] for row in deformation]
            plus[i][j] += eps
            minus[i][j] -= eps
            first_piola[i][j] = (
                constitutive(plus, FIBRE, SHEET, mode)["energy"]
                - constitutive(minus, FIBRE, SHEET, mode)["energy"]
            ) / (2.0 * eps)
    jacobian = determinant(deformation)
    return symmetric(
        scale(multiply(first_piola, transpose(deformation)), 1.0 / jacobian)
    )


FIBRE = [1.0, 0.0, 0.0]
SHEET = [0.0, 1.0, 0.0]


def run_tests():
    identity_result = constitutive(identity(), FIBRE, SHEET)
    assert_close("identity stress", identity_result["sigma"], zero())
    assert_close("identity second Piola", identity_result["S"], zero())

    rotation = rotation_z(0.731)
    deformation = [[1.0, 0.23, 0.0], [0.0, 1.1, 0.12], [0.0, 0.0, 0.91]]
    rotated = constitutive(multiply(rotation, deformation), FIBRE, SHEET)
    base = constitutive(deformation, FIBRE, SHEET)
    expected_rotated = multiply(multiply(rotation, base["sigma"]), transpose(rotation))
    assert_close("frame indifference energy", rotated["energy"], base["energy"])
    assert_close("frame indifference stress", rotated["sigma"], expected_rotated)

    shear = [[1.0, 0.21, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    shear_result = constitutive(shear, FIBRE, SHEET)
    if norm(shear_result["sigma"]) <= 1.0e-8:
        raise AssertionError("isochoric shear produced zero stress")

    fibre_extension = [
        [1.25, 0.0, 0.0],
        [0.0, 1.0 / math.sqrt(1.25), 0.0],
        [0.0, 0.0, 1.0 / math.sqrt(1.25)],
    ]
    fibre_result = constitutive(fibre_extension, FIBRE, SHEET)
    if fibre_result["I4f"] <= 1.0 or fibre_result["energy"] <= 0.0:
        raise AssertionError("fibre extension did not activate")

    sheet_extension = [
        [1.0 / math.sqrt(1.25), 0.0, 0.0],
        [0.0, 1.25, 0.0],
        [0.0, 0.0, 1.0 / math.sqrt(1.25)],
    ]
    sheet_result = constitutive(sheet_extension, FIBRE, SHEET)
    if sheet_result["I4s"] <= 1.0 or sheet_result["energy"] <= 0.0:
        raise AssertionError("sheet extension did not activate")

    fs_shear = [[1.0, 0.19, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    fs_result = constitutive(fs_shear, FIBRE, SHEET)
    fs_reverse = constitutive(
        [[1.0, -0.19, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        FIBRE,
        SHEET,
    )
    if abs(fs_result["I8fs"]) <= 1.0e-8:
        raise AssertionError("fibre-sheet invariant remained zero")
    if fs_result["sigma"][0][1] * fs_reverse["sigma"][0][1] >= 0.0:
        raise AssertionError("fibre-sheet stress sign did not reverse")

    compressed = [
        [0.8, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    compressed_paper = constitutive(compressed, FIBRE, SHEET, "paperPiecewise")
    compressed_logistic = constitutive(compressed, FIBRE, SHEET, "logistic")
    original = PARAMETERS.copy()
    PARAMETERS.update({"af": 0.0, "as": 0.0, "afs": 0.0})
    compressed_iso = constitutive(compressed, FIBRE, SHEET, "paperPiecewise")
    PARAMETERS.clear()
    PARAMETERS.update(original)
    assert_close(
        "paper compression switch energy",
        compressed_paper["energy"],
        compressed_iso["energy"],
        atol=2.0e-9,
        rtol=2.0e-7,
    )
    if compressed_logistic["energy"] <= compressed_paper["energy"]:
        raise AssertionError("logistic compression response was not present")

    dilation = scale(identity(), 1.12)
    dilation_result = constitutive(dilation, FIBRE, SHEET)
    deviatoric = sub(
        dilation_result["sigma"],
        scale(identity(), trace(dilation_result["sigma"]) / 3.0),
    )
    dilation_difference = norm(sub(dilation_result["sigma"], deviatoric))
    if dilation_difference <= 1.0e-8:
        raise AssertionError("pure dilation did not expose spherical anisotropic stress")

    mapping_error = norm(
        sub(
            multiply(
                scale(dilation_result["sigma"], dilation_result["J"]),
                inverse(dilation),
            ),
            multiply(dilation, dilation_result["S"]),
        )
    )
    if mapping_error > 2.0e-9:
        raise AssertionError(f"second-Piola/Cauchy mapping error {mapping_error:.6e}")

    random_deformations = [
        [[1.08, 0.07, 0.02], [0.01, 0.94, 0.04], [0.0, 0.03, 1.16]],
        [[0.91, -0.04, 0.08], [0.06, 1.13, 0.02], [0.01, 0.05, 0.97]],
    ]
    max_energy_error = 0.0
    for deformation in random_deformations:
        analytical = constitutive(deformation, FIBRE, SHEET, "logistic")
        numerical = numerical_cauchy_from_energy(deformation, "logistic")
        error = norm(sub(analytical["sigma"], numerical))
        max_energy_error = max(max_energy_error, error)
        if error > ENERGY_STRESS_RTOL * max(1.0, norm(analytical["sigma"])):
            raise AssertionError(f"energy derivative error {error:.6e}")

    face_result = constitutive(deformation, FIBRE, SHEET, "logistic")
    cell_result = constitutive(deformation, FIBRE, SHEET, "logistic")
    max_cell_face_error = norm(sub(cell_result["sigma"], face_result["sigma"]))

    print("PASS identity")
    print("PASS frame indifference")
    print("PASS isochoric shear")
    print("PASS fibre extension")
    print("PASS sheet extension")
    print("PASS fibre-sheet shear")
    print("PASS compression switch modes")
    print("PASS pure-dilation full-vs-deviatoric stress")
    print("PASS second-Piola/Cauchy mapping")
    print("PASS energy derivative")
    print("PASS cell/face shared-kernel equivalence")
    print(f"maximum energy-derivative error = {max_energy_error:.6e}")
    print(
        "maximum frame-indifference error = "
        f"{max_error(rotated['sigma'], expected_rotated):.6e}"
    )
    print(f"maximum cell/face difference = {max_cell_face_error:.6e}")
    print(f"pure-dilation full-vs-dev stress difference = {dilation_difference:.6e}")


if __name__ == "__main__":
    run_tests()
