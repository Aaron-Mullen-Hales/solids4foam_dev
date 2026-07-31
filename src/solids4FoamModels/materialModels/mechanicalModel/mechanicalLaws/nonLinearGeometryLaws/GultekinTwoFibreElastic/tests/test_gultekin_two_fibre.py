#!/usr/bin/env python3
"""Independent material-point checks for GultekinTwoFibreElastic."""

from __future__ import annotations

import math
import sys

import numpy as np


MU = 10.0
K1 = 500.0
K2 = 2.0
ATOL = 2.0e-10
RTOL = 2.0e-9
ENERGY_STRESS_RTOL = 3.0e-6

ROOT_TWO = math.sqrt(2.0)
M1 = np.array([1.0 / ROOT_TWO, 1.0 / ROOT_TWO, 0.0])
M2 = np.array([-1.0 / ROOT_TWO, 1.0 / ROOT_TWO, 0.0])
IDENTITY = np.eye(3)


def positive_part(value: float) -> float:
    return max(value, 0.0)


def deviator(tensor: np.ndarray) -> np.ndarray:
    return tensor - np.trace(tensor) * IDENTITY / 3.0


def assert_close(
    name: str,
    actual: np.ndarray | float,
    expected: np.ndarray | float,
    atol: float = ATOL,
    rtol: float = RTOL,
) -> None:
    if not np.allclose(actual, expected, atol=atol, rtol=rtol):
        difference = np.max(np.abs(np.asarray(actual) - np.asarray(expected)))
        raise AssertionError(f"{name}: maximum error {difference:.6e}")


def invariant_data(
    deformation: np.ndarray,
    direction: np.ndarray,
    split: bool,
) -> tuple[float, float, np.ndarray]:
    jacobian = float(np.linalg.det(deformation))
    right_cauchy_green = deformation.T @ deformation
    raw_invariant = float(direction @ right_cauchy_green @ direction)
    invariant = jacobian ** (-2.0 / 3.0) * raw_invariant if split else raw_invariant
    current_direction = deformation @ direction
    return raw_invariant, invariant, current_direction


def strain_energy(
    deformation: np.ndarray,
    split: bool = False,
    tension_only: bool = False,
    use_second_family: bool = True,
) -> float:
    jacobian = float(np.linalg.det(deformation))
    if not np.isfinite(jacobian) or jacobian <= 0.0:
        raise ValueError(f"invalid J={jacobian}")

    right_cauchy_green = deformation.T @ deformation
    first_invariant_bar = jacobian ** (-2.0 / 3.0) * np.trace(right_cauchy_green)
    energy = 0.5 * MU * (first_invariant_bar - 3.0)

    directions = (M1, M2) if use_second_family else (M1,)
    for direction in directions:
        _, invariant, _ = invariant_data(deformation, direction, split)
        fibre_strain = invariant - 1.0
        if tension_only:
            fibre_strain = positive_part(fibre_strain)
        energy += K1 / (2.0 * K2) * (
            math.exp(K2 * fibre_strain * fibre_strain) - 1.0
        )

    return energy


def analytical_stress(
    deformation: np.ndarray,
    split: bool = False,
    tension_only: bool = False,
    use_second_family: bool = True,
) -> dict[str, np.ndarray | float]:
    jacobian = float(np.linalg.det(deformation))
    if not np.isfinite(jacobian) or jacobian <= 0.0:
        raise ValueError(f"invalid J={jacobian}")

    left_cauchy_green = deformation @ deformation.T
    jacobian_minus_two_thirds = jacobian ** (-2.0 / 3.0)
    matrix_stress = (
        MU
        / jacobian
        * deviator(jacobian_minus_two_thirds * left_cauchy_green)
    )

    family_stresses: list[np.ndarray] = []
    invariants: list[float] = []
    directions = (M1, M2) if use_second_family else (M1,)
    for direction in directions:
        _, invariant, current_direction = invariant_data(
            deformation, direction, split
        )
        invariants.append(invariant)
        fibre_strain = invariant - 1.0
        if tension_only:
            fibre_strain = positive_part(fibre_strain)

        current_dyad = np.outer(current_direction, current_direction)
        coefficient = (
            2.0
            * K1
            * fibre_strain
            * math.exp(K2 * fibre_strain * fibre_strain)
            / jacobian
        )
        if split:
            family_stresses.append(
                coefficient
                * jacobian_minus_two_thirds
                * deviator(current_dyad)
            )
        else:
            family_stresses.append(coefficient * current_dyad)

    if not use_second_family:
        family_stresses.append(np.zeros((3, 3)))
        invariants.append(0.0)

    passive = matrix_stress + family_stresses[0] + family_stresses[1]
    return {
        "J": jacobian,
        "I4": invariants[0],
        "I6": invariants[1],
        "sigma_iso": matrix_stress,
        "sigma_fibre1": family_stresses[0],
        "sigma_fibre2": family_stresses[1],
        "sigma_passive": passive,
    }


def numerical_cauchy_from_energy(
    deformation: np.ndarray,
    split: bool,
    step: float = 1.0e-6,
) -> np.ndarray:
    first_piola = np.zeros((3, 3))
    for row in range(3):
        for column in range(3):
            perturbation = np.zeros((3, 3))
            perturbation[row, column] = step
            energy_plus = strain_energy(deformation + perturbation, split=split)
            energy_minus = strain_energy(deformation - perturbation, split=split)
            first_piola[row, column] = (
                energy_plus - energy_minus
            ) / (2.0 * step)

    jacobian = float(np.linalg.det(deformation))
    cauchy = first_piola @ deformation.T / jacobian
    return 0.5 * (cauchy + cauchy.T)


def isochoric_extension_along(direction: np.ndarray, stretch: float) -> np.ndarray:
    direction_dyad = np.outer(direction, direction)
    return (
        stretch * direction_dyad
        + stretch ** (-0.5) * (IDENTITY - direction_dyad)
    )


def rotation_z(angle: float) -> np.ndarray:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return np.array(
        [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )


def effective_implicit_stiffness(
    specified_stiffness: float,
    imp_k_coefficient: float,
) -> float:
    """Return the numerical solver stiffness; it is not constitutive data."""
    return specified_stiffness * imp_k_coefficient


def mixed_pressure_compressibility_term(
    pressure: float,
    bulk_modulus: float,
) -> float:
    """Return the physical -p/K term owned by the mixed pressure equation."""
    if not math.isfinite(bulk_modulus) or bulk_modulus <= 0.0:
        raise ValueError("bulk_modulus must be finite and positive")
    return -pressure / bulk_modulus


def run_tests() -> None:
    tests_run = 0

    identity = analytical_stress(IDENTITY)
    assert_close("A.J", identity["J"], 1.0)
    assert_close("A.I4", identity["I4"], 1.0)
    assert_close("A.I6", identity["I6"], 1.0)
    assert_close("A.passive stress", identity["sigma_passive"], np.zeros((3, 3)))
    tests_run += 1
    print("PASS A: undeformed configuration")

    rotation = rotation_z(0.731)
    rotated = analytical_stress(rotation)
    assert_close("B.rotation invariants", [rotated["I4"], rotated["I6"]], [1.0, 1.0])
    assert_close("B.rotation stress", rotated["sigma_passive"], np.zeros((3, 3)))
    shear = IDENTITY.copy()
    shear[0, 1] = 0.23
    rotated_shear = analytical_stress(rotation @ shear)
    base_shear = analytical_stress(shear)
    assert_close(
        "B.objectivity",
        rotated_shear["sigma_passive"],
        rotation @ base_shear["sigma_passive"] @ rotation.T,
    )
    tests_run += 1
    print("PASS B: rigid rotation and objectivity")

    stretch = 1.2
    uniaxial = np.diag([stretch, stretch ** (-0.5), stretch ** (-0.5)])
    uniaxial_result = analytical_stress(uniaxial)
    assert_close("C.J", uniaxial_result["J"], 1.0)
    if not np.all(np.isfinite(uniaxial_result["sigma_passive"])):
        raise AssertionError("C: non-finite stress")
    tests_run += 1
    print("PASS C: isochoric uniaxial extension")

    fibre_extension = isochoric_extension_along(M1, 1.18)
    fibre_result = analytical_stress(fibre_extension)
    assert_close("D.I4", fibre_result["I4"], 1.18**2)
    if np.linalg.norm(fibre_result["sigma_fibre1"]) <= 0.0:
        raise AssertionError("D: family-1 stress was not activated")
    tests_run += 1
    print("PASS D: extension along fibre family 1")

    transverse_extension = np.diag([1.12 ** (-0.5), 1.12 ** (-0.5), 1.12])
    transverse_result = analytical_stress(transverse_extension)
    if not (transverse_result["I4"] < 1.0 and transverse_result["I6"] < 1.0):
        raise AssertionError("E: expected compression of both fibre invariants")
    if np.linalg.norm(transverse_result["sigma_fibre1"]) <= 0.0:
        raise AssertionError("E: default compressed-fibre response is missing")
    tests_run += 1
    print("PASS E: extension transverse to both fibre families")

    simple_shear = IDENTITY.copy()
    simple_shear[0, 1] = 0.35
    shear_result = analytical_stress(simple_shear)
    assert_close("F.J", shear_result["J"], 1.0)
    if not np.all(np.isfinite(shear_result["sigma_passive"])):
        raise AssertionError("F: non-finite simple-shear stress")
    tests_run += 1
    print("PASS F: simple shear")

    def swapped_total_stress(deformation: np.ndarray) -> np.ndarray:
        original_m1 = M1.copy()
        original_m2 = M2.copy()
        jacobian = np.linalg.det(deformation)
        left_cauchy_green = deformation @ deformation.T
        total = MU / jacobian * deviator(
            jacobian ** (-2.0 / 3.0) * left_cauchy_green
        )
        for direction in (original_m2, original_m1):
            _, invariant, current = invariant_data(deformation, direction, False)
            strain = invariant - 1.0
            total += (
                2.0
                * K1
                * strain
                * math.exp(K2 * strain * strain)
                / jacobian
                * np.outer(current, current)
            )
        return total

    assert_close(
        "G.interchange",
        shear_result["sigma_passive"],
        swapped_total_stress(simple_shear),
    )
    tests_run += 1
    print("PASS G: fibre-family interchange")

    one_family = analytical_stress(simple_shear, use_second_family=False)
    assert_close("H.family-2 stress", one_family["sigma_fibre2"], np.zeros((3, 3)))
    tests_run += 1
    print("PASS H: disabled second family")

    split_isochoric = analytical_stress(simple_shear, split=True)
    unsplit_isochoric = analytical_stress(simple_shear, split=False)
    assert_close(
        "I.J=1 energy",
        strain_energy(simple_shear, split=True),
        strain_energy(simple_shear, split=False),
    )
    assert_close(
        "I.J=1 pressure-equivalent stress",
        deviator(split_isochoric["sigma_passive"]),
        deviator(unsplit_isochoric["sigma_passive"]),
    )
    dilated_shear = 1.07 * simple_shear
    split_dilated = analytical_stress(dilated_shear, split=True)
    unsplit_dilated = analytical_stress(dilated_shear, split=False)
    if np.allclose(
        deviator(split_dilated["sigma_passive"]),
        deviator(unsplit_dilated["sigma_passive"]),
        atol=1.0e-7,
        rtol=1.0e-7,
    ):
        raise AssertionError("I: expected split/unsplit difference for J != 1")
    tests_run += 1
    print("PASS I: split/unsplit pressure-equivalence at J=1 and difference at J!=1")

    transverse_tension_only = analytical_stress(
        transverse_extension, tension_only=True
    )
    assert_close(
        "J.compressed family 1",
        transverse_tension_only["sigma_fibre1"],
        np.zeros((3, 3)),
    )
    assert_close(
        "J.compressed family 2",
        transverse_tension_only["sigma_fibre2"],
        np.zeros((3, 3)),
    )
    fibre_tension_only = analytical_stress(fibre_extension, tension_only=True)
    assert_close(
        "J.stretched family 1",
        fibre_tension_only["sigma_fibre1"],
        fibre_result["sigma_fibre1"],
    )
    tests_run += 1
    print("PASS J: tension-only option")

    general_deformation = np.array(
        [[1.08, 0.11, -0.03], [0.02, 0.97, 0.07], [0.01, -0.04, 1.03]]
    )
    for split in (False, True):
        analytical = analytical_stress(general_deformation, split=split)[
            "sigma_passive"
        ]
        numerical = numerical_cauchy_from_energy(general_deformation, split=split)
        assert_close(
            f"K.energy-stress split={split}",
            analytical,
            numerical,
            atol=2.0e-6,
            rtol=ENERGY_STRESS_RTOL,
        )
    tests_run += 1
    print("PASS K: energy-stress consistency for split and unsplit forms")

    cell_result = analytical_stress(general_deformation, split=True)
    face_result = analytical_stress(general_deformation.copy(), split=True)
    for key in (
        "J",
        "I4",
        "I6",
        "sigma_iso",
        "sigma_fibre1",
        "sigma_fibre2",
        "sigma_passive",
    ):
        assert_close(f"L.{key}", cell_result[key], face_result[key])
    tests_run += 1
    print("PASS L: homogeneous cell-versus-face consistency")

    baseline_stiffness = effective_implicit_stiffness(2010.0, 1.0)
    doubled_stiffness = effective_implicit_stiffness(4020.0, 1.0)
    assert_close("M.stiffness ratio", doubled_stiffness / baseline_stiffness, 2.0)
    assert_close(
        "M.constitutive independence",
        analytical_stress(simple_shear)["sigma_passive"],
        shear_result["sigma_passive"],
    )
    tests_run += 1
    print("PASS M: explicit implicitShearModulus sensitivity")

    pressure = 37.0
    bulk_modulus = 1000.0
    doubled_bulk_modulus = 2.0 * bulk_modulus
    assert_close(
        "N.pressure residual ratio",
        mixed_pressure_compressibility_term(pressure, doubled_bulk_modulus)
        / mixed_pressure_compressibility_term(pressure, bulk_modulus),
        0.5,
    )
    assert_close(
        "N.constitutive bulk independence",
        analytical_stress(simple_shear)["sigma_passive"],
        shear_result["sigma_passive"],
    )
    for invalid_bulk_modulus in (0.0, -1.0, math.nan, math.inf):
        try:
            mixed_pressure_compressibility_term(
                pressure,
                invalid_bulk_modulus,
            )
        except ValueError:
            pass
        else:
            raise AssertionError(
                f"N: accepted invalid bulk modulus {invalid_bulk_modulus}"
            )
    tests_run += 1
    print("PASS N: finite mixed bulk-modulus separation")

    print(
        f"PASS: {tests_run} material-level checks; "
        f"atol={ATOL:.1e}, rtol={RTOL:.1e}, "
        f"energy-stress rtol={ENERGY_STRESS_RTOL:.1e}"
    )


if __name__ == "__main__":
    try:
        run_tests()
    except Exception as error:  # noqa: BLE001 - test runner must report any failure
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
