#!/usr/bin/env python3
"""Focused algebra and source-contract tests for the Gultekin FV solid."""

from __future__ import annotations

import math
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = ROOT / "src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispSolid"
NEW_DIR = ROOT / "src/solids4FoamModels/solidModels/nonLinGeomTotalLagTotalDispGultekinSolid"
LAW_DIR = ROOT / (
    "src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/"
    "nonLinearGeometryLaws/GultekinTwoFibreElastic"
)

BASE_H = (BASE_DIR / "nonLinGeomTotalLagTotalDispSolid.H").read_text()
BASE_C = (BASE_DIR / "nonLinGeomTotalLagTotalDispSolid.C").read_text()
NEW_H = (NEW_DIR / "nonLinGeomTotalLagTotalDispGultekinSolid.H").read_text()
NEW_C = (NEW_DIR / "nonLinGeomTotalLagTotalDispGultekinSolid.C").read_text()
LAW_C = (LAW_DIR / "GultekinTwoFibreElastic.C").read_text()
LAW_H = (LAW_DIR / "GultekinTwoFibreElastic.H").read_text()
SOLID_MODEL_C = (
    ROOT / "src/solids4FoamModels/solidModels/solidModel/solidModel.C"
).read_text()


def g_gultekin(j: float) -> float:
    return (j - 1.0) / j


def dg_gultekin(j: float) -> float:
    return 1.0 / (j * j)


def diag_add(*tensors: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(sum(values) for values in zip(*tensors))  # type: ignore[return-value]


def diag_dev(tensor: tuple[float, float, float]) -> tuple[float, float, float]:
    mean = sum(tensor) / 3.0
    return tuple(value - mean for value in tensor)  # type: ignore[return-value]


class GultekinMixedSolidTests(unittest.TestCase):
    K = 5.0e6
    J_VALUES = (1.0, 1.1, 1.2, 1.4, 1.6)

    def test_A_scalar_volumetric_law(self) -> None:
        for j in self.J_VALUES:
            p = -self.K * g_gultekin(j)
            self.assertAlmostEqual(-p, self.K * (j - 1.0) / j, places=9)

    def test_B_pressure_residual_closure(self) -> None:
        for j in self.J_VALUES:
            p = -self.K * g_gultekin(j)
            residual = -p / self.K - g_gultekin(j)
            self.assertLessEqual(abs(residual), 4.0 * math.ulp(1.0))

    def test_C_constraint_derivative(self) -> None:
        for j in self.J_VALUES:
            h = 1.0e-6 * j
            centred = (g_gultekin(j + h) - g_gultekin(j - h)) / (2.0 * h)
            self.assertAlmostEqual(centred, dg_gultekin(j), delta=2.0e-10)

    def test_D_full_fibre_stress_retention(self) -> None:
        sigma_iso = (-1.0, 0.5, 0.5)
        sigma_fibre = (0.0, 6.0, 3.0)
        p = -2.0
        pressure = (-p, -p, -p)

        sigma_new = diag_add(sigma_iso, sigma_fibre, pressure)
        sigma_wrong = diag_add(diag_dev(diag_add(sigma_iso, sigma_fibre)), pressure)

        self.assertEqual(sigma_new, (1.0, 8.5, 5.5))
        self.assertEqual(sigma_wrong, (-2.0, 5.5, 2.5))
        self.assertNotEqual(sigma_new, sigma_wrong)

    def test_E_existing_model_regression_contract(self) -> None:
        self.assertRegex(
            BASE_H,
            r"retainFullPassiveStressInMixedSplit\(\) const\s*\{\s*return false;",
        )
        self.assertIn("return 0.5*(pow(J, 2.0) - 1.0)/J;", BASE_C)
        self.assertIn("return 0.5*(J - 1.0/J);", BASE_C)
        self.assertIn("sigma() = dev(sigma());", BASE_C)
        self.assertIn('return "0.5*(J^2 - 1)/J";', BASE_C)

    def test_F_pressure_counted_once_in_new_path(self) -> None:
        self.assertRegex(
            NEW_H,
            r"retainFullPassiveStressInMixedSplit\(\) const\s*\{\s*return true;",
        )
        isolated_branch = re.search(
            r"if \(retainFullPassiveStressInMixedSplit\(\)\)\s*\{(?P<body>.*?)\n\s*\}",
            BASE_C,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(isolated_branch)
        body = isolated_branch.group("body")  # type: ignore[union-attr]
        self.assertEqual(body.count("sigma() = sigma() - p()*I;"), 1)
        self.assertIn("return;", body)
        self.assertNotIn("sigma() = sigma() - p()*I;", NEW_C)
        self.assertIn("sigma[cellI] = result.sigmaPassive;", LAW_C)

    def test_G_written_stress_is_spatial_cauchy_in_pa(self) -> None:
        self.assertIn('"sigma"', SOLID_MODEL_C)
        self.assertIn("IOobject::AUTO_WRITE", SOLID_MODEL_C)
        self.assertIn("dimForce/dimArea", SOLID_MODEL_C)
        self.assertIn("passive Cauchy stress", LAW_H)
        self.assertIn("sigma = sigmaPassive - p I", NEW_C)

    def test_material_unsplit_contract_is_unchanged(self) -> None:
        self.assertIn("result.sigmaIso = (mu_.value()/result.J)*dev(bBar);", LAW_C)
        self.assertIn(": coefficient*currentFibreDyad;", LAW_C)
        self.assertIn(
            "result.sigmaIso + result.sigmaFibre1 + result.sigmaFibre2;",
            LAW_C,
        )
        self.assertIn("sigmaPatch[faceI] = result.sigmaPassive;", LAW_C)

    def test_runtime_registration_and_constraint(self) -> None:
        self.assertIn(
            'TypeName("nonLinGeomTotalLagTotalDispGultekinSolid")', NEW_H
        )
        self.assertIn("addToRunTimeSelectionTable", NEW_C)
        self.assertIn("return (J - 1.0)/J;", NEW_C)
        self.assertIn("return 1.0/sqr(J);", NEW_C)
        self.assertIn("- tConstraint()", BASE_C)
        self.assertIn("anisotropicSplit false;", NEW_C)
        self.assertIn("fibresTensionOnly false;", NEW_C)
        self.assertIn("useSecondFibreFamily true;", NEW_C)

    def test_no_finite_element_assembly_added(self) -> None:
        forbidden = (
            "elementStiffness",
            "shapeFunction",
            "HuWashizu",
            "finiteElementQuadrature",
            "Q1Element",
            "P0Element",
        )
        for token in forbidden:
            self.assertNotIn(token, NEW_C)


if __name__ == "__main__":
    unittest.main(verbosity=2)
