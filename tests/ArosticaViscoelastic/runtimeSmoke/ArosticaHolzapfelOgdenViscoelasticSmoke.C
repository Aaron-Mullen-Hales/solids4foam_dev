/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalLaw.H"
#include "solidModel.H"

#include <cmath>

using namespace Foam;


namespace
{

class stressTestSolidModel
:
    public solidModel
{
public:

    stressTestSolidModel(Time& runTime)
    :
        solidModel("arosticaViscoSmoke", runTime)
    {}

    virtual nonLinearGeometry::nonLinearType nonLinGeom() const
    {
        return nonLinearGeometry::TOTAL_LAGRANGIAN;
    }

    virtual volVectorField& solutionD()
    {
        return D();
    }

    virtual bool evolve()
    {
        return true;
    }

    virtual tmp<vectorField> tractionBoundarySnGrad
    (
        const vectorField&,
        const scalarField&,
        const fvPatch& patch
    ) const
    {
        return tmp<vectorField>
        (
            new vectorField(patch.size(), vector::zero)
        );
    }
};


scalar maximumDifference
(
    const volSymmTensorField& a,
    const volSymmTensorField& b
)
{
    scalar result = 0.0;

    forAll(a, cellI)
    {
        result = max(result, mag(a[cellI] - b[cellI]));
    }

    forAll(a.boundaryField(), patchI)
    {
        const symmTensorField& aPatch = a.boundaryField()[patchI];
        const symmTensorField& bPatch = b.boundaryField()[patchI];

        forAll(aPatch, faceI)
        {
            result = max(result, mag(aPatch[faceI] - bPatch[faceI]));
        }
    }

    return result;
}


scalar maximumDifference
(
    const surfaceSymmTensorField& a,
    const surfaceSymmTensorField& b
)
{
    scalar result = 0.0;

    forAll(a, faceI)
    {
        result = max(result, mag(a[faceI] - b[faceI]));
    }

    forAll(a.boundaryField(), patchI)
    {
        const symmTensorField& aPatch = a.boundaryField()[patchI];
        const symmTensorField& bPatch = b.boundaryField()[patchI];

        forAll(aPatch, faceI)
        {
            result = max(result, mag(aPatch[faceI] - bPatch[faceI]));
        }
    }

    return result;
}


scalar maximumCellFaceDifference
(
    const volSymmTensorField& cellStress,
    const surfaceSymmTensorField& faceStress
)
{
    scalar result = 0.0;

    forAll(faceStress, faceI)
    {
        result = max(result, mag(faceStress[faceI] - cellStress[0]));
    }

    forAll(faceStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = faceStress.boundaryField()[patchI];
        forAll(patch, faceI)
        {
            result = max(result, mag(patch[faceI] - cellStress[0]));
        }
    }

    return result;
}


scalar maximumCellFaceBoundaryDifference
(
    const volSymmTensorField& cellStress,
    const surfaceSymmTensorField& faceStress
)
{
    scalar result = 0.0;

    forAll(faceStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = faceStress.boundaryField()[patchI];
        forAll(patch, faceI)
        {
            result = max(result, mag(patch[faceI] - cellStress[0]));
        }
    }

    return result;
}

}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    stressTestSolidModel solid(runTime);
    dynamicFvMesh& mesh = solid.mesh();

    IOdictionary lawDictionary
    (
        IOobject
        (
            "lawProperties.visco",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<mechanicalLaw> law = mechanicalLaw::NewNonLinGeomMechLaw
    (
        "viscoSmoke",
        mesh,
        lawDictionary,
        nonLinearGeometry::TOTAL_LAGRANGIAN
    );

    if (law->type() != "ArosticaHolzapfelOgdenViscoelastic")
    {
        FatalErrorInFunction
            << "Unexpected runtime type " << law->type()
            << abort(FatalError);
    }

    volTensorField& gradD = const_cast<volTensorField&>
    (
        mesh.lookupObject<volTensorField>("grad(D)")
    );
    gradD = tensor
    (
        0.0, 0.15, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );

    surfaceTensorField gradDf
    (
        IOobject
        (
            "grad(D)f", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "gradDf", dimless,
            tensor
            (
                0.0, 0.15, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        )
    );

    auto setTrialGradient =
        [&]
        (
            const scalar xx,
            const scalar xy,
            const scalar xz,
            const scalar yx,
            const scalar yy,
            const scalar yz,
            const scalar zx,
            const scalar zy,
            const scalar zz
        )
        {
            const tensor value(xx, xy, xz, yx, yy, yz, zx, zy, zz);
            gradD = value;
            gradDf = value;
        };

    volSymmTensorField cellStress
    (
        IOobject
        (
            "cellStress", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );
    surfaceSymmTensorField faceStress
    (
        IOobject
        (
            "faceStress", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    law->correct(cellStress);
    law->correct(faceStress);
    const scalar initialCellFaceDifference =
        maximumCellFaceDifference(cellStress, faceStress);

    const volSymmTensorField stressA(cellStress);
    const surfaceSymmTensorField stressAf(faceStress);

    const volSymmTensorField& E =
        mesh.lookupObject<volSymmTensorField>("ArosticaE_viscoSmoke");
    const volSymmTensorField& EOld =
        mesh.lookupObject<volSymmTensorField>("ArosticaEOld_viscoSmoke");
    const volSymmTensorField& EOldOld =
        mesh.lookupObject<volSymmTensorField>
        (
            "ArosticaEOldOld_viscoSmoke"
        );
    const surfaceSymmTensorField& Ef =
        mesh.lookupObject<surfaceSymmTensorField>("ArosticaEf_viscoSmoke");
    const surfaceSymmTensorField& EfOld =
        mesh.lookupObject<surfaceSymmTensorField>
        (
            "ArosticaEfOld_viscoSmoke"
        );
    const surfaceSymmTensorField& EfOldOld =
        mesh.lookupObject<surfaceSymmTensorField>
        (
            "ArosticaEfOldOld_viscoSmoke"
        );
    const volSymmTensorField& Edot =
        mesh.lookupObject<volSymmTensorField>("ArosticaEdot_viscoSmoke");
    const surfaceSymmTensorField& Edotf =
        mesh.lookupObject<surfaceSymmTensorField>
        (
            "ArosticaEdotf_viscoSmoke"
        );

    scalar firstStepEdotError = 0.0;
    forAll(Edot, cellI)
    {
        firstStepEdotError = max
        (
            firstStepEdotError,
            mag(Edot[cellI] - E[cellI]/runTime.deltaTValue())
        );
    }
    forAll(Edot.boundaryField(), patchI)
    {
        const symmTensorField& edotPatch = Edot.boundaryField()[patchI];
        const symmTensorField& ePatch = E.boundaryField()[patchI];
        forAll(edotPatch, faceI)
        {
            firstStepEdotError = max
            (
                firstStepEdotError,
                mag(edotPatch[faceI] - ePatch[faceI]/runTime.deltaTValue())
            );
        }
    }
    forAll(Edotf, faceI)
    {
        firstStepEdotError = max
        (
            firstStepEdotError,
            mag(Edotf[faceI] - Ef[faceI]/runTime.deltaTValue())
        );
    }
    forAll(Edotf.boundaryField(), patchI)
    {
        const symmTensorField& edotPatch = Edotf.boundaryField()[patchI];
        const symmTensorField& ePatch = Ef.boundaryField()[patchI];
        forAll(edotPatch, faceI)
        {
            firstStepEdotError = max
            (
                firstStepEdotError,
                mag(edotPatch[faceI] - ePatch[faceI]/runTime.deltaTValue())
            );
        }
    }

    if (firstStepEdotError > 1e-10)
    {
        FatalErrorInFunction
            << "First-step Euler fallback is incorrect; maximum Edot error = "
            << firstStepEdotError << abort(FatalError);
    }

    // updateTotalFields() is the sole accepted-history commit.  Capture its
    // input and verify both the shift and the same-time double-call guard.
    const volSymmTensorField EOldBeforeCommit(EOld);
    const volSymmTensorField EOldOldBeforeCommit(EOldOld);
    const surfaceSymmTensorField EfOldBeforeCommit(EfOld);
    const surfaceSymmTensorField EfOldOldBeforeCommit(EfOldOld);
    const volSymmTensorField EAccepted(E);
    const surfaceSymmTensorField EfAccepted(Ef);

    setTrialGradient
    (
        0.0, 0.11, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    law->correct(cellStress);
    law->correct(faceStress);
    const scalar trialBCellFaceDifference =
        maximumCellFaceDifference(cellStress, faceStress);
    const volSymmTensorField stressB(cellStress);
    const surfaceSymmTensorField stressBf(faceStress);

    const scalar distinctTrialDifference =
        max
        (
            maximumDifference(stressA, stressB),
            maximumDifference(stressAf, stressBf)
        );

    if (distinctTrialDifference <= 1e-10)
    {
        FatalErrorInFunction
            << "Distinct trial B did not change stress"
            << abort(FatalError);
    }

    setTrialGradient
    (
        0.0, 0.15, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    law->correct(cellStress);
    law->correct(faceStress);
    const scalar repeatedTrialDifference =
        max
        (
            maximumDifference(cellStress, stressA),
            maximumDifference(faceStress, stressAf)
        );

    if (repeatedTrialDifference > 1e-10)
    {
        FatalErrorInFunction
            << "Repeated trial A changed stress; maximum difference = "
            << repeatedTrialDifference << abort(FatalError);
    }

    // Commit A, then verify that a rejected C trial does not alter the next
    // evaluation of B. updateTotalFields() is the only history commit.
    law->updateTotalFields();

    const scalar firstCellHistoryShift =
        max
        (
            maximumDifference(EOldOld, EOldBeforeCommit),
            maximumDifference(EOld, EAccepted)
        );
    const scalar firstFaceHistoryShift =
        max
        (
            maximumDifference(EfOldOld, EfOldOldBeforeCommit),
            maximumDifference(EfOld, EfAccepted)
        );

    if (firstCellHistoryShift > 1e-10 || firstFaceHistoryShift > 1e-10)
    {
        FatalErrorInFunction
            << "Accepted history shift is incorrect; cell difference = "
            << firstCellHistoryShift << ", face difference = "
            << firstFaceHistoryShift << abort(FatalError);
    }

    const volSymmTensorField EOldAfterCommit(EOld);
    const volSymmTensorField EOldOldAfterCommit(EOldOld);
    const surfaceSymmTensorField EfOldAfterCommit(EfOld);
    const surfaceSymmTensorField EfOldOldAfterCommit(EfOldOld);

    // A second finalisation call at this same physical time must be a no-op.
    law->updateTotalFields();
    const scalar doubleCellHistoryShift =
        max
        (
            maximumDifference(EOld, EOldAfterCommit),
            maximumDifference(EOldOld, EOldOldAfterCommit)
        );
    const scalar doubleFaceHistoryShift =
        max
        (
            maximumDifference(EfOld, EfOldAfterCommit),
            maximumDifference(EfOldOld, EfOldOldAfterCommit)
        );

    if (doubleCellHistoryShift > 1e-10 || doubleFaceHistoryShift > 1e-10)
    {
        FatalErrorInFunction
            << "A second same-time history update changed accepted state; "
            << "cell difference = " << doubleCellHistoryShift
            << ", face difference = " << doubleFaceHistoryShift
            << abort(FatalError);
    }

    scalar maximumTrialHistoryChange = 0.0;
    setTrialGradient
    (
        0.0, 0.11, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    law->correct(cellStress);
    law->correct(faceStress);
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EOld, EOldAfterCommit),
            maximumDifference(EOldOld, EOldOldAfterCommit)
        )
    );
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EfOld, EfOldAfterCommit),
            maximumDifference(EfOldOld, EfOldOldAfterCommit)
        )
    );
    const volSymmTensorField acceptedStressB(cellStress);
    const surfaceSymmTensorField acceptedStressBf(faceStress);

    setTrialGradient
    (
        0.0, 0.23, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    law->correct(cellStress);
    law->correct(faceStress);
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EOld, EOldAfterCommit),
            maximumDifference(EOldOld, EOldOldAfterCommit)
        )
    );
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EfOld, EfOldAfterCommit),
            maximumDifference(EfOldOld, EfOldOldAfterCommit)
        )
    );
    setTrialGradient
    (
        0.0, 0.11, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    law->correct(cellStress);
    law->correct(faceStress);
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EOld, EOldAfterCommit),
            maximumDifference(EOldOld, EOldOldAfterCommit)
        )
    );
    maximumTrialHistoryChange = max
    (
        maximumTrialHistoryChange,
        max
        (
            maximumDifference(EfOld, EfOldAfterCommit),
            maximumDifference(EfOldOld, EfOldOldAfterCommit)
        )
    );
    const scalar rejectedTrialDifference =
        max
        (
            maximumDifference(cellStress, acceptedStressB),
            maximumDifference(faceStress, acceptedStressBf)
        );

    if (rejectedTrialDifference > 1e-10)
    {
        FatalErrorInFunction
            << "Rejected trial changed later stress; maximum difference = "
            << rejectedTrialDifference << abort(FatalError);
    }

    const scalar cellFaceDifference =
        maximumCellFaceDifference(cellStress, faceStress);
    const scalar cellFaceBoundaryDifference =
        maximumCellFaceBoundaryDifference(cellStress, faceStress);

    // Advance a controlled, non-zero history through unequal time steps. The
    // expected coefficients are evaluated independently here from the
    // OpenFOAM time values; the law's Edot field is then compared for cells
    // and direct faces before the single accepted-history commit.
    scalar variableDtEdotError = 0.0;
    scalar variableDtCellFaceError = 0.0;
    const scalarField variableDts
    (
        List<scalar>({0.010, 0.015, 0.008, 0.020, 0.011})
    );
    const scalarField variableGammas
    (
        List<scalar>({0.010, 0.025, 0.018, 0.035, 0.012})
    );

    forAll(variableDts, stepI)
    {
        runTime.setDeltaT(variableDts[stepI]);
        runTime++;

        setTrialGradient
        (
            0.0, variableGammas[stepI], 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        );
        law->correct(cellStress);
        law->correct(faceStress);

        const scalar deltaT = runTime.deltaTValue();
        const scalar deltaT0 = runTime.deltaT0Value();
        const scalar coefft = 1.0 + deltaT/(deltaT + deltaT0);
        const scalar coefft00 =
            deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        const scalar coeffCurrent = coefft/deltaT;
        const scalar coeffOld = -(coefft + coefft00)/deltaT;
        const scalar coeffOldOld = coefft00/deltaT;

        const tensor trialGradient
        (
            0.0, variableGammas[stepI], 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        );
        const tensor trialF(tensor::I + trialGradient.T());
        const symmTensor trialE
        (
            0.5*(symm(trialF.T() & trialF) - symmTensor::I)
        );

        forAll(Edot, cellI)
        {
            const symmTensor expected
            (
                coeffCurrent*trialE
              + coeffOld*EOld[cellI]
              + coeffOldOld*EOldOld[cellI]
            );
            variableDtEdotError = max
            (
                variableDtEdotError,
                mag(Edot[cellI] - expected)
            );
        }
        forAll(Edotf, faceI)
        {
            const symmTensor expected
            (
                coeffCurrent*trialE
              + coeffOld*EfOld[faceI]
              + coeffOldOld*EfOldOld[faceI]
            );
            variableDtEdotError = max
            (
                variableDtEdotError,
                mag(Edotf[faceI] - expected)
            );
            variableDtCellFaceError = max
            (
                variableDtCellFaceError,
                mag(Edotf[faceI] - Edot[0])
            );
        }

        law->updateTotalFields();
    }

    if (variableDtEdotError > 1e-10)
    {
        FatalErrorInFunction
            << "Variable-deltaT backward Edot mismatch; maximum error = "
            << variableDtEdotError << abort(FatalError);
    }

    IOdictionary zeroDictionary
    (
        IOobject
        (
            "lawProperties.viscoZero", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    autoPtr<mechanicalLaw> zeroLaw = mechanicalLaw::NewNonLinGeomMechLaw
    (
        "viscoZero",
        mesh,
        zeroDictionary,
        nonLinearGeometry::TOTAL_LAGRANGIAN
    );

    IOdictionary elasticDictionary
    (
        IOobject
        (
            "lawProperties.elastic", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    autoPtr<mechanicalLaw> elasticLaw = mechanicalLaw::NewNonLinGeomMechLaw
    (
        "elasticReference",
        mesh,
        elasticDictionary,
        nonLinearGeometry::TOTAL_LAGRANGIAN
    );

    volSymmTensorField zeroStress(cellStress);
    volSymmTensorField elasticStress(cellStress);
    surfaceSymmTensorField zeroStressf(faceStress);
    surfaceSymmTensorField elasticStressf(faceStress);
    zeroLaw->correct(zeroStress);
    zeroLaw->correct(zeroStressf);
    elasticLaw->correct(elasticStress);
    elasticLaw->correct(elasticStressf);

    const scalar zeroEtaDifference =
        max
        (
            maximumDifference(zeroStress, elasticStress),
            maximumDifference(zeroStressf, elasticStressf)
        );

    if (zeroEtaDifference > 1e-10)
    {
        FatalErrorInFunction
            << "eta=0 does not reproduce elasticity; maximum difference = "
            << zeroEtaDifference << abort(FatalError);
    }

    Info<< "PASS: Aróstica viscoelastic runtime smoke" << nl
        << "    runtime type = " << law->type() << nl
        << "    distinct-trial maximum difference = "
        << distinctTrialDifference << nl
        << "    initial cell/face maximum difference = "
        << initialCellFaceDifference << nl
        << "    trial-B cell/face maximum difference = "
        << trialBCellFaceDifference << nl
        << "    repeated-trial maximum difference = "
        << repeatedTrialDifference << nl
        << "    rejected-trial maximum difference = "
        << rejectedTrialDifference << nl
        << "    first cell history-shift difference = "
        << firstCellHistoryShift << nl
        << "    first face history-shift difference = "
        << firstFaceHistoryShift << nl
        << "    same-time double-update cell difference = "
        << doubleCellHistoryShift << nl
        << "    same-time double-update face difference = "
        << doubleFaceHistoryShift << nl
        << "    maximum trial history change = "
        << maximumTrialHistoryChange << nl
        << "    first-step Edot maximum error = "
        << firstStepEdotError << nl
        << "    variable-deltaT Edot maximum error = "
        << variableDtEdotError << nl
        << "    variable-deltaT cell/face Edot difference = "
        << variableDtCellFaceError << nl
        << "    cell/face maximum difference = " << cellFaceDifference << nl
        << "    cell/face boundary maximum difference = "
        << cellFaceBoundaryDifference << nl
        << "    zero-eta maximum difference = " << zeroEtaDifference
        << endl;

    return 0;
}

// ************************************************************************* //
