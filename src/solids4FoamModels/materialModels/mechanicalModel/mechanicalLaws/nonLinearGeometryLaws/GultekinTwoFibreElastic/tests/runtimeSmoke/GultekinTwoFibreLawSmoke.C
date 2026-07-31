/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    GultekinTwoFibreLawSmoke

Description
    Instantiates a nonlinear mechanical law from a one-cell test dictionary.
    For GultekinTwoFibreElastic, calls both correct() overloads, compares the
    returned passive Cauchy stress with an independent analytical value, and
    optionally checks sensitivity to a second numerical stiffness dictionary.

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
        solidModel("stressTestSolid", runTime)
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

struct StressTestResult
{
    symmTensor cellStress;
    symmTensor faceStress;
    scalar bulkModulus;
    scalar effectiveImplicitStiffness;
};


symmTensor fibreStress
(
    const tensor& F,
    const symmTensor& C,
    const vector& fibre,
    const scalar J,
    const scalar k1,
    const scalar k2
)
{
    const scalar invariant = fibre & (C & fibre);
    const scalar strain = invariant - 1.0;
    const vector currentFibre(F & fibre);

    return
        (2.0*k1*strain*std::exp(k2*sqr(strain))/J)
       *symm(currentFibre*currentFibre);
}


symmTensor expectedPassiveStress
(
    const dictionary& dict,
    const tensor& F
)
{
    if (dict.lookupOrDefault<Switch>("anisotropicSplit", false))
    {
        FatalErrorInFunction
            << "The compiled one-cell fixture expects anisotropicSplit false"
            << abort(FatalError);
    }

    if (dict.lookupOrDefault<Switch>("fibresTensionOnly", false))
    {
        FatalErrorInFunction
            << "The compiled one-cell fixture expects fibresTensionOnly false"
            << abort(FatalError);
    }

#ifdef OPENFOAM_COM
    const dimensionedScalar mu("mu", dict);
    const dimensionedScalar k1("k1", dict);
    const dimensionedScalar k2("k2", dict);
#else
    const dimensionedScalar mu(dict.lookup("mu"));
    const dimensionedScalar k1(dict.lookup("k1"));
    const dimensionedScalar k2(dict.lookup("k2"));
#endif
    const scalar J = det(F);
    const scalar Jm23 = std::pow(J, -2.0/3.0);
    const symmTensor C(symm(F.T() & F));
    const symmTensor bBar(Jm23*symm(F & F.T()));
    const scalar inverseSqrtTwo = 1.0/std::sqrt(2.0);
    const vector fibre1(inverseSqrtTwo, inverseSqrtTwo, 0.0);
    const vector fibre2(inverseSqrtTwo, -inverseSqrtTwo, 0.0);

    symmTensor stress((mu.value()/J)*dev(bBar));

    stress += fibreStress
    (
        F, C, fibre1, J, k1.value(), k2.value()
    );

    if (dict.lookupOrDefault<Switch>("useSecondFibreFamily", true))
    {
        stress += fibreStress
        (
            F, C, fibre2, J, k1.value(), k2.value()
        );
    }

    return stress;
}


void checkStress
(
    const symmTensor& actual,
    const symmTensor& expected,
    const string& location,
    const scalar tolerance,
    scalar& maximumError
)
{
    const scalar error = mag(actual - expected);
    maximumError = max(maximumError, error);

    if (!std::isfinite(error) || error > tolerance)
    {
        FatalErrorInFunction
            << "Stress mismatch at " << location << nl
            << "    actual = " << actual << nl
            << "    expected = " << expected << nl
            << "    error = " << error << nl
            << "    tolerance = " << tolerance
            << abort(FatalError);
    }
}


void checkUniformScalarField
(
    const volScalarField& field,
    const scalar expected,
    const word& description,
    const scalar tolerance
)
{
    scalar maximumError = 0.0;
    label valueCount = 0;

    forAll(field, cellI)
    {
        maximumError = max(maximumError, mag(field[cellI] - expected));
        ++valueCount;
    }

    forAll(field.boundaryField(), patchI)
    {
        const scalarField& patch = field.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            maximumError = max
            (
                maximumError,
                mag(patch[faceI] - expected)
            );
            ++valueCount;
        }
    }

    if
    (
        !valueCount
     || !std::isfinite(maximumError)
     || maximumError > tolerance
    )
    {
        FatalErrorInFunction
            << description << " field mismatch" << nl
            << "    expected = " << expected << nl
            << "    maximum error = " << maximumError << nl
            << "    tolerance = " << tolerance << nl
            << "    checked values = " << valueCount
            << abort(FatalError);
    }

    Info<< "PASS: " << description << " API" << nl
        << "    expected value = " << expected << nl
        << "    checked values = " << valueCount << nl
        << "    maximum error = " << maximumError << endl;
}


StressTestResult runStressTest
(
    mechanicalLaw& law,
    const dictionary& dict,
    const tensor& F,
    volSymmTensorField& cellStress,
    surfaceSymmTensorField& faceStress
)
{
    cellStress = dimensionedSymmTensor
    (
        "zero",
        dimPressure,
        symmTensor::zero
    );
    faceStress = dimensionedSymmTensor
    (
        "zero",
        dimPressure,
        symmTensor::zero
    );

    law.correct(cellStress);
    law.correct(faceStress);

    const symmTensor expected(expectedPassiveStress(dict, F));
    const scalar tolerance = 1e-9*max(scalar(1.0), mag(expected));
    scalar maximumCellError = 0.0;
    scalar maximumFaceError = 0.0;
    label cellValueCount = 0;
    label faceValueCount = 0;
    symmTensor representativeFaceStress(symmTensor::zero);

    forAll(cellStress, cellI)
    {
        checkStress
        (
            cellStress[cellI],
            expected,
            "cell " + Foam::name(cellI),
            tolerance,
            maximumCellError
        );
        ++cellValueCount;
    }

    forAll(cellStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = cellStress.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            checkStress
            (
                patch[faceI],
                expected,
                "cell-field boundary face " + Foam::name(faceI),
                tolerance,
                maximumCellError
            );
            ++cellValueCount;
        }
    }

    forAll(faceStress, faceI)
    {
        if (!faceValueCount)
        {
            representativeFaceStress = faceStress[faceI];
        }

        checkStress
        (
            faceStress[faceI],
            expected,
            "internal surface face " + Foam::name(faceI),
            tolerance,
            maximumFaceError
        );
        ++faceValueCount;
    }

    forAll(faceStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = faceStress.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            if (!faceValueCount)
            {
                representativeFaceStress = patch[faceI];
            }

            checkStress
            (
                patch[faceI],
                expected,
                "surface-field boundary face " + Foam::name(faceI),
                tolerance,
                maximumFaceError
            );
            ++faceValueCount;
        }
    }

    if (!cellValueCount || !faceValueCount)
    {
        FatalErrorInFunction
            << "The one-cell test did not exercise both field locations"
            << abort(FatalError);
    }

    if (!dict.found("implicitShearModulus"))
    {
        FatalIOErrorInFunction(dict)
            << "The compiled test requires an explicit dimensioned "
            << "implicitShearModulus entry"
            << exit(FatalIOError);
    }

    if (!dict.found("bulkModulus"))
    {
        FatalIOErrorInFunction(dict)
            << "The compiled test requires an explicit dimensioned "
            << "bulkModulus entry"
            << exit(FatalIOError);
    }

#ifdef OPENFOAM_COM
    const dimensionedScalar specifiedStiffness
    (
        "implicitShearModulus",
        dict
    );
    const dimensionedScalar specifiedBulkModulus
    (
        "bulkModulus",
        dict
    );
#else
    const dimensionedScalar specifiedStiffness
    (
        dict.lookup("implicitShearModulus")
    );
    const dimensionedScalar specifiedBulkModulus
    (
        dict.lookup("bulkModulus")
    );
#endif
    const scalar impKcoeff =
        dict.lookupOrDefault<scalar>("impKcoeff", 1.0);
    const scalar expectedEffectiveStiffness =
        impKcoeff*specifiedStiffness.value();
    const tmp<volScalarField> tShearModulus = law.shearModulus();
    const tmp<volScalarField> tImpK = law.impK();
    const tmp<volScalarField> tBulkModulus = law.bulkModulus();
    const volScalarField& shearModulus = tShearModulus();
    const volScalarField& impK = tImpK();
    const volScalarField& bulkModulus = tBulkModulus();
    const scalar stiffnessTolerance =
        1e-12*max(scalar(1.0), mag(expectedEffectiveStiffness));
    const scalar bulkTolerance =
        1e-12*max(scalar(1.0), mag(specifiedBulkModulus.value()));

    if (bulkModulus.dimensions() != dimPressure)
    {
        FatalErrorInFunction
            << "bulkModulus() returned dimensions "
            << bulkModulus.dimensions() << "; expected " << dimPressure
            << abort(FatalError);
    }

    checkUniformScalarField
    (
        shearModulus,
        expectedEffectiveStiffness,
        "shearModulus",
        stiffnessTolerance
    );
    checkUniformScalarField
    (
        impK,
        expectedEffectiveStiffness,
        "impK",
        stiffnessTolerance
    );
    checkUniformScalarField
    (
        bulkModulus,
        specifiedBulkModulus.value(),
        "bulkModulus",
        bulkTolerance
    );

    Info<< "PASS: both correct() overloads returned passive Cauchy stress" << nl
        << "    evaluated cell values = " << cellValueCount << nl
        << "    evaluated face values = " << faceValueCount << nl
        << "    maximum cell error = " << maximumCellError << nl
        << "    maximum face error = " << maximumFaceError << nl
        << "    nonzero registered p and pf were ignored" << nl
        << "    specified implicitShearModulus = "
        << specifiedStiffness << nl
        << "    impKcoeff = " << impKcoeff << nl
        << "    effective numerical stiffness = "
        << expectedEffectiveStiffness << nl
        << "    specified bulkModulus = "
        << specifiedBulkModulus << endl;

    StressTestResult result;
    result.cellStress = cellStress[0];
    result.faceStress = representativeFaceStress;
    result.bulkModulus = bulkModulus[0];
    result.effectiveImplicitStiffness = expectedEffectiveStiffness;
    return result;
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "dict",
        "word",
        "Dictionary name in the constant directory"
    );
    argList::addOption
    (
        "compareDict",
        "word",
        "Second Gultekin dictionary for numerical-stiffness sensitivity"
    );
    argList::addOption
    (
        "compareBulkDict",
        "word",
        "Second Gultekin dictionary for bulk-modulus sensitivity"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    stressTestSolidModel solid(runTime);
    dynamicFvMesh& mesh = solid.mesh();

    word dictionaryName("lawProperties.gtf");
    word comparisonDictionaryName(word::null);
    word bulkComparisonDictionaryName(word::null);

#ifdef OPENFOAM_COM
    args.readIfPresent("dict", dictionaryName);
    args.readIfPresent("compareDict", comparisonDictionaryName);
    args.readIfPresent("compareBulkDict", bulkComparisonDictionaryName);
#else
    args.optionReadIfPresent("dict", dictionaryName);
    args.optionReadIfPresent("compareDict", comparisonDictionaryName);
    args.optionReadIfPresent
    (
        "compareBulkDict",
        bulkComparisonDictionaryName
    );
#endif

    IOdictionary lawDictionary
    (
        IOobject
        (
            dictionaryName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word requestedType(lawDictionary.lookup("type"));
    autoPtr<mechanicalLaw> law = mechanicalLaw::NewNonLinGeomMechLaw
    (
        "runtimeSmoke",
        mesh,
        lawDictionary,
        nonLinearGeometry::TOTAL_LAGRANGIAN
    );

    if (law->type() != requestedType)
    {
        FatalErrorInFunction
            << "Requested " << requestedType << " but constructed "
            << law->type() << abort(FatalError);
    }

    Info<< "PASS: instantiated runtime mechanical law " << law->type()
        << " from " << dictionaryName << endl;

    if (requestedType == "GultekinTwoFibreElastic")
    {
        const scalar shearAmount = 0.2;
        const tensor prescribedGradD
        (
            0.0, 0.0, 0.0,
            shearAmount, 0.0, 0.0,
            0.0, 0.0, 0.0
        );
        const tensor prescribedF
        (
            1.0, shearAmount, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        );

        volTensorField& gradD = const_cast<volTensorField&>
        (
            mesh.lookupObject<volTensorField>("grad(D)")
        );
        gradD = prescribedGradD;
        surfaceTensorField gradDf
        (
            IOobject
            (
                "grad(D)f",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("gradDf", dimless, prescribedGradD)
        );
        volScalarField pressure
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("p", dimPressure, 12345.0)
        );
        surfaceScalarField facePressure
        (
            IOobject
            (
                "pf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("pf", dimPressure, -6789.0)
        );
        volSymmTensorField cellStress
        (
            IOobject
            (
                "cellStressTest",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            (
                "zero",
                dimPressure,
                symmTensor::zero
            )
        );
        surfaceSymmTensorField faceStress
        (
            IOobject
            (
                "faceStressTest",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            (
                "zero",
                dimPressure,
                symmTensor::zero
            )
        );

        const StressTestResult primaryResult = runStressTest
        (
            law(),
            lawDictionary,
            prescribedF,
            cellStress,
            faceStress
        );

        if (comparisonDictionaryName.size())
        {
            law.clear();

            IOdictionary comparisonDictionary
            (
                IOobject
                (
                    comparisonDictionaryName,
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const word comparisonType
            (
                comparisonDictionary.lookup("type")
            );

            if (comparisonType != requestedType)
            {
                FatalIOErrorInFunction(comparisonDictionary)
                    << "Comparison type " << comparisonType
                    << " differs from primary type " << requestedType
                    << exit(FatalIOError);
            }

            law = mechanicalLaw::NewNonLinGeomMechLaw
            (
                "runtimeSmokeSensitivity",
                mesh,
                comparisonDictionary,
                nonLinearGeometry::TOTAL_LAGRANGIAN
            );

            const StressTestResult comparisonResult = runStressTest
            (
                law(),
                comparisonDictionary,
                prescribedF,
                cellStress,
                faceStress
            );
            const scalar cellStressDifference = mag
            (
                comparisonResult.cellStress - primaryResult.cellStress
            );
            const scalar faceStressDifference = mag
            (
                comparisonResult.faceStress - primaryResult.faceStress
            );
            const scalar stressTolerance =
                1e-12*max(scalar(1.0), mag(primaryResult.cellStress));
            const scalar bulkTolerance =
                1e-12*max(scalar(1.0), mag(primaryResult.bulkModulus));

            if
            (
                cellStressDifference > stressTolerance
             || faceStressDifference > stressTolerance
             || comparisonResult.effectiveImplicitStiffness
             == primaryResult.effectiveImplicitStiffness
             || mag
                (
                    comparisonResult.bulkModulus
                  - primaryResult.bulkModulus
                ) > bulkTolerance
            )
            {
                FatalErrorInFunction
                    << "Unexpected implicitShearModulus sensitivity" << nl
                    << "    primary effective stiffness = "
                    << primaryResult.effectiveImplicitStiffness << nl
                    << "    comparison effective stiffness = "
                    << comparisonResult.effectiveImplicitStiffness << nl
                    << "    cell stress difference = "
                    << cellStressDifference << nl
                    << "    face stress difference = "
                    << faceStressDifference << nl
                    << "    primary bulkModulus = "
                    << primaryResult.bulkModulus << nl
                    << "    comparison bulkModulus = "
                    << comparisonResult.bulkModulus << nl
                    << "    stress tolerance = " << stressTolerance
                    << abort(FatalError);
            }

            Info<< "PASS: implicitShearModulus sensitivity" << nl
                << "    primary effective stiffness = "
                << primaryResult.effectiveImplicitStiffness << nl
                << "    comparison effective stiffness = "
                << comparisonResult.effectiveImplicitStiffness << nl
                << "    stiffness ratio = "
                << comparisonResult.effectiveImplicitStiffness
                  /primaryResult.effectiveImplicitStiffness << nl
                << "    bulkModulus difference = "
                << comparisonResult.bulkModulus
                 - primaryResult.bulkModulus << nl
                << "    cell stress difference = "
                << cellStressDifference << nl
                << "    face stress difference = "
                << faceStressDifference << endl;
        }

        if (bulkComparisonDictionaryName.size())
        {
            law.clear();

            IOdictionary bulkComparisonDictionary
            (
                IOobject
                (
                    bulkComparisonDictionaryName,
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const word bulkComparisonType
            (
                bulkComparisonDictionary.lookup("type")
            );

            if (bulkComparisonType != requestedType)
            {
                FatalIOErrorInFunction(bulkComparisonDictionary)
                    << "Comparison type " << bulkComparisonType
                    << " differs from primary type " << requestedType
                    << exit(FatalIOError);
            }

            law = mechanicalLaw::NewNonLinGeomMechLaw
            (
                "runtimeSmokeBulkSensitivity",
                mesh,
                bulkComparisonDictionary,
                nonLinearGeometry::TOTAL_LAGRANGIAN
            );

            const StressTestResult bulkComparisonResult = runStressTest
            (
                law(),
                bulkComparisonDictionary,
                prescribedF,
                cellStress,
                faceStress
            );
            const scalar cellStressDifference = mag
            (
                bulkComparisonResult.cellStress - primaryResult.cellStress
            );
            const scalar faceStressDifference = mag
            (
                bulkComparisonResult.faceStress - primaryResult.faceStress
            );
            const scalar stressTolerance =
                1e-12*max(scalar(1.0), mag(primaryResult.cellStress));
            const scalar stiffnessTolerance =
                1e-12*max
                (
                    scalar(1.0),
                    mag(primaryResult.effectiveImplicitStiffness)
                );
            const scalar bulkRatio =
                bulkComparisonResult.bulkModulus/primaryResult.bulkModulus;

#ifdef OPENFOAM_COM
            const dimensionedScalar primaryBulk("bulkModulus", lawDictionary);
            const dimensionedScalar comparisonBulk
            (
                "bulkModulus",
                bulkComparisonDictionary
            );
#else
            const dimensionedScalar primaryBulk
            (
                lawDictionary.lookup("bulkModulus")
            );
            const dimensionedScalar comparisonBulk
            (
                bulkComparisonDictionary.lookup("bulkModulus")
            );
#endif
            const scalar expectedBulkRatio =
                comparisonBulk.value()/primaryBulk.value();
            const scalar ratioTolerance =
                1e-12*max(scalar(1.0), mag(expectedBulkRatio));

            if
            (
                cellStressDifference > stressTolerance
             || faceStressDifference > stressTolerance
             || mag
                (
                    bulkComparisonResult.effectiveImplicitStiffness
                  - primaryResult.effectiveImplicitStiffness
                ) > stiffnessTolerance
             || mag(bulkRatio - expectedBulkRatio) > ratioTolerance
             || mag(expectedBulkRatio - 1.0) <= ratioTolerance
            )
            {
                FatalErrorInFunction
                    << "Unexpected bulkModulus sensitivity" << nl
                    << "    primary bulkModulus = "
                    << primaryResult.bulkModulus << nl
                    << "    comparison bulkModulus = "
                    << bulkComparisonResult.bulkModulus << nl
                    << "    bulkModulus ratio = " << bulkRatio << nl
                    << "    expected ratio = " << expectedBulkRatio << nl
                    << "    primary effective stiffness = "
                    << primaryResult.effectiveImplicitStiffness << nl
                    << "    comparison effective stiffness = "
                    << bulkComparisonResult.effectiveImplicitStiffness << nl
                    << "    cell stress difference = "
                    << cellStressDifference << nl
                    << "    face stress difference = "
                    << faceStressDifference
                    << abort(FatalError);
            }

            Info<< "PASS: bulkModulus sensitivity" << nl
                << "    primary bulkModulus = "
                << primaryResult.bulkModulus << nl
                << "    comparison bulkModulus = "
                << bulkComparisonResult.bulkModulus << nl
                << "    bulkModulus ratio = " << bulkRatio << nl
                << "    effective stiffness difference = "
                << bulkComparisonResult.effectiveImplicitStiffness
                 - primaryResult.effectiveImplicitStiffness << nl
                << "    cell stress difference = "
                << cellStressDifference << nl
                << "    face stress difference = "
                << faceStressDifference << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
