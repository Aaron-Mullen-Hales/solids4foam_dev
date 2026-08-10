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

\*---------------------------------------------------------------------------*/

#include "GultekinTwoFibreElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//making everything belong to the FOAM namespace
//defining it so that it is a valid nonlineargeom mechanical model
namespace Foam
{
    defineTypeNameAndDebug(GultekinTwoFibreElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GultekinTwoFibreElastic, nonLinGeomMechLaw
    );
}

//local helper functions
namespace
{

using namespace Foam;

//finding the bulk moduluse and returning error if not found or less than equal zero
dimensionedScalar requiredBulkModulus(const dictionary& dict)
{
    if (!dict.found("bulkModulus"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry bulkModulus." << nl
            << "GultekinTwoFibreElastic requires an explicit finite mixed "
            << "bulk modulus with pressure dimensions."
            << exit(FatalIOError);
    }

    const dimensionedScalar result(dict.lookup("bulkModulus"));

    if (result.dimensions() != dimPressure)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus has dimensions " << result.dimensions()
            << "; expected " << dimPressure
            << exit(FatalIOError);
    }

    if (!std::isfinite(result.value()) || result.value() <= 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus must be finite and positive; value = "
            << result.value()
            << exit(FatalIOError);
    }

    return result;
}


word requiredFieldName(const dictionary& dict, const word& keyword)
{
    if (!dict.found(keyword))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry " << keyword
            << exit(FatalIOError);
    }

    return word(dict.lookup(keyword));
}


word optionalSecondFieldName
(
    const dictionary& dict,
    const word& keyword
)
{
    const Switch useSecond
    (
        dict.lookupOrDefault<Switch>("useSecondFibreFamily", true)
    );

    if (useSecond)
    {
        return requiredFieldName(dict, keyword);
    }

    return dict.lookupOrDefault<word>(keyword, word::null);
}


template<class FieldType>
IOobject findReferenceFieldIOobject
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject io
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

#ifdef FOAMEXTEND
    bool ok = io.headerOk();
#elif defined(OPENFOAM_ORG)
    bool ok = io.typeHeaderOk<FieldType>(true);
#else
    bool ok = io.typeHeaderOk<FieldType>(true, false, false);
#endif

    if (!ok)
    {
        io.instance() = "0";

#ifdef FOAMEXTEND
        ok = io.headerOk();
#elif defined(OPENFOAM_ORG)
        ok = io.typeHeaderOk<FieldType>(true);
#else
        ok = io.typeHeaderOk<FieldType>(true, false, false);
#endif
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Cannot find required reference field " << fieldName
            << " in either " << mesh.time().timeName() << " or 0"
            << exit(FatalError);
    }

    io.readOpt() = IOobject::MUST_READ;

    return io;
}


bool finiteVector(const vector& value)
{
    return
        std::isfinite(value.x())
     && std::isfinite(value.y())
     && std::isfinite(value.z());
}


bool finiteSymmTensor(const symmTensor& value)
{
    for (direction cmpt = 0; cmpt < symmTensor::nComponents; ++cmpt)
    {
        if (!std::isfinite(value[cmpt]))
        {
            return false;
        }
    }

    return true;
}


template<class VectorFieldType>
void validateReferenceField
(
    const VectorFieldType& field,
    const scalar unitTolerance
)
{
    if (field.dimensions() != dimless)
    {
        FatalErrorInFunction
            << "Reference fibre field " << field.name()
            << " has dimensions " << field.dimensions()
            << "; expected " << dimless
            << abort(FatalError);
    }

    scalar minMagnitude = GREAT;
    scalar maxMagnitude = -GREAT;
    scalar maxUnitError = 0.0;
    label nonFiniteCount = 0;
    label zeroCount = 0;
    label valueCount = 0;

    forAll(field, valueI)
    {
        const vector& value = field[valueI];
        ++valueCount;

        if (!finiteVector(value))
        {
            ++nonFiniteCount;
            continue;
        }

        const scalar magnitude = mag(value);
        minMagnitude = min(minMagnitude, magnitude);
        maxMagnitude = max(maxMagnitude, magnitude);
        maxUnitError = max(maxUnitError, mag(magnitude - 1.0));

        if (magnitude <= VSMALL)
        {
            ++zeroCount;
        }
    }

    forAll(field.boundaryField(), patchI)
    {
        const vectorField& patch = field.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            const vector& value = patch[faceI];
            ++valueCount;

            if (!finiteVector(value))
            {
                ++nonFiniteCount;
                continue;
            }

            const scalar magnitude = mag(value);
            minMagnitude = min(minMagnitude, magnitude);
            maxMagnitude = max(maxMagnitude, magnitude);
            maxUnitError = max(maxUnitError, mag(magnitude - 1.0));

            if (magnitude <= VSMALL)
            {
                ++zeroCount;
            }
        }
    }

    reduce(minMagnitude, minOp<scalar>());
    reduce(maxMagnitude, maxOp<scalar>());
    reduce(maxUnitError, maxOp<scalar>());
    reduce(nonFiniteCount, sumOp<label>());
    reduce(zeroCount, sumOp<label>());
    reduce(valueCount, sumOp<label>());

    if (!valueCount)
    {
        FatalErrorInFunction
            << "Reference fibre field " << field.name() << " is empty"
            << abort(FatalError);
    }

    Info<< "Reference fibre field " << field.name() << " magnitude:" << nl
        << "    minimum = " << minMagnitude << nl
        << "    maximum = " << maxMagnitude << nl
        << "    maximum unit-length error = " << maxUnitError << endl;

    if (nonFiniteCount)
    {
        FatalErrorInFunction
            << "Reference fibre field " << field.name() << " contains "
            << nonFiniteCount << " NaN or Inf values"
            << abort(FatalError);
    }

    if (zeroCount)
    {
        FatalErrorInFunction
            << "Reference fibre field " << field.name() << " contains "
            << zeroCount << " zero-magnitude values"
            << abort(FatalError);
    }

    if (maxUnitError > unitTolerance)
    {
        FatalErrorInFunction
            << "Reference fibre field " << field.name()
            << " is not unit length" << nl
            << "    minimum magnitude = " << minMagnitude << nl
            << "    maximum magnitude = " << maxMagnitude << nl
            << "    maximum unit-length error = " << maxUnitError << nl
            << "    fibreUnitTolerance = " << unitTolerance << nl
            << "The source data have not been normalised automatically."
            << abort(FatalError);
    }
}


volScalarField* newScalarDiagnostic
(
    const word& name,
    const fvMesh& mesh
)
{
    return new volScalarField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    );
}


volSymmTensorField* newStressDiagnostic
(
    const word& name,
    const fvMesh& mesh
)
{
    return new volSymmTensorField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::GultekinTwoFibreElastic::ConstitutiveResult
Foam::GultekinTwoFibreElastic::evaluateConstitutive
(
    const tensor& F,
    const vector& fibre1,
    const vector& fibre2
) const
{
    ConstitutiveResult result;
    result.J = det(F);
    result.I4 = 0.0;
    result.I6 = 0.0;
    result.lambda4 = 0.0;
    result.lambda6 = 0.0;
    result.sigmaIso = symmTensor::zero;
    result.sigmaFibre1 = symmTensor::zero;
    result.sigmaFibre2 = symmTensor::zero;
    result.sigmaPassive = symmTensor::zero;
    result.exponentClipped = false;
    result.maximumUnclippedExponent = 0.0;

    if (!std::isfinite(result.J) || result.J <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid deformation Jacobian J = " << result.J
            << "; J must be finite and positive"
            << abort(FatalError);
    }

    if (!finiteVector(fibre1) || mag(fibre1) <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid fibre-family-1 direction " << fibre1
            << abort(FatalError);
    }

    if
    (
        useSecondFibreFamily_
     && (!finiteVector(fibre2) || mag(fibre2) <= VSMALL)
    )
    {
        FatalErrorInFunction
            << "Invalid fibre-family-2 direction " << fibre2
            << abort(FatalError);
    }

    const tensor FT(F.T());
    const symmTensor C(symm(FT & F));
    const symmTensor b(symm(F & FT));
    const scalar Jm23 = std::pow(result.J, -2.0/3.0);
    const symmTensor bBar(Jm23*b);

    result.sigmaIso = (mu_.value()/result.J)*dev(bBar);

    const scalar rawI4 = fibre1 & (C & fibre1);

    if (!std::isfinite(rawI4) || rawI4 < 0.0)
    {
        FatalErrorInFunction
            << "Invalid fibre-family-1 invariant I4 = " << rawI4
            << abort(FatalError);
    }

    result.I4 = anisotropicSplit_ ? Jm23*rawI4 : rawI4;

    if (!std::isfinite(result.I4) || result.I4 < 0.0)
    {
        FatalErrorInFunction
            << "Invalid constitutive fibre-family-1 invariant = "
            << result.I4 << abort(FatalError);
    }

    result.lambda4 = std::sqrt(result.I4);

    scalar fibreStrain = result.I4 - 1.0;
    if (fibresTensionOnly_ && fibreStrain < 0.0)
    {
        fibreStrain = 0.0;
    }

    scalar exponent = k2_.value()*sqr(fibreStrain);
    result.maximumUnclippedExponent = exponent;

    if (!std::isfinite(exponent))
    {
        FatalErrorInFunction
            << "Non-finite fibre-family-1 exponential argument " << exponent
            << abort(FatalError);
    }

    if (exponent > exponentLimit_)
    {
        if (!clipExponent_)
        {
            FatalErrorInFunction
                << "Fibre-family-1 exponential argument " << exponent
                << " exceeds exponentLimit " << exponentLimit_ << nl
                << "Set clipExponent true to clip with a warning."
                << abort(FatalError);
        }

        exponent = exponentLimit_;
        result.exponentClipped = true;
    }

    if (fibreStrain != 0.0)
    {
        const vector currentFibre(F & fibre1);
        const symmTensor currentFibreDyad(symm(currentFibre*currentFibre));
        const scalar coefficient =
            2.0*k1_.value()*fibreStrain*std::exp(exponent)/result.J;

        result.sigmaFibre1 =
            anisotropicSplit_
          ? coefficient*Jm23*dev(currentFibreDyad)
          : coefficient*currentFibreDyad;
    }

    if (useSecondFibreFamily_)
    {
        const scalar rawI6 = fibre2 & (C & fibre2);

        if (!std::isfinite(rawI6) || rawI6 < 0.0)
        {
            FatalErrorInFunction
                << "Invalid fibre-family-2 invariant I6 = " << rawI6
                << abort(FatalError);
        }

        result.I6 = anisotropicSplit_ ? Jm23*rawI6 : rawI6;

        if (!std::isfinite(result.I6) || result.I6 < 0.0)
        {
            FatalErrorInFunction
                << "Invalid constitutive fibre-family-2 invariant = "
                << result.I6 << abort(FatalError);
        }

        result.lambda6 = std::sqrt(result.I6);

        fibreStrain = result.I6 - 1.0;
        if (fibresTensionOnly_ && fibreStrain < 0.0)
        {
            fibreStrain = 0.0;
        }

        exponent = k2_.value()*sqr(fibreStrain);
        result.maximumUnclippedExponent =
            max(result.maximumUnclippedExponent, exponent);

        if (!std::isfinite(exponent))
        {
            FatalErrorInFunction
                << "Non-finite fibre-family-2 exponential argument "
                << exponent << abort(FatalError);
        }

        if (exponent > exponentLimit_)
        {
            if (!clipExponent_)
            {
                FatalErrorInFunction
                    << "Fibre-family-2 exponential argument " << exponent
                    << " exceeds exponentLimit " << exponentLimit_ << nl
                    << "Set clipExponent true to clip with a warning."
                    << abort(FatalError);
            }

            exponent = exponentLimit_;
            result.exponentClipped = true;
        }

        if (fibreStrain != 0.0)
        {
            const vector currentFibre(F & fibre2);
            const symmTensor currentFibreDyad
            (
                symm(currentFibre*currentFibre)
            );
            const scalar coefficient =
                2.0*k1_.value()*fibreStrain*std::exp(exponent)/result.J;

            result.sigmaFibre2 =
                anisotropicSplit_
              ? coefficient*Jm23*dev(currentFibreDyad)
              : coefficient*currentFibreDyad;
        }
    }

    result.sigmaPassive =
        result.sigmaIso + result.sigmaFibre1 + result.sigmaFibre2;

    if
    (
        !finiteSymmTensor(result.sigmaIso)
     || !finiteSymmTensor(result.sigmaFibre1)
     || !finiteSymmTensor(result.sigmaFibre2)
     || !finiteSymmTensor(result.sigmaPassive)
    )
    {
        FatalErrorInFunction
            << "Non-finite Cauchy stress generated for J = " << result.J
            << ", I4 = " << result.I4
            << ", I6 = " << result.I6
            << abort(FatalError);
    }

    return result;
}


void Foam::GultekinTwoFibreElastic::validateMaterialParameters() const
{
    if (rho_.dimensions() != dimDensity)
    {
        FatalErrorInFunction
            << "rho has dimensions " << rho_.dimensions()
            << "; expected " << dimDensity
            << abort(FatalError);
    }

    if
    (
        mu_.dimensions() != dimPressure
     || k1_.dimensions() != dimPressure
     || bulkModulus_.dimensions() != dimPressure
     || unscaledImplicitShearModulus_.dimensions() != dimPressure
     || implicitShearModulus_.dimensions() != dimPressure
    )
    {
        FatalErrorInFunction
            << "mu, k1, bulkModulus and both implicit stiffness values must "
            << "have pressure dimensions " << dimPressure
            << abort(FatalError);
    }

    if (k2_.dimensions() != dimless)
    {
        FatalErrorInFunction
            << "k2 has dimensions " << k2_.dimensions()
            << "; expected " << dimless
            << abort(FatalError);
    }

    if (!std::isfinite(rho_.value()) || rho_.value() <= 0.0)
    {
        FatalErrorInFunction
            << "rho must be finite and positive; rho = " << rho_.value()
            << abort(FatalError);
    }

    if (!std::isfinite(mu_.value()) || mu_.value() < 0.0)
    {
        FatalErrorInFunction
            << "mu must be finite and non-negative; mu = " << mu_.value()
            << abort(FatalError);
    }

    if (!std::isfinite(k1_.value()) || k1_.value() < 0.0)
    {
        FatalErrorInFunction
            << "k1 must be finite and non-negative; k1 = " << k1_.value()
            << abort(FatalError);
    }

    if
    (
        !std::isfinite(bulkModulus_.value())
     || bulkModulus_.value() <= 0.0
    )
    {
        FatalErrorInFunction
            << "bulkModulus must be finite and positive; value = "
            << bulkModulus_.value()
            << abort(FatalError);
    }

    if (!std::isfinite(k2_.value()) || k2_.value() <= 0.0)
    {
        FatalErrorInFunction
            << "k2 must be finite and positive; k2 = " << k2_.value()
            << abort(FatalError);
    }

    if
    (
        !std::isfinite(unscaledImplicitShearModulus_.value())
     || unscaledImplicitShearModulus_.value() <= 0.0
    )
    {
        FatalErrorInFunction
            << "implicitShearModulus must be finite and positive; value = "
            << unscaledImplicitShearModulus_.value()
            << abort(FatalError);
    }

    if
    (
        !std::isfinite(implicitShearModulus_.value())
     || implicitShearModulus_.value() <= 0.0
    )
    {
        FatalErrorInFunction
            << "The numerical implicit stiffness must be finite and positive; "
            << "value = " << implicitShearModulus_.value()
            << abort(FatalError);
    }

    if (!std::isfinite(impKcoeff_) || impKcoeff_ <= 0.0)
    {
        FatalErrorInFunction
            << "impKcoeff must be finite and positive; impKcoeff = "
            << impKcoeff_ << abort(FatalError);
    }

    if (!std::isfinite(fibreUnitTolerance_) || fibreUnitTolerance_ <= 0.0)
    {
        FatalErrorInFunction
            << "fibreUnitTolerance must be finite and positive; value = "
            << fibreUnitTolerance_ << abort(FatalError);
    }

    if
    (
        !std::isfinite(exponentLimit_)
     || exponentLimit_ <= 0.0
     || exponentLimit_ > 650.0
    )
    {
        FatalErrorInFunction
            << "exponentLimit must be finite and in (0, 650]; value = "
            << exponentLimit_ << abort(FatalError);
    }
}


void Foam::GultekinTwoFibreElastic::validateFibreFields() const
{
    if (fibreField1Name_ == faceFibreField1Name_)
    {
        FatalErrorInFunction
            << "fibreField1 and faceFibreField1 must name distinct volume "
            << "and surface fields; both are " << fibreField1Name_
            << abort(FatalError);
    }

    if (useSecondFibreFamily_)
    {
        if (fibreField1Name_ == fibreField2Name_)
        {
            FatalErrorInFunction
                << "fibreField1 and fibreField2 must be distinct"
                << abort(FatalError);
        }

        if (faceFibreField1Name_ == faceFibreField2Name_)
        {
            FatalErrorInFunction
                << "faceFibreField1 and faceFibreField2 must be distinct"
                << abort(FatalError);
        }

        if (fibreField2Name_ == faceFibreField2Name_)
        {
            FatalErrorInFunction
                << "fibreField2 and faceFibreField2 must name distinct volume "
                << "and surface fields; both are " << fibreField2Name_
                << abort(FatalError);
        }
    }

    validateReferenceField(fibre1_, fibreUnitTolerance_);
    validateReferenceField(faceFibre1_, fibreUnitTolerance_);

    if (useSecondFibreFamily_)
    {
        validateReferenceField(fibre2Ptr_(), fibreUnitTolerance_);
        validateReferenceField(faceFibre2Ptr_(), fibreUnitTolerance_);
    }
}


void Foam::GultekinTwoFibreElastic::makeDiagnostics()
{
    if (!writeDiagnostics_)
    {
        return;
    }

    diagnosticJPtr_.set(newScalarDiagnostic("GTF_J", mesh()));
    diagnosticI4Ptr_.set(newScalarDiagnostic("GTF_I4", mesh()));
    diagnosticI6Ptr_.set(newScalarDiagnostic("GTF_I6", mesh()));
    diagnosticLambda4Ptr_.set(newScalarDiagnostic("GTF_lambda4", mesh()));
    diagnosticLambda6Ptr_.set(newScalarDiagnostic("GTF_lambda6", mesh()));
    diagnosticSigmaIsoPtr_.set(newStressDiagnostic("GTF_sigmaIso", mesh()));
    diagnosticSigmaFibre1Ptr_.set
    (
        newStressDiagnostic("GTF_sigmaFibre1", mesh())
    );
    diagnosticSigmaFibre2Ptr_.set
    (
        newStressDiagnostic("GTF_sigmaFibre2", mesh())
    );
    diagnosticSigmaPassivePtr_.set
    (
        newStressDiagnostic("GTF_sigmaPassive", mesh())
    );
}


void Foam::GultekinTwoFibreElastic::reportExponentClipping
(
    const word& location,
    label clippedCount,
    scalar maximumUnclippedExponent
) const
{
    reduce(clippedCount, sumOp<label>());
    reduce(maximumUnclippedExponent, maxOp<scalar>());

    if (clippedCount)
    {
        WarningInFunction
            << "Clipped the fibre exponential argument at " << location
            << nl
            << "    clipped values = " << clippedCount << nl
            << "    maximum unclipped argument = "
            << maximumUnclippedExponent << nl
            << "    configured exponentLimit = " << exponentLimit_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GultekinTwoFibreElastic::GultekinTwoFibreElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu_(dict.lookup("mu")),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2")),
    bulkModulus_(requiredBulkModulus(dict)),
    anisotropicSplit_
    (
        dict.lookupOrDefault<Switch>("anisotropicSplit", false)
    ),
    fibresTensionOnly_
    (
        dict.lookupOrDefault<Switch>("fibresTensionOnly", false)
    ),
    useSecondFibreFamily_
    (
        dict.lookupOrDefault<Switch>("useSecondFibreFamily", true)
    ),
    fibreField1Name_(requiredFieldName(dict, "fibreField1")),
    fibreField2Name_(optionalSecondFieldName(dict, "fibreField2")),
    faceFibreField1Name_(requiredFieldName(dict, "faceFibreField1")),
    faceFibreField2Name_
    (
        optionalSecondFieldName(dict, "faceFibreField2")
    ),
    fibre1_
    (
        findReferenceFieldIOobject<volVectorField>(fibreField1Name_, mesh),
        mesh
    ),
    fibre2Ptr_(),
    faceFibre1_
    (
        findReferenceFieldIOobject<surfaceVectorField>
        (
            faceFibreField1Name_, mesh
        ),
        mesh
    ),
    faceFibre2Ptr_(),
    impKcoeff_(dict.lookupOrDefault<scalar>("impKcoeff", 1.0)),
    implicitShearModulusSpecified_(dict.found("implicitShearModulus")),
    unscaledImplicitShearModulus_
    (
        "implicitShearModulus",
        dimPressure,
        0.0
    ),
    implicitShearModulus_("implicitShearModulus", dimPressure, 0.0),
    fibreUnitTolerance_
    (
        dict.lookupOrDefault<scalar>("fibreUnitTolerance", 1e-6)
    ),
    clipExponent_(dict.lookupOrDefault<Switch>("clipExponent", true)),
    exponentLimit_(dict.lookupOrDefault<scalar>("exponentLimit", 650.0)),
    writeDiagnostics_
    (
        dict.lookupOrDefault<Switch>("writeDiagnostics", false)
    ),
    diagnosticJPtr_(),
    diagnosticI4Ptr_(),
    diagnosticI6Ptr_(),
    diagnosticLambda4Ptr_(),
    diagnosticLambda6Ptr_(),
    diagnosticSigmaIsoPtr_(),
    diagnosticSigmaFibre1Ptr_(),
    diagnosticSigmaFibre2Ptr_(),
    diagnosticSigmaPassivePtr_()
{
    if (useSecondFibreFamily_)
    {
        fibre2Ptr_.set
        (
            new volVectorField
            (
                findReferenceFieldIOobject<volVectorField>
                (
                    fibreField2Name_, mesh
                ),
                mesh
            )
        );

        faceFibre2Ptr_.set
        (
            new surfaceVectorField
            (
                findReferenceFieldIOobject<surfaceVectorField>
                (
                    faceFibreField2Name_, mesh
                ),
                mesh
            )
        );
    }

    if (implicitShearModulusSpecified_)
    {
        unscaledImplicitShearModulus_ =
            dimensionedScalar(dict.lookup("implicitShearModulus"));
    }
    else
    {
        const scalar familyCount = useSecondFibreFamily_ ? 2.0 : 1.0;
        unscaledImplicitShearModulus_ =
            mu_ + 2.0*familyCount*k1_;
    }

    implicitShearModulus_ =
        impKcoeff_*unscaledImplicitShearModulus_;

    validateMaterialParameters();
    validateFibreFields();
    makeDiagnostics();

    F().storeOldTime();
    Ff().storeOldTime();

    Info<< "GultekinTwoFibreElastic options:" << nl
        << "    anisotropicSplit = " << anisotropicSplit_ << nl
        << "    fibresTensionOnly = " << fibresTensionOnly_ << nl
        << "    useSecondFibreFamily = " << useSecondFibreFamily_ << nl
        << "    fibreField1 = " << fibreField1Name_ << nl
        << "    faceFibreField1 = " << faceFibreField1Name_ << nl;

    if (useSecondFibreFamily_)
    {
        Info<< "    fibreField2 = " << fibreField2Name_ << nl
            << "    faceFibreField2 = " << faceFibreField2Name_ << nl;
    }

    const dimensionedScalar reciprocalBulkModulus
    (
        "reciprocalBulkModulus",
        dimless/dimPressure,
        1.0/bulkModulus_.value()
    );

    Info<< "    bulkModulus = " << bulkModulus_ << nl
        << "    reciprocalBulkModulus = " << reciprocalBulkModulus << nl;

    if (mu_.value() > 0.0)
    {
        Info<< "    bulk-to-matrix-shear ratio = "
            << bulkModulus_.value()/mu_.value() << nl;
    }
    else
    {
        Info<< "    bulk-to-matrix-shear ratio = undefined (mu = 0)" << nl;
    }

    Info<< "    implicitShearModulus source = "
        <<
        (
            implicitShearModulusSpecified_
          ? "dictionary"
          : "default estimate"
        )
        << nl
        << "    implicitShearModulus = "
        << unscaledImplicitShearModulus_ << nl
        << "    impKcoeff = " << impKcoeff_ << nl
        << "    effective implicit stiffness = "
        << implicitShearModulus_ << nl
        << "    clipExponent = " << clipExponent_ << nl
        << "    exponentLimit = " << exponentLimit_ << nl
        << "    writeDiagnostics = " << writeDiagnostics_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GultekinTwoFibreElastic::~GultekinTwoFibreElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GultekinTwoFibreElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::GultekinTwoFibreElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            implicitShearModulus_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::GultekinTwoFibreElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulusLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            bulkModulus_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::GultekinTwoFibreElastic::shearModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "shearModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            implicitShearModulus_
        )
    );
}


void Foam::GultekinTwoFibreElastic::correct(volSymmTensorField& sigma)
{
    dimensionedScalar zeroBulkModulus("K", dimPressure, 0.0);
    if (updateF(sigma, mu_, zeroBulkModulus))
    {
        return;
    }

    const volTensorField& deformationGradient = F();
    label clippedCount = 0;
    scalar maximumUnclippedExponent = 0.0;

    forAll(sigma, cellI)
    {
        const vector fibre2 =
            useSecondFibreFamily_ ? fibre2Ptr_()[cellI] : vector::zero;
        const ConstitutiveResult result = evaluateConstitutive
        (
            deformationGradient[cellI], fibre1_[cellI], fibre2
        );

        sigma[cellI] = result.sigmaPassive;
        clippedCount += result.exponentClipped ? 1 : 0;
        maximumUnclippedExponent = max
        (
            maximumUnclippedExponent,
            result.maximumUnclippedExponent
        );

        if (writeDiagnostics_)
        {
            diagnosticJPtr_()[cellI] = result.J;
            diagnosticI4Ptr_()[cellI] = result.I4;
            diagnosticI6Ptr_()[cellI] = result.I6;
            diagnosticLambda4Ptr_()[cellI] = result.lambda4;
            diagnosticLambda6Ptr_()[cellI] = result.lambda6;
            diagnosticSigmaIsoPtr_()[cellI] = result.sigmaIso;
            diagnosticSigmaFibre1Ptr_()[cellI] = result.sigmaFibre1;
            diagnosticSigmaFibre2Ptr_()[cellI] = result.sigmaFibre2;
            diagnosticSigmaPassivePtr_()[cellI] = result.sigmaPassive;
        }
    }

    forAll(sigma.boundaryField(), patchI)
    {
        symmTensorField& sigmaPatch = sigma.boundaryFieldRef()[patchI];
        const tensorField& deformationGradientPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& fibre1Patch = fibre1_.boundaryField()[patchI];

        forAll(sigmaPatch, faceI)
        {
            const vector fibre2 =
                useSecondFibreFamily_
              ? fibre2Ptr_().boundaryField()[patchI][faceI]
              : vector::zero;
            const ConstitutiveResult result = evaluateConstitutive
            (
                deformationGradientPatch[faceI],
                fibre1Patch[faceI],
                fibre2
            );

            sigmaPatch[faceI] = result.sigmaPassive;
            clippedCount += result.exponentClipped ? 1 : 0;
            maximumUnclippedExponent = max
            (
                maximumUnclippedExponent,
                result.maximumUnclippedExponent
            );

            if (writeDiagnostics_)
            {
                diagnosticJPtr_().boundaryFieldRef()[patchI][faceI] = result.J;
                diagnosticI4Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.I4;
                diagnosticI6Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.I6;
                diagnosticLambda4Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.lambda4;
                diagnosticLambda6Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.lambda6;
                diagnosticSigmaIsoPtr_().boundaryFieldRef()[patchI][faceI] =
                    result.sigmaIso;
                diagnosticSigmaFibre1Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.sigmaFibre1;
                diagnosticSigmaFibre2Ptr_().boundaryFieldRef()[patchI][faceI] =
                    result.sigmaFibre2;
                diagnosticSigmaPassivePtr_().boundaryFieldRef()[patchI][faceI] =
                    result.sigmaPassive;
            }
        }
    }

    reportExponentClipping
    (
        "cells and cell boundaries",
        clippedCount,
        maximumUnclippedExponent
    );
}


void Foam::GultekinTwoFibreElastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    dimensionedScalar zeroBulkModulus("K", dimPressure, 0.0);
    if (updateF(sigma, mu_, zeroBulkModulus))
    {
        return;
    }

    const surfaceTensorField& deformationGradient = Ff();
    label clippedCount = 0;
    scalar maximumUnclippedExponent = 0.0;

    forAll(sigma, faceI)
    {
        const vector fibre2 =
            useSecondFibreFamily_ ? faceFibre2Ptr_()[faceI] : vector::zero;
        const ConstitutiveResult result = evaluateConstitutive
        (
            deformationGradient[faceI], faceFibre1_[faceI], fibre2
        );

        sigma[faceI] = result.sigmaPassive;
        clippedCount += result.exponentClipped ? 1 : 0;
        maximumUnclippedExponent = max
        (
            maximumUnclippedExponent,
            result.maximumUnclippedExponent
        );
    }

    forAll(sigma.boundaryField(), patchI)
    {
        symmTensorField& sigmaPatch = sigma.boundaryFieldRef()[patchI];
        const tensorField& deformationGradientPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& fibre1Patch =
            faceFibre1_.boundaryField()[patchI];

        forAll(sigmaPatch, faceI)
        {
            const vector fibre2 =
                useSecondFibreFamily_
              ? faceFibre2Ptr_().boundaryField()[patchI][faceI]
              : vector::zero;
            const ConstitutiveResult result = evaluateConstitutive
            (
                deformationGradientPatch[faceI],
                fibre1Patch[faceI],
                fibre2
            );

            sigmaPatch[faceI] = result.sigmaPassive;
            clippedCount += result.exponentClipped ? 1 : 0;
            maximumUnclippedExponent = max
            (
                maximumUnclippedExponent,
                result.maximumUnclippedExponent
            );
        }
    }

    reportExponentClipping
    (
        "faces and face boundaries",
        clippedCount,
        maximumUnclippedExponent
    );
}


void Foam::GultekinTwoFibreElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}


// ************************************************************************* //
