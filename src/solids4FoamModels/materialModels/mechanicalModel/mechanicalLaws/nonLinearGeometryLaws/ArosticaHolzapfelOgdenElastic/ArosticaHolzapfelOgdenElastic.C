/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ArosticaHolzapfelOgdenElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ArosticaHolzapfelOgdenElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, ArosticaHolzapfelOgdenElastic, nonLinGeomMechLaw
    );
}

namespace
{

using namespace Foam;


word requiredFieldName
(
    const dictionary& dict,
    const word& keyword
)
{
    if (!dict.found(keyword))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry " << keyword
            << exit(FatalIOError);
    }

    return word(dict.lookup(keyword));
}


template<class FieldType>
IOobject requiredReferenceFieldIOobject
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
            << "Cannot find required Aróstica reference field " << fieldName
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
    for (direction cmptI = 0; cmptI < symmTensor::nComponents; ++cmptI)
    {
        if (!std::isfinite(value[cmptI]))
        {
            return false;
        }
    }

    return true;
}


void validateTriadValue
(
    const vector& f,
    const vector& s,
    const vector& n,
    const scalar tolerance,
    const word& location,
    const label valueI
)
{
    if
    (
        !finiteVector(f)
     || !finiteVector(s)
     || !finiteVector(n)
     || mag(f) <= VSMALL
     || mag(s) <= VSMALL
     || mag(n) <= VSMALL
     || mag(mag(f) - 1.0) > tolerance
     || mag(mag(s) - 1.0) > tolerance
     || mag(mag(n) - 1.0) > tolerance
     || mag(f & s) > tolerance
     || mag(f & n) > tolerance
     || mag(s & n) > tolerance
    )
    {
        FatalErrorInFunction
            << "Invalid supplied Aróstica " << location << " triad at "
            << valueI << nl
            << "    f0 = " << f << ", s0 = " << s << ", n0 = " << n
            << nl
            << "    magnitudes = " << mag(f) << ", " << mag(s) << ", "
            << mag(n) << nl
            << "    dot products = " << (f & s) << ", " << (f & n)
            << ", " << (s & n)
            << exit(FatalError);
    }
}


template<class FieldType>
void validateTriadField
(
    const FieldType& fibre,
    const FieldType& sheet,
    const FieldType& normal,
    const scalar tolerance,
    const word& location
)
{
    if
    (
        fibre.dimensions() != dimless
     || sheet.dimensions() != dimless
     || normal.dimensions() != dimless
    )
    {
        FatalErrorInFunction
            << "Aróstica " << location << " triad fields must be dimensionless"
            << exit(FatalError);
    }

    forAll(fibre, valueI)
    {
        validateTriadValue
        (
            fibre[valueI],
            sheet[valueI],
            normal[valueI],
            tolerance,
            location,
            valueI
        );
    }

    forAll(fibre.boundaryField(), patchI)
    {
        const auto& fibrePatch = fibre.boundaryField()[patchI];
        const auto& sheetPatch = sheet.boundaryField()[patchI];
        const auto& normalPatch = normal.boundaryField()[patchI];

        forAll(fibrePatch, faceI)
        {
            validateTriadValue
            (
                fibrePatch[faceI],
                sheetPatch[faceI],
                normalPatch[faceI],
                tolerance,
                location + " boundary patch " + Foam::name(patchI),
                faceI
            );
        }
    }
}


void updateTriadDiagnostics
(
    const vector& f,
    const vector& s,
    const vector& n,
    scalar& minF,
    scalar& maxF,
    scalar& minS,
    scalar& maxS,
    scalar& minN,
    scalar& maxN,
    scalar& maxFS,
    scalar& maxFN,
    scalar& maxSN,
    scalar& minDet,
    scalar& maxDet
)
{
    minF = min(minF, mag(f));
    maxF = max(maxF, mag(f));
    minS = min(minS, mag(s));
    maxS = max(maxS, mag(s));
    minN = min(minN, mag(n));
    maxN = max(maxN, mag(n));
    maxFS = max(maxFS, mag(f & s));
    maxFN = max(maxFN, mag(f & n));
    maxSN = max(maxSN, mag(s & n));

    const scalar determinant = (f ^ s) & n;
    minDet = min(minDet, determinant);
    maxDet = max(maxDet, determinant);
}


template<class FieldType>
void reportTriadDiagnostics
(
    const FieldType& fibre,
    const FieldType& sheet,
    const FieldType& normal,
    const word& location
)
{
    scalar minF = GREAT;
    scalar maxF = -GREAT;
    scalar minS = GREAT;
    scalar maxS = -GREAT;
    scalar minN = GREAT;
    scalar maxN = -GREAT;
    scalar maxFS = 0.0;
    scalar maxFN = 0.0;
    scalar maxSN = 0.0;
    scalar minDet = GREAT;
    scalar maxDet = -GREAT;

    forAll(fibre, valueI)
    {
        updateTriadDiagnostics
        (
            fibre[valueI],
            sheet[valueI],
            normal[valueI],
            minF,
            maxF,
            minS,
            maxS,
            minN,
            maxN,
            maxFS,
            maxFN,
            maxSN,
            minDet,
            maxDet
        );
    }

    forAll(fibre.boundaryField(), patchI)
    {
        const auto& fibrePatch = fibre.boundaryField()[patchI];
        const auto& sheetPatch = sheet.boundaryField()[patchI];
        const auto& normalPatch = normal.boundaryField()[patchI];

        forAll(fibrePatch, faceI)
        {
            updateTriadDiagnostics
            (
                fibrePatch[faceI],
                sheetPatch[faceI],
                normalPatch[faceI],
                minF,
                maxF,
                minS,
                maxS,
                minN,
                maxN,
                maxFS,
                maxFN,
                maxSN,
                minDet,
                maxDet
            );
        }
    }

    Info<< "Aróstica " << location << " triad diagnostics:" << nl
        << "    |f0| range = " << minF << " to " << maxF << nl
        << "    |s0| range = " << minS << " to " << maxS << nl
        << "    |n0| range = " << minN << " to " << maxN << nl
        << "    max |f0.s0| = " << maxFS << nl
        << "    max |f0.n0| = " << maxFN << nl
        << "    max |s0.n0| = " << maxSN << nl
        << "    determinant range = " << minDet << " to " << maxDet
        << endl;
}


scalar validatePositiveFinite
(
    const dictionary& dict,
    const word& keyword
)
{
    const scalar value = dict.lookupOrDefault<scalar>(keyword, -GREAT);

    if (!std::isfinite(value) || value <= 0.0)
    {
        FatalIOErrorInFunction(dict)
            << keyword << " must be finite and positive; value = " << value
            << exit(FatalIOError);
    }

    return value;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::ArosticaHolzapfelOgdenElastic::CompressionSwitchResult
Foam::ArosticaHolzapfelOgdenElastic::compressionSwitchValue
(
    const scalar I4
) const
{
    CompressionSwitchResult result;

    if (compressionSwitch_ == "paperPiecewise")
    {
        if (I4 > 1.0)
        {
            result.chi = I4;
            result.dChi = 1.0;
        }
        else
        {
            result.chi = 0.0;
            result.dChi = 0.0;
        }
    }
    else if (compressionSwitch_ == "logistic")
    {
        const scalar exponent = -compressionSwitchK_*(I4 - 1.0);
        const scalar logistic = 1.0/(1.0 + std::exp(exponent));

        result.chi = logistic;
        result.dChi = compressionSwitchK_*logistic*(1.0 - logistic);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown compressionSwitch " << compressionSwitch_ << nl
            << "Valid options are paperPiecewise and logistic"
            << exit(FatalError);
    }

    return result;
}


Foam::ArosticaHolzapfelOgdenElastic::ConstitutiveResult
Foam::ArosticaHolzapfelOgdenElastic::evaluateConstitutive
(
    const tensor& F,
    const vector& f0,
    const vector& s0,
    const vector& n0
) const
{
    ConstitutiveResult result;
    result.J = det(F);

    if (!std::isfinite(result.J) || result.J <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid deformation Jacobian J = " << result.J
            << "; J must be finite and positive"
            << abort(FatalError);
    }

    if
    (
        !finiteVector(f0)
     || !finiteVector(s0)
     || !finiteVector(n0)
    )
    {
        FatalErrorInFunction
            << "Non-finite Aróstica reference triad supplied to constitutive "
            << "kernel"
            << abort(FatalError);
    }

    const symmTensor C(symm(F.T() & F));
    const scalar traceC = tr(C);
    const scalar Jm23 = std::pow(result.J, -2.0/3.0);
    const symmTensor Cinv(inv(C));

    result.I1bar = Jm23*traceC;
    result.I4f = f0 & (C & f0);
    result.I4s = s0 & (C & s0);
    result.I8fs = f0 & (C & s0);

    if
    (
        !std::isfinite(result.I1bar)
     || !std::isfinite(result.I4f)
     || !std::isfinite(result.I4s)
     || !std::isfinite(result.I8fs)
    )
    {
        FatalErrorInFunction
            << "Non-finite Aróstica invariant generated for J = "
            << result.J
            << abort(FatalError);
    }

    const scalar W1 =
        a_.value()/(2.0*b_.value())
       *(std::exp(b_.value()*(result.I1bar - 3.0)) - 1.0);

    const CompressionSwitchResult fibreSwitch =
        compressionSwitchValue(result.I4f);
    const CompressionSwitchResult sheetSwitch =
        compressionSwitchValue(result.I4s);

    const scalar fibreExponent =
        std::exp(bf_.value()*sqr(result.I4f - 1.0));
    const scalar sheetExponent =
        std::exp(bs_.value()*sqr(result.I4s - 1.0));
    const scalar fibreSheetExponent =
        std::exp(bfs_.value()*sqr(result.I8fs));

    result.energy =
        W1
      + af_.value()/(2.0*bf_.value())
       *fibreSwitch.chi*(fibreExponent - 1.0)
      + as_.value()/(2.0*bs_.value())
       *sheetSwitch.chi*(sheetExponent - 1.0)
      + afs_.value()/(2.0*bfs_.value())
       *(fibreSheetExponent - 1.0);

    const symmTensor dI1bar
    (
        Jm23
       *
        (
            symmTensor::I
          - (traceC/3.0)*Cinv
        )
    );

    symmTensor S
    (
        a_.value()*std::exp(b_.value()*(result.I1bar - 3.0))*dI1bar
    );

    const scalar qf = result.I4f - 1.0;
    const scalar qs = result.I4s - 1.0;
    const scalar dWf =
        af_.value()/(2.0*bf_.value())
       *
        (
            fibreSwitch.dChi*(fibreExponent - 1.0)
          + fibreSwitch.chi*fibreExponent*(2.0*bf_.value()*qf)
        );
    const scalar dWs =
        as_.value()/(2.0*bs_.value())
       *
        (
            sheetSwitch.dChi*(sheetExponent - 1.0)
          + sheetSwitch.chi*sheetExponent*(2.0*bs_.value()*qs)
        );

    S += 2.0*dWf*symm(f0*f0);
    S += 2.0*dWs*symm(s0*s0);

    const scalar dW8 =
        afs_.value()*result.I8fs*fibreSheetExponent;
    S += 2.0*dW8*symm(f0*s0);

    result.S = S;
    result.sigma = symm(F & S & F.T())/result.J;

    if
    (
        !finiteSymmTensor(result.S)
     || !finiteSymmTensor(result.sigma)
     || !std::isfinite(result.energy)
    )
    {
        FatalErrorInFunction
            << "Non-finite Aróstica stress or energy generated for J = "
            << result.J
            << abort(FatalError);
    }

    return result;
}


void Foam::ArosticaHolzapfelOgdenElastic::validateMaterialParameters() const
{
    if
    (
        rho_.dimensions() != dimDensity
     || a_.dimensions() != dimPressure
     || af_.dimensions() != dimPressure
     || as_.dimensions() != dimPressure
     || afs_.dimensions() != dimPressure
     || bulkModulus_.dimensions() != dimPressure
    )
    {
        FatalErrorInFunction
            << "rho, a, af, as, afs and bulkModulus must have dimensions "
            << "density or pressure as appropriate"
            << abort(FatalError);
    }

    if
    (
        b_.dimensions() != dimless
     || bf_.dimensions() != dimless
     || bs_.dimensions() != dimless
     || bfs_.dimensions() != dimless
    )
    {
        FatalErrorInFunction
            << "b, bf, bs and bfs must be dimensionless"
            << abort(FatalError);
    }

    if
    (
        !std::isfinite(rho_.value()) || rho_.value() <= 0.0
     || !std::isfinite(a_.value()) || a_.value() < 0.0
     || !std::isfinite(af_.value()) || af_.value() < 0.0
     || !std::isfinite(as_.value()) || as_.value() < 0.0
     || !std::isfinite(afs_.value()) || afs_.value() < 0.0
     || !std::isfinite(b_.value()) || b_.value() <= 0.0
     || !std::isfinite(bf_.value()) || bf_.value() <= 0.0
     || !std::isfinite(bs_.value()) || bs_.value() <= 0.0
     || !std::isfinite(bfs_.value()) || bfs_.value() <= 0.0
     || !std::isfinite(bulkModulus_.value())
     || bulkModulus_.value() <= 0.0
     || !std::isfinite(impKcoeff_) || impKcoeff_ <= 0.0
     || !std::isfinite(compressionSwitchK_)
     || compressionSwitchK_ <= 0.0
     || !std::isfinite(directionTolerance_)
     || directionTolerance_ <= 0.0
     || !std::isfinite(tangentEps_)
     || tangentEps_ <= 0.0
    )
    {
        FatalErrorInFunction
            << "Aróstica material values must be finite and positive where "
            << "required"
            << abort(FatalError);
    }

    if
    (
        compressionSwitch_ != "paperPiecewise"
     && compressionSwitch_ != "logistic"
    )
    {
        FatalErrorInFunction
            << "Unknown compressionSwitch " << compressionSwitch_ << nl
            << "Valid options are paperPiecewise and logistic"
            << abort(FatalError);
    }
}


void Foam::ArosticaHolzapfelOgdenElastic::validateReferenceFields() const
{
    validateTriadField
    (
        f0_, s0_, n0_, directionTolerance_, "cell"
    );
    validateTriadField
    (
        f0f_, s0f_, n0f_, directionTolerance_, "face"
    );

    reportTriadDiagnostics(f0_, s0_, n0_, "cell");
    reportTriadDiagnostics(f0f_, s0f_, n0f_, "face");
}


void Foam::ArosticaHolzapfelOgdenElastic::calculateRepresentativeShearModulus()
{
    if (!f0_.size())
    {
        FatalErrorInFunction
            << "Cannot calculate representative shear modulus without cells"
            << exit(FatalError);
    }

    const vector axes[3] = {f0_[0], s0_[0], n0_[0]};
    const scalar perturbation = 1.0e-6;
    scalar sampledModulus[3];

    const label pairA[3] = {0, 1, 2};
    const label pairB[3] = {1, 2, 0};

    for (label sampleI = 0; sampleI < 3; ++sampleI)
    {
        const vector& axisA = axes[pairA[sampleI]];
        const vector& axisB = axes[pairB[sampleI]];
        const tensor Ftest
        (
            tensor::I
          + perturbation*(axisA*axisB + axisB*axisA)
        );
        const ConstitutiveResult result =
            evaluateConstitutive(Ftest, axes[0], axes[1], axes[2]);

        sampledModulus[sampleI] =
            mag(axisA & (result.sigma & axisB))
           /(2.0*perturbation);

        if
        (
            !std::isfinite(sampledModulus[sampleI])
         || sampledModulus[sampleI] <= VSMALL
        )
        {
            FatalErrorInFunction
                << "Non-positive sampled local shear modulus "
                << sampledModulus[sampleI]
                << " for pair " << pairA[sampleI] << "-"
                << pairB[sampleI]
                << exit(FatalError);
        }
    }

    const scalar selected =
        max(sampledModulus[0], max(sampledModulus[1], sampledModulus[2]));

    muRepresentative_ =
        dimensionedScalar("muRepresentative", dimPressure, selected);

    Info<< "Aróstica local shear-modulus samples:" << nl
        << "    f0-s0 = " << sampledModulus[0] << nl
        << "    s0-n0 = " << sampledModulus[1] << nl
        << "    n0-f0 = " << sampledModulus[2] << nl
        << "    selected representative = " << muRepresentative_ << nl
        << "    impKcoeff = " << impKcoeff_ << endl;
}


void Foam::ArosticaHolzapfelOgdenElastic::fillTangent
(
    mat66& tangent,
    const tensor& F,
    const vector& f0,
    const vector& s0,
    const vector& n0
) const
{
    tangent.clear();

    const label tensorComponent[6] =
    {
        tensor::XX,
        tensor::YY,
        tensor::ZZ,
        tensor::XY,
        tensor::YZ,
        tensor::XZ
    };

    for (label componentI = 0; componentI < symmTensor::nComponents; ++componentI)
    {
        tensor Fplus(F);
        tensor Fminus(F);
        Fplus[tensorComponent[componentI]] += tangentEps_;
        Fminus[tensorComponent[componentI]] -= tangentEps_;

        const ConstitutiveResult plus =
            evaluateConstitutive(Fplus, f0, s0, n0);
        const ConstitutiveResult minus =
            evaluateConstitutive(Fminus, f0, s0, n0);
        const symmTensor derivative
        (
            (plus.sigma - minus.sigma)/(2.0*tangentEps_)
        );

        tangent(symmTensor::XX, componentI) = derivative. xx();
        tangent(symmTensor::YY, componentI) = derivative. yy();
        tangent(symmTensor::ZZ, componentI) = derivative. zz();
        tangent(symmTensor::XY, componentI) = derivative. xy();
        tangent(symmTensor::YZ, componentI) = derivative.yz();
        tangent(symmTensor::XZ, componentI) = derivative.xz();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ArosticaHolzapfelOgdenElastic::ArosticaHolzapfelOgdenElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    a_(dict.lookup("a")),
    b_(dict.lookup("b")),
    af_(dict.lookup("af")),
    bf_(dict.lookup("bf")),
    as_(dict.lookup("as")),
    bs_(dict.lookup("bs")),
    afs_(dict.lookup("afs")),
    bfs_(dict.lookup("bfs")),
    bulkModulus_(dict.lookup("bulkModulus")),
    impKcoeff_(validatePositiveFinite(dict, "impKcoeff")),
    compressionSwitch_(word(dict.lookup("compressionSwitch"))),
    compressionSwitchK_(validatePositiveFinite(dict, "compressionSwitchK")),
    directionTolerance_
    (
        dict.lookupOrDefault<scalar>("directionTolerance", 1.0e-6)
    ),
    tangentEps_
    (
        dict.lookupOrDefault<scalar>("tangentEps", 1.0e-7)
    ),
    fibreFieldName_(requiredFieldName(dict, "fibreField")),
    sheetFieldName_(requiredFieldName(dict, "sheetField")),
    sheetNormalFieldName_(requiredFieldName(dict, "sheetNormalField")),
    faceFibreFieldName_(requiredFieldName(dict, "faceFibreField")),
    faceSheetFieldName_(requiredFieldName(dict, "faceSheetField")),
    faceSheetNormalFieldName_
    (
        requiredFieldName(dict, "faceSheetNormalField")
    ),
    f0_
    (
        requiredReferenceFieldIOobject<volVectorField>
        (
            fibreFieldName_, mesh
        ),
        mesh
    ),
    s0_
    (
        requiredReferenceFieldIOobject<volVectorField>
        (
            sheetFieldName_, mesh
        ),
        mesh
    ),
    n0_
    (
        requiredReferenceFieldIOobject<volVectorField>
        (
            sheetNormalFieldName_, mesh
        ),
        mesh
    ),
    f0f_
    (
        requiredReferenceFieldIOobject<surfaceVectorField>
        (
            faceFibreFieldName_, mesh
        ),
        mesh
    ),
    s0f_
    (
        requiredReferenceFieldIOobject<surfaceVectorField>
        (
            faceSheetFieldName_, mesh
        ),
        mesh
    ),
    n0f_
    (
        requiredReferenceFieldIOobject<surfaceVectorField>
        (
            faceSheetNormalFieldName_, mesh
        ),
        mesh
    ),
    muRepresentative_("muRepresentative", dimPressure, 0.0)
{
    validateMaterialParameters();
    validateReferenceFields();
    calculateRepresentativeShearModulus();

    Info<< "ArosticaHolzapfelOgdenElastic options:" << nl
        << "    compressionSwitch = " << compressionSwitch_ << nl
        << "    compressionSwitchK = " << compressionSwitchK_ << nl
        << "    bulkModulus = " << bulkModulus_ << nl
        << "    impKcoeff = " << impKcoeff_ << nl
        << "    returns complete non-volumetric stress; no -pI or U(J)"
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ArosticaHolzapfelOgdenElastic::~ArosticaHolzapfelOgdenElastic()
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::ArosticaHolzapfelOgdenElastic::rho() const
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


Foam::tmp<Foam::volScalarField>
Foam::ArosticaHolzapfelOgdenElastic::impK() const
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
            impKcoeff_*muRepresentative_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::ArosticaHolzapfelOgdenElastic::bulkModulus() const
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
Foam::ArosticaHolzapfelOgdenElastic::shearModulus() const
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
            muRepresentative_
        )
    );
}


void Foam::ArosticaHolzapfelOgdenElastic::correct(volSymmTensorField& sigma)
{
    dimensionedScalar zeroBulkModulus("K", dimPressure, 0.0);

    if (updateF(sigma, muRepresentative_, zeroBulkModulus))
    {
        return;
    }

    const volTensorField& deformationGradient = F();

    forAll(sigma, cellI)
    {
        sigma[cellI] =
            evaluateConstitutive
            (
                deformationGradient[cellI],
                f0_[cellI],
                s0_[cellI],
                n0_[cellI]
            ).sigma;
    }

    forAll(sigma.boundaryField(), patchI)
    {
        symmTensorField& sigmaPatch = sigma.boundaryFieldRef()[patchI];
        const tensorField& deformationGradientPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& f0Patch = f0_.boundaryField()[patchI];
        const vectorField& s0Patch = s0_.boundaryField()[patchI];
        const vectorField& n0Patch = n0_.boundaryField()[patchI];

        forAll(sigmaPatch, faceI)
        {
            sigmaPatch[faceI] =
                evaluateConstitutive
                (
                    deformationGradientPatch[faceI],
                    f0Patch[faceI],
                    s0Patch[faceI],
                    n0Patch[faceI]
                ).sigma;
        }
    }
}


void Foam::ArosticaHolzapfelOgdenElastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    dimensionedScalar zeroBulkModulus("K", dimPressure, 0.0);

    if (updateF(sigma, muRepresentative_, zeroBulkModulus))
    {
        return;
    }

    const surfaceTensorField& deformationGradient = Ff();

    forAll(sigma, faceI)
    {
        sigma[faceI] =
            evaluateConstitutive
            (
                deformationGradient[faceI],
                f0f_[faceI],
                s0f_[faceI],
                n0f_[faceI]
            ).sigma;
    }

    forAll(sigma.boundaryField(), patchI)
    {
        symmTensorField& sigmaPatch = sigma.boundaryFieldRef()[patchI];
        const tensorField& deformationGradientPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& f0Patch = f0f_.boundaryField()[patchI];
        const vectorField& s0Patch = s0f_.boundaryField()[patchI];
        const vectorField& n0Patch = n0f_.boundaryField()[patchI];

        forAll(sigmaPatch, faceI)
        {
            sigmaPatch[faceI] =
                evaluateConstitutive
                (
                    deformationGradientPatch[faceI],
                    f0Patch[faceI],
                    s0Patch[faceI],
                    n0Patch[faceI]
                ).sigma;
        }
    }
}


void Foam::ArosticaHolzapfelOgdenElastic::materialTangentField
(
    List<mat66>& matTan
) const
{
    passiveMaterialTangentField(matTan);
}


void Foam::ArosticaHolzapfelOgdenElastic::passiveMaterialTangentField
(
    List<mat66>& matTan
) const
{
    const surfaceTensorField& deformationGradient = Ff();
    matTan.resize(mesh().nFaces());

    forAll(matTan, faceI)
    {
        matTan[faceI].clear();
    }

    const label nInternalFaces = mesh().nInternalFaces();
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        fillTangent
        (
            matTan[faceI],
            deformationGradient[faceI],
            f0f_[faceI],
            s0f_[faceI],
            n0f_[faceI]
        );
    }

    forAll(deformationGradient.boundaryField(), patchI)
    {
        const tensorField& deformationGradientPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& f0Patch = f0f_.boundaryField()[patchI];
        const vectorField& s0Patch = s0f_.boundaryField()[patchI];
        const vectorField& n0Patch = n0f_.boundaryField()[patchI];
        const label start = mesh().boundaryMesh()[patchI].start();

        forAll(deformationGradientPatch, faceI)
        {
            fillTangent
            (
                matTan[start + faceI],
                deformationGradientPatch[faceI],
                f0Patch[faceI],
                s0Patch[faceI],
                n0Patch[faceI]
            );
        }
    }
}


bool Foam::ArosticaHolzapfelOgdenElastic::passiveNominalTangentField
(
    List<tensor>& tangent
) const
{
    const surfaceTensorField& deformationGradient = Ff();
    tangent.setSize(9*mesh().nFaces(), tensor::zero);

    auto fillFace =
    [this, &tangent]
    (
        const label faceI,
        const tensor& F,
        const vector& f0,
        const vector& s0,
        const vector& n0
    )
    {
        for (label componentI = 0; componentI < 9; ++componentI)
        {
            tensor Fplus(F);
            tensor Fminus(F);
            Fplus[componentI] += tangentEps_;
            Fminus[componentI] -= tangentEps_;

            const tensor Pplus
            (
                Fplus & evaluateConstitutive(Fplus, f0, s0, n0).S
            );
            const tensor Pminus
            (
                Fminus & evaluateConstitutive(Fminus, f0, s0, n0).S
            );

            // Ppassive = J*sigma*F^-T = F*S for the passive
            // second-Piola stress returned by the constitutive kernel.
            tangent[9*faceI + componentI] =
                (Pplus - Pminus)/(2.0*tangentEps_);
        }
    };

    for (label faceI = 0; faceI < mesh().nInternalFaces(); ++faceI)
    {
        fillFace
        (
            faceI,
            deformationGradient[faceI],
            f0f_[faceI],
            s0f_[faceI],
            n0f_[faceI]
        );
    }

    forAll(deformationGradient.boundaryField(), patchI)
    {
        const tensorField& FPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& fPatch = f0f_.boundaryField()[patchI];
        const vectorField& sPatch = s0f_.boundaryField()[patchI];
        const vectorField& nPatch = n0f_.boundaryField()[patchI];
        const label start = mesh().boundaryMesh()[patchI].start();

        forAll(FPatch, faceI)
        {
            fillFace
            (
                start + faceI,
                FPatch[faceI],
                fPatch[faceI],
                sPatch[faceI],
                nPatch[faceI]
            );
        }
    }

    return true;
}


void Foam::ArosticaHolzapfelOgdenElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}


// ************************************************************************* //
