/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ArosticaHolzapfelOgdenViscoelastic.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ArosticaHolzapfelOgdenViscoelastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, ArosticaHolzapfelOgdenViscoelastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::ArosticaHolzapfelOgdenViscoelastic::fieldName
(
    const word& prefix
) const
{
    return "Arostica" + prefix + '_' + name();
}


template<class FieldType>
bool Foam::ArosticaHolzapfelOgdenViscoelastic::fieldPresent
(
    const word& fieldNameValue
) const
{
    IOobject io
    (
        fieldNameValue,
        mesh().time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

#ifdef OPENFOAM_ORG
    return io.typeHeaderOk<FieldType>(true);
#else
    return io.typeHeaderOk<FieldType>(true, false, false);
#endif
}


Foam::word Foam::ArosticaHolzapfelOgdenViscoelastic::readDdtScheme() const
{
    if (!mesh().schemesDict().found("ddtSchemes"))
    {
        FatalErrorInFunction
            << "The fvSchemes dictionary has no ddtSchemes sub-dictionary"
            << exit(FatalError);
    }

    const dictionary& ddtSchemes =
        mesh().schemesDict().subDict("ddtSchemes");

    if (!ddtSchemes.found("default"))
    {
        FatalErrorInFunction
            << "The ddtSchemes dictionary has no default entry"
            << exit(FatalError);
    }

    const word scheme(ddtSchemes.lookup("default"));

    if (scheme != "backward" && scheme != "Euler")
    {
        FatalErrorInFunction
            << "Unsupported Aróstica material viscosity ddt scheme "
            << scheme << ". Supported schemes are backward and Euler"
            << exit(FatalError);
    }

    return scheme;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::temporalCoefficients
(
    scalar& current,
    scalar& old,
    scalar& oldOld
) const
{
    const scalar deltaT = mesh().time().deltaTValue();

    if (!std::isfinite(deltaT) || deltaT <= SMALL)
    {
        FatalErrorInFunction
            << "Invalid deltaT for Aróstica material viscosity: " << deltaT
            << exit(FatalError);
    }

    if (ddtScheme_ == "Euler" || !haveSecondHistory_)
    {
        current = 1.0/deltaT;
        old = -current;
        oldOld = 0.0;
        return;
    }

    const scalar deltaT0 = mesh().time().deltaT0Value();

    if (!std::isfinite(deltaT0) || deltaT0 <= SMALL)
    {
        FatalErrorInFunction
            << "The backward Aróstica material derivative requires a "
            << "positive previous deltaT; got " << deltaT0
            << exit(FatalError);
    }

    const scalar coefft = 1.0 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    current = coefft/deltaT;
    old = -coefft0/deltaT;
    oldOld = coefft00/deltaT;
}


Foam::symmTensor Foam::ArosticaHolzapfelOgdenViscoelastic::strain
(
    const tensor& F
)
{
    return 0.5*(symm(F.T() & F) - symmTensor::I);
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::initialiseHistoryIfRequired()
{
    if (haveHistory_)
    {
        return;
    }

    // If no explicit history was written, the current E fields are the
    // accepted state read at startup (zero on a fresh undeformed case).
    EOld_ = E_;
    EOldOld_ = E_;
    EfOld_ = Ef_;
    EfOldOld_ = Ef_;
    haveHistory_ = true;
}


Foam::symmTensor
Foam::ArosticaHolzapfelOgdenViscoelastic::viscousStress
(
    const tensor& F,
    const symmTensor& EOld,
    const symmTensor& EOldOld,
    const scalar current,
    const scalar old,
    const scalar oldOld
) const
{
    const symmTensor Edot
    (
        current*strain(F) + old*EOld + oldOld*EOldOld
    );
    const symmTensor Sviscous(eta_.value()*Edot);
    const scalar J = det(F);

    if (!std::isfinite(J) || J <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid deformation Jacobian J = " << J
            << " in Aróstica viscous stress"
            << exit(FatalError);
    }

    return symm(F & Sviscous & F.T())/J;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::updateViscousFields()
{
    initialiseHistoryIfRequired();

    scalar current = 0.0;
    scalar old = 0.0;
    scalar oldOld = 0.0;
    temporalCoefficients(current, old, oldOld);

    const volTensorField& deformationGradient = F();
    forAll(E_, cellI)
    {
        E_[cellI] = strain(deformationGradient[cellI]);
        Edot_[cellI] =
            current*E_[cellI] + old*EOld_[cellI] + oldOld*EOldOld_[cellI];
        Sviscous_[cellI] = eta_.value()*Edot_[cellI];
        sigmaViscous_[cellI] =
            viscousStress
            (
                deformationGradient[cellI],
                EOld_[cellI],
                EOldOld_[cellI],
                current,
                old,
                oldOld
            );
    }

    forAll(E_.boundaryField(), patchI)
    {
        symmTensorField& EPatch = E_.boundaryFieldRef()[patchI];
        symmTensorField& EdotPatch = Edot_.boundaryFieldRef()[patchI];
        symmTensorField& SPatch = Sviscous_.boundaryFieldRef()[patchI];
        symmTensorField& sigmaPatch =
            sigmaViscous_.boundaryFieldRef()[patchI];
        const tensorField& FPatch =
            deformationGradient.boundaryField()[patchI];
        const symmTensorField& EOldPatch = EOld_.boundaryField()[patchI];
        const symmTensorField& EOldOldPatch =
            EOldOld_.boundaryField()[patchI];

        forAll(EPatch, faceI)
        {
            EPatch[faceI] = strain(FPatch[faceI]);
            EdotPatch[faceI] =
                current*EPatch[faceI]
              + old*EOldPatch[faceI]
              + oldOld*EOldOldPatch[faceI];
            SPatch[faceI] = eta_.value()*EdotPatch[faceI];
            sigmaPatch[faceI] =
                viscousStress
                (
                    FPatch[faceI],
                    EOldPatch[faceI],
                    EOldOldPatch[faceI],
                    current,
                    old,
                    oldOld
                );
        }
    }

    const surfaceTensorField& faceDeformationGradient = Ff();
    forAll(Ef_, faceI)
    {
        Ef_[faceI] = strain(faceDeformationGradient[faceI]);
        Edotf_[faceI] =
            current*Ef_[faceI]
          + old*EfOld_[faceI]
          + oldOld*EfOldOld_[faceI];
        Sviscousf_[faceI] = eta_.value()*Edotf_[faceI];
        sigmaViscousf_[faceI] =
            viscousStress
            (
                faceDeformationGradient[faceI],
                EfOld_[faceI],
                EfOldOld_[faceI],
                current,
                old,
                oldOld
            );
    }

    forAll(Ef_.boundaryField(), patchI)
    {
        symmTensorField& EPatch = Ef_.boundaryFieldRef()[patchI];
        symmTensorField& EdotPatch = Edotf_.boundaryFieldRef()[patchI];
        symmTensorField& SPatch = Sviscousf_.boundaryFieldRef()[patchI];
        symmTensorField& sigmaPatch =
            sigmaViscousf_.boundaryFieldRef()[patchI];
        const tensorField& FPatch =
            faceDeformationGradient.boundaryField()[patchI];
        const symmTensorField& EOldPatch = EfOld_.boundaryField()[patchI];
        const symmTensorField& EOldOldPatch =
            EfOldOld_.boundaryField()[patchI];

        forAll(EPatch, faceI)
        {
            EPatch[faceI] = strain(FPatch[faceI]);
            EdotPatch[faceI] =
                current*EPatch[faceI]
              + old*EOldPatch[faceI]
              + oldOld*EOldOldPatch[faceI];
            SPatch[faceI] = eta_.value()*EdotPatch[faceI];
            sigmaPatch[faceI] =
                viscousStress
                (
                    FPatch[faceI],
                    EOldPatch[faceI],
                    EOldOldPatch[faceI],
                    current,
                    old,
                    oldOld
                );
        }
    }
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::fillViscoelasticTangent
(
    mat66& tangent,
    const tensor& F,
    const vector& f0,
    const vector& s0,
    const vector& n0,
    const symmTensor& EOld,
    const symmTensor& EOldOld,
    const scalar current,
    const scalar old,
    const scalar oldOld
) const
{
    tangent.clear();

    const label tensorComponent[6] =
    {
        tensor::XX, tensor::YY, tensor::ZZ,
        tensor::XY, tensor::YZ, tensor::XZ
    };

    for (label componentI = 0; componentI < symmTensor::nComponents; ++componentI)
    {
        tensor Fplus(F);
        tensor Fminus(F);
        Fplus[tensorComponent[componentI]] += tangentEps_;
        Fminus[tensorComponent[componentI]] -= tangentEps_;

        const symmTensor sigmaPlus =
            evaluateConstitutive(Fplus, f0, s0, n0).sigma
          + viscousStress
            (
                Fplus, EOld, EOldOld, current, old, oldOld
            );
        const symmTensor sigmaMinus =
            evaluateConstitutive(Fminus, f0, s0, n0).sigma
          + viscousStress
            (
                Fminus, EOld, EOldOld, current, old, oldOld
            );
        const symmTensor derivative
        (
            (sigmaPlus - sigmaMinus)/(2.0*tangentEps_)
        );

        tangent(symmTensor::XX, componentI) = derivative.xx();
        tangent(symmTensor::YY, componentI) = derivative.yy();
        tangent(symmTensor::ZZ, componentI) = derivative.zz();
        tangent(symmTensor::XY, componentI) = derivative.xy();
        tangent(symmTensor::YZ, componentI) = derivative.yz();
        tangent(symmTensor::XZ, componentI) = derivative.xz();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ArosticaHolzapfelOgdenViscoelastic::ArosticaHolzapfelOgdenViscoelastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    ArosticaHolzapfelOgdenElastic(name, mesh, dict, nonLinGeom),
    eta_(dict.lookup("eta")),
    writeDiagnostics_(dict.lookupOrDefault<Switch>("writeDiagnostics", false)),
    ddtScheme_(readDdtScheme()),
    E_
    (
        IOobject
        (
            fieldName("E"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    EOld_
    (
        IOobject
        (
            fieldName("EOld"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    EOldOld_
    (
        IOobject
        (
            fieldName("EOldOld"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Ef_
    (
        IOobject
        (
            fieldName("Ef"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    EfOld_
    (
        IOobject
        (
            fieldName("EfOld"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    EfOldOld_
    (
        IOobject
        (
            fieldName("EfOldOld"), mesh.time().timeName(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Edot_
    (
        IOobject
        (
            fieldName("Edot"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless/dimTime, symmTensor::zero)
    ),
    Sviscous_
    (
        IOobject
        (
            fieldName("Sviscous"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    sigmaViscous_
    (
        IOobject
        (
            fieldName("sigmaViscous"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    Edotf_
    (
        IOobject
        (
            fieldName("Edotf"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless/dimTime, symmTensor::zero)
    ),
    Sviscousf_
    (
        IOobject
        (
            fieldName("Sviscousf"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    sigmaViscousf_
    (
        IOobject
        (
            fieldName("sigmaViscousf"), mesh.time().timeName(), mesh,
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    haveHistory_
    (
        fieldPresent<volSymmTensorField>(fieldName("EOld"))
     && fieldPresent<surfaceSymmTensorField>(fieldName("EfOld"))
    ),
    haveSecondHistory_
    (
        fieldPresent<volSymmTensorField>(fieldName("EOldOld"))
     && fieldPresent<surfaceSymmTensorField>(fieldName("EfOldOld"))
    ),
    lastHistoryUpdateTimeIndex_(-1)
{
    if
    (
        eta_.dimensions() != dimPressure*dimTime
     || !std::isfinite(eta_.value())
     || eta_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "eta must have dimensions pressure*time and be finite and "
            << "non-negative; got " << eta_ << exit(FatalError);
    }

    if (writeDiagnostics_)
    {
        Edot_.writeOpt() = IOobject::AUTO_WRITE;
        Sviscous_.writeOpt() = IOobject::AUTO_WRITE;
        sigmaViscous_.writeOpt() = IOobject::AUTO_WRITE;
        Edotf_.writeOpt() = IOobject::AUTO_WRITE;
        Sviscousf_.writeOpt() = IOobject::AUTO_WRITE;
        sigmaViscousf_.writeOpt() = IOobject::AUTO_WRITE;
    }

    Info<< "ArosticaHolzapfelOgdenViscoelastic options:" << nl
        << "    eta = " << eta_ << nl
        << "    ddt scheme = " << ddtScheme_ << nl
        << "    explicit accepted history = " << haveHistory_ << nl
        << "    second history = " << haveSecondHistory_ << nl
        << "    complete viscous stress retained; no dev projection"
        << endl;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ArosticaHolzapfelOgdenViscoelastic::~ArosticaHolzapfelOgdenViscoelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ArosticaHolzapfelOgdenViscoelastic::correct
(
    volSymmTensorField& sigma
)
{
    ArosticaHolzapfelOgdenElastic::correct(sigma);
    updateViscousFields();
    sigma += sigmaViscous_;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    ArosticaHolzapfelOgdenElastic::correct(sigma);
    updateViscousFields();
    sigma += sigmaViscousf_;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::materialTangentField
(
    List<mat66>& matTan
) const
{
    scalar current = 0.0;
    scalar old = 0.0;
    scalar oldOld = 0.0;
    temporalCoefficients(current, old, oldOld);

    const surfaceTensorField& deformationGradient = faceDeformationGradient();
    matTan.resize(mesh().nFaces());

    forAll(matTan, faceI)
    {
        matTan[faceI].clear();
    }

    for (label faceI = 0; faceI < mesh().nInternalFaces(); ++faceI)
    {
        fillViscoelasticTangent
        (
            matTan[faceI], deformationGradient[faceI],
            f0f_[faceI], s0f_[faceI], n0f_[faceI],
            EfOld_[faceI], EfOldOld_[faceI], current, old, oldOld
        );
    }

    forAll(deformationGradient.boundaryField(), patchI)
    {
        const tensorField& FPatch =
            deformationGradient.boundaryField()[patchI];
        const vectorField& fPatch = f0f_.boundaryField()[patchI];
        const vectorField& sPatch = s0f_.boundaryField()[patchI];
        const vectorField& nPatch = n0f_.boundaryField()[patchI];
        const symmTensorField& oldPatch = EfOld_.boundaryField()[patchI];
        const symmTensorField& oldOldPatch =
            EfOldOld_.boundaryField()[patchI];
        const label start = mesh().boundaryMesh()[patchI].start();

        forAll(FPatch, faceI)
        {
            fillViscoelasticTangent
            (
                matTan[start + faceI], FPatch[faceI],
                fPatch[faceI], sPatch[faceI], nPatch[faceI],
                oldPatch[faceI], oldOldPatch[faceI],
                current, old, oldOld
            );
        }
    }
}


bool Foam::ArosticaHolzapfelOgdenViscoelastic::viscousNominalTangentField
(
    List<tensor>& tangent
) const
{
    scalar current = 0.0;
    scalar old = 0.0;
    scalar oldOld = 0.0;
    temporalCoefficients(current, old, oldOld);

    const surfaceTensorField& deformationGradient = faceDeformationGradient();
    tangent.setSize(9*mesh().nFaces(), tensor::zero);

    forAll(deformationGradient, faceI)
    {
        const tensor& F = deformationGradient[faceI];
        const symmTensor E = strain(F);
        const symmTensor& EOld = EfOld_[faceI];
        const symmTensor& EOldOld = EfOldOld_[faceI];

        const symmTensor Sviscous
        (
            eta_.value()
           *(
                current*E
              + old*EOld
              + oldOld*EOldOld
            )
        );

        for (label componentI = 0; componentI < 9; ++componentI)
        {
            tensor dF(tensor::zero);
            dF[componentI] = 1.0;

            const symmTensor dE
            (
                0.5*symm
                (
                    (F.T() & dF) + (dF.T() & F)
                )
            );
            const symmTensor dSviscous
            (
                eta_.value()*current*dE
            );

            // Sviscous is symmetric, hence the nominal stress is exactly
            // Pviscous = J*sigmaViscous*F^-T = F*Sviscous.  This includes
            // both the current push-forward derivative and dEdot/dF.
            tangent[9*faceI + componentI] =
                (dF & Sviscous) + (F & dSviscous);
        }
    }

    return true;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::updateTotalFields()
{
    const label timeIndex = mesh().time().timeIndex();

    if (lastHistoryUpdateTimeIndex_ == timeIndex)
    {
        return;
    }

    initialiseHistoryIfRequired();

    EOldOld_ = EOld_;
    EOld_ = E_;
    EfOldOld_ = EfOld_;
    EfOld_ = Ef_;
    haveSecondHistory_ = true;
    lastHistoryUpdateTimeIndex_ = timeIndex;
}


void Foam::ArosticaHolzapfelOgdenViscoelastic::setRestart()
{
    ArosticaHolzapfelOgdenElastic::setRestart();
    E_.writeOpt() = IOobject::AUTO_WRITE;
    EOld_.writeOpt() = IOobject::AUTO_WRITE;
    EOldOld_.writeOpt() = IOobject::AUTO_WRITE;
    Ef_.writeOpt() = IOobject::AUTO_WRITE;
    EfOld_.writeOpt() = IOobject::AUTO_WRITE;
    EfOldOld_.writeOpt() = IOobject::AUTO_WRITE;
}


// ************************************************************************* //
