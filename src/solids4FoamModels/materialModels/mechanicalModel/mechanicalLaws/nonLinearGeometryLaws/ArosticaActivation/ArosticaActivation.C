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

#include "ArosticaActivation.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ArosticaActivation, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, ArosticaActivation, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ArosticaActivation::validateActivationParameters() const
{
    if
    (
        sigma0_.dimensions() != dimPressure
     || !std::isfinite(sigma0_.value())
     || sigma0_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "sigma0 must have dimensions of pressure and be finite and "
            << "non-negative; got " << sigma0_ << exit(FatalError);
    }

    if
    (
        gamma_.dimensions() != dimTime
     || !std::isfinite(gamma_.value())
     || gamma_.value() <= 0.0
    )
    {
        FatalErrorInFunction
            << "gamma must have dimensions of time and be finite and "
            << "positive; got " << gamma_ << exit(FatalError);
    }

    const dimensionSet rateDims(dimless/dimTime);

    if
    (
        alphaMin_.dimensions() != rateDims
     || !std::isfinite(alphaMin_.value())
    )
    {
        FatalErrorInFunction
            << "alphaMin must have dimensions of inverse time and be "
            << "finite; got " << alphaMin_ << exit(FatalError);
    }

    if
    (
        alphaMax_.dimensions() != rateDims
     || !std::isfinite(alphaMax_.value())
    )
    {
        FatalErrorInFunction
            << "alphaMax must have dimensions of inverse time and be "
            << "finite; got " << alphaMax_ << exit(FatalError);
    }

    if
    (
        tSys_.dimensions() != dimTime
     || !std::isfinite(tSys_.value())
     || tSys_.value() < 0.0
    )
    {
        FatalErrorInFunction
            << "tSys must have dimensions of time and be finite and "
            << "non-negative; got " << tSys_ << exit(FatalError);
    }

    if
    (
        tDias_.dimensions() != dimTime
     || !std::isfinite(tDias_.value())
     || tDias_.value() <= tSys_.value()
    )
    {
        FatalErrorInFunction
            << "tDias must have dimensions of time, be finite and be "
            << "greater than tSys; got tDias = " << tDias_
            << ", tSys = " << tSys_ << exit(FatalError);
    }
}


Foam::scalar Foam::ArosticaActivation::sPlus(const scalar deltaT) const
{
    return 0.5*(1.0 + std::tanh(deltaT/gamma_.value()));
}


Foam::scalar Foam::ArosticaActivation::sMinus(const scalar deltaT) const
{
    return 0.5*(1.0 - std::tanh(deltaT/gamma_.value()));
}


Foam::scalar Foam::ArosticaActivation::activationFunction
(
    const scalar t
) const
{
    // Eq. (6): f(t) = S+(t - tSys) S-(t - tDias)
    const scalar f = sPlus(t - tSys_.value())*sMinus(t - tDias_.value());

    return alphaMax_.value()*f + alphaMin_.value()*(1.0 - f);
}


Foam::scalar Foam::ArosticaActivation::tauDerivative
(
    const scalar t,
    const scalar tau
) const
{
    // Eq. (5): dtau/dt = -|a(t)| tau + sigma0 |a(t)|_+
    const scalar a = activationFunction(t);

    return -mag(a)*tau + sigma0_.value()*max(a, 0.0);
}


Foam::scalar Foam::ArosticaActivation::integrateTau(const scalar tEnd) const
{
    if (tEnd <= 0.0)
    {
        return 0.0;
    }

    // tau(t) is spatially uniform and depends only on time, so it is
    // integrated here with a fixed sub-step fine enough to resolve the
    // gamma_ transition time scale, independent of the outer solver time
    // step
    const scalar subStep = min(gamma_.value()/20.0, 1.0e-4);

    const label nSteps = max(label(tEnd/subStep) + 1, 1);

    if (nSteps > 20000000)
    {
        FatalErrorInFunction
            << "Excessive number of sub-steps (" << nSteps << ") required "
            << "to integrate the Aróstica active tension to t = " << tEnd
            << "; check the gamma parameter" << exit(FatalError);
    }

    const scalar h = tEnd/nSteps;

    scalar t = 0.0;
    scalar tau = 0.0;

    for (label stepI = 0; stepI < nSteps; ++stepI)
    {
        const scalar k1 = tauDerivative(t, tau);
        const scalar k2 = tauDerivative(t + 0.5*h, tau + 0.5*h*k1);
        const scalar k3 = tauDerivative(t + 0.5*h, tau + 0.5*h*k2);
        const scalar k4 = tauDerivative(t + h, tau + h*k3);

        tau += (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        t += h;

        if (!std::isfinite(tau))
        {
            FatalErrorInFunction
                << "Non-finite active tension tau = " << tau
                << " encountered while integrating the Aróstica activation "
                << "model at t = " << t << exit(FatalError);
        }
    }

    return tau;
}


Foam::scalar Foam::ArosticaActivation::currentTau() const
{
    const label timeIndex = mesh().time().timeIndex();

    if (timeIndex != cachedTimeIndex_)
    {
        tau_ = integrateTau(mesh().time().value());
        cachedTimeIndex_ = timeIndex;

        Ta_ = dimensionedScalar("Ta", sigma0_.dimensions(), tau_);
    }

    return tau_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ArosticaActivation::ArosticaActivation
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    passiveMechLawPtr_
    (
        mechanicalLaw::NewNonLinGeomMechLaw
        (
            word(dict.subDict("passiveMechanicalLaw").lookup("type")),
            mesh,
            dict.subDict("passiveMechanicalLaw"),
            nonLinGeom
        )
    ),
    f0_
    (
        IOobject
        (
            "f0",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    f0f0_("f0f0", sqr(f0_)),
    f0f_
    (
        IOobject
        (
            "f0f",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    f0f0f_("f0f0f", sqr(f0f_)),
    sigma0_(dict.lookup("sigma0")),
    gamma_(dict.lookup("gamma")),
    alphaMin_(dict.lookup("alphaMin")),
    alphaMax_(dict.lookup("alphaMax")),
    tSys_(dict.lookup("tSys")),
    tDias_(dict.lookup("tDias")),
    Ta_
    (
        IOobject
        (
            "Ta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ta", sigma0_.dimensions(), 0.0)
    ),
    cachedTimeIndex_(-1),
    tau_(0.0)
{
    validateActivationParameters();

    Info<< "ArosticaActivation active stress model (Bestel-Clement-Sorine):"
        << nl
        << "    sigma0 = " << sigma0_ << nl
        << "    gamma = " << gamma_ << nl
        << "    alphaMin = " << alphaMin_ << nl
        << "    alphaMax = " << alphaMax_ << nl
        << "    tSys = " << tSys_ << nl
        << "    tDias = " << tDias_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ArosticaActivation::~ArosticaActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ArosticaActivation::impK() const
{
    return passiveMechLawPtr_->impK();
}


void Foam::ArosticaActivation::materialTangentField(List<mat66>& matTan) const
{
    passiveMechLawPtr_->materialTangentField(matTan);
}


Foam::tmp<Foam::volScalarField> Foam::ArosticaActivation::bulkModulus() const
{
    return passiveMechLawPtr_->bulkModulus();
}


Foam::tmp<Foam::volScalarField> Foam::ArosticaActivation::shearModulus() const
{
    return passiveMechLawPtr_->shearModulus();
}


bool Foam::ArosticaActivation::hasActiveStress() const
{
    return mag(currentTau()) > SMALL;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::ArosticaActivation::activeCauchyStress() const
{
    // Take a reference to the deformation gradient maintained by the
    // passive law; see the equivalent note in electroMechanicalLaw. The
    // passive law is bound as a const reference so that the public F()
    // const accessor is selected.
    const mechanicalLaw& passiveLaw = passiveMechLawPtr_();
    const volTensorField& F = passiveLaw.F();

    const volScalarField J(det(F));

    const dimensionedScalar Ta("Ta", sigma0_.dimensions(), currentTau());

    // Convert active 2nd Piola-Kirchhoff stress tau(t) f0 x f0 to Cauchy
    // stress
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaActive",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            symm(F & (Ta*f0f0_) & F.T())/J
        )
    );
}


Foam::tmp<Foam::surfaceSymmTensorField>
Foam::ArosticaActivation::activeCauchyStressf() const
{
    // See the note in activeCauchyStress()
    const mechanicalLaw& passiveLaw = passiveMechLawPtr_();
    const surfaceTensorField& F = passiveLaw.Ff();

    const surfaceScalarField J(det(F));

    const dimensionedScalar Ta("Ta", sigma0_.dimensions(), currentTau());

    // tau(t) is spatially uniform, so no interpolation from cells to faces
    // is required
    return tmp<surfaceSymmTensorField>
    (
        new surfaceSymmTensorField
        (
            IOobject
            (
                "sigmaActivef",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (1.0/J)*symm(F & (Ta*f0f0f_) & F.T())
        )
    );
}


void Foam::ArosticaActivation::correct(volSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if (hasActiveStress())
    {
        const tmp<volSymmTensorField> tSigmaActive = activeCauchyStress();
        sigma += tSigmaActive();
    }
}


void Foam::ArosticaActivation::correctStressComponents
(
    volSymmTensorField& sigmaToProject,
    volSymmTensorField& sigmaPreserved
)
{
    passiveMechLawPtr_->correctStressComponents
    (
        sigmaToProject, sigmaPreserved
    );

    if (hasActiveStress())
    {
        const tmp<volSymmTensorField> tSigmaActive = activeCauchyStress();
        sigmaPreserved += tSigmaActive();
    }
}


void Foam::ArosticaActivation::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if (hasActiveStress())
    {
        const tmp<surfaceSymmTensorField> tSigmaActive = activeCauchyStressf();
        sigma += tSigmaActive();
    }
}


void Foam::ArosticaActivation::correctStressComponents
(
    surfaceSymmTensorField& sigmaToProject,
    surfaceSymmTensorField& sigmaPreserved
)
{
    passiveMechLawPtr_->correctStressComponents
    (
        sigmaToProject, sigmaPreserved
    );

    if (hasActiveStress())
    {
        const tmp<surfaceSymmTensorField> tSigmaActive =
            activeCauchyStressf();
        sigmaPreserved += tSigmaActive();
    }
}


void Foam::ArosticaActivation::updateFf()
{
    passiveMechLawPtr_->updateFf();
}


void Foam::ArosticaActivation::updateTotalFields()
{
    passiveMechLawPtr_->updateTotalFields();
}


void Foam::ArosticaActivation::setRestart()
{
    passiveMechLawPtr_->setRestart();
}


// ************************************************************************* //
