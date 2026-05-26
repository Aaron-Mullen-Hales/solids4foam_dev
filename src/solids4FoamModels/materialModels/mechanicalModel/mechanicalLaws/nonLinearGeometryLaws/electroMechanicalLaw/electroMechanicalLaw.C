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

#include "electroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, electroMechanicalLaw, nonLinGeomMechLaw
    );
}

namespace
{

const Foam::dictionary& electroMechanicalFibreDict
(
    const Foam::dictionary& dict
)
{
    if
    (
        dict.lookupOrDefault<Foam::Switch>
        (
            "uniformFibreField",
            Foam::Switch(false)
        )
    )
    {
        return dict;
    }

    if (dict.found("passiveMechanicalLaw"))
    {
        const Foam::dictionary& passiveDict =
            dict.subDict("passiveMechanicalLaw");

        if
        (
            passiveDict.lookupOrDefault<Foam::Switch>
            (
                "uniformFibreField",
                Foam::Switch(false)
            )
        )
        {
            return passiveDict;
        }
    }

    return dict;
}


Foam::IOobject findElectroMechanicalFibreFieldIOobject
(
    const Foam::word& fieldName,
    const Foam::fvMesh& mesh
)
{
    Foam::IOobject io
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::READ_IF_PRESENT,
        Foam::IOobject::NO_WRITE
    );

#ifdef FOAMEXTEND
    if (!io.headerOk())
#elif defined(OPENFOAM_ORG)
    if (!io.typeHeaderOk<Foam::volVectorField>(true))
#else
    if (!io.typeHeaderOk<Foam::volVectorField>(true, false, false))
#endif
    {
        io.instance() = "0";

#ifdef FOAMEXTEND
        if (!io.headerOk())
#elif defined(OPENFOAM_ORG)
        if (!io.typeHeaderOk<Foam::volVectorField>(true))
#else
        if (!io.typeHeaderOk<Foam::volVectorField>(true, false, false))
#endif
        {
            FatalErrorInFunction
                << "Cannot find required fibre field " << fieldName
                << " in either " << mesh.time().timeName() << " or 0"
                << Foam::exit(Foam::FatalError);
        }
    }

    io.readOpt() = Foam::IOobject::MUST_READ;

    return io;
}


Foam::tmp<Foam::volVectorField> makeElectroMechanicalF0
(
    const Foam::fvMesh& mesh,
    const Foam::dictionary& dict
)
{
    const Foam::dictionary& fibreDict = electroMechanicalFibreDict(dict);

    if
    (
        fibreDict.lookupOrDefault<Foam::Switch>
        (
            "uniformFibreField",
            Foam::Switch(false)
        )
    )
    {
        return Foam::tmp<Foam::volVectorField>
        (
            new Foam::volVectorField
            (
                Foam::IOobject
                (
                    "f0",
                    mesh.time().timeName(),
                    mesh,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::NO_WRITE
                ),
                mesh,
                Foam::dimensionedVector
                (
                    "f0",
                    Foam::dimless,
                    fibreDict.lookup("f0")
                )
            )
        );
    }

    return Foam::tmp<Foam::volVectorField>
    (
        new Foam::volVectorField
        (
            findElectroMechanicalFibreFieldIOobject("f0", mesh),
            mesh
        )
    );
}


Foam::tmp<Foam::surfaceVectorField> makeElectroMechanicalF0f
(
    const Foam::fvMesh& mesh,
    const Foam::dictionary& dict
)
{
    const Foam::dictionary& fibreDict = electroMechanicalFibreDict(dict);

    if
    (
        fibreDict.lookupOrDefault<Foam::Switch>
        (
            "uniformFibreField",
            Foam::Switch(false)
        )
    )
    {
        return Foam::tmp<Foam::surfaceVectorField>
        (
            new Foam::surfaceVectorField
            (
                Foam::IOobject
                (
                    "f0f",
                    mesh.time().timeName(),
                    mesh,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::NO_WRITE
                ),
                mesh,
                Foam::dimensionedVector
                (
                    "f0",
                    Foam::dimless,
                    fibreDict.lookup("f0")
                )
            )
        );
    }

    return Foam::tmp<Foam::surfaceVectorField>
    (
        new Foam::surfaceVectorField
        (
            Foam::IOobject
            (
                "f0f",
                mesh.time().timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE
            ),
            Foam::fvc::interpolate
            (
                Foam::volVectorField
                (
                    findElectroMechanicalFibreFieldIOobject("f0", mesh),
                    mesh
                )
            )
        )
    );
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::electroMechanicalLaw::electroMechanicalLaw
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
    f0_(makeElectroMechanicalF0(mesh, dict)),
    f0f0_("f0f0", sqr(f0_)),
    f0f_(makeElectroMechanicalF0f(mesh, dict)),
    f0f0f_("f0f0f", sqr(f0f_)),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    useFieldTa_(false),
    fieldTaChecked_(false)
{
    if (rampTime_ < 0.0)
    {
        FatalErrorIn("electroMechanicalLaw::electroMechanicalLaw(...)")
            << "rampTime should be greater than or equal to zero"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalLaw::~electroMechanicalLaw()
{}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::impK() const
{
    return passiveMechLawPtr_->impK();
}


void Foam::electroMechanicalLaw::materialTangentField(List<mat66>& matTan) const
{
    passiveMechLawPtr_->materialTangentField(matTan);
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::bulkModulus() const
{
    return passiveMechLawPtr_->bulkModulus();
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::shearModulus() const
{
    return passiveMechLawPtr_->shearModulus();
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Lazy check for field-based active tension
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the deformation gradient
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    if (useFieldTa_)
    {
        // Field-based active tension from the coupling model
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
        sigma += symm(F & (Ta*f0f0_) & F.T())/J;
    }
    else
    {
        // Constant active tension with optional ramp
        dimensionedScalar currentTa = Ta_;
        if (mesh().time().value() < rampTime_)
        {
            currentTa = (mesh().time().value()/rampTime_)*Ta_;
        }

        sigma += symm(F & (currentTa*f0f0_) & F.T())/J;
    }
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Lazy check for field-based active tension (same as vol variant)
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the deformation gradient
    const surfaceTensorField& F = this->Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    if (useFieldTa_)
    {
        // Interpolate field-based active tension to faces
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        const surfaceScalarField Taf(fvc::interpolate(Ta));

        // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
        sigma += J*symm(F & (Taf*f0f0f_) & F.T());
    }
    else
    {
        // Constant active tension with optional ramp
        dimensionedScalar currentTa = Ta_;
        if (mesh().time().value() < rampTime_)
        {
            currentTa = (mesh().time().value()/rampTime_)*Ta_;
        }

        sigma += J*symm(F & (currentTa*f0f0f_) & F.T());
    }
}


// ************************************************************************* //
