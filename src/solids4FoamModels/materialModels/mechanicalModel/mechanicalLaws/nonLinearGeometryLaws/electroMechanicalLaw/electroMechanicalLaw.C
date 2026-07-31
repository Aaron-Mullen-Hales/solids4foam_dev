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
#include "../GuccioneElastic/GuccioneElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "ZoneIDs.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, electroMechanicalLaw, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace
{

using namespace Foam;

bool finiteVector(const Foam::vector& v)
{
    return
        std::isfinite(v.x())
     && std::isfinite(v.y())
     && std::isfinite(v.z());
}


void setBoundaryMaskToOwner(Foam::volScalarField& mask)
{
    const Foam::fvMesh& mesh = mask.mesh();

    forAll(mask.boundaryField(), patchI)
    {
        Foam::scalarField& patch = mask.boundaryFieldRef()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            patch[faceI] = mask[faceCells[faceI]];
        }
    }
}


void setBoundaryMaskToOwner
(
    Foam::surfaceScalarField& faceMask,
    const Foam::volScalarField& cellMask
)
{
    const Foam::fvMesh& mesh = cellMask.mesh();

    forAll(faceMask.boundaryField(), patchI)
    {
        Foam::scalarField& patch = faceMask.boundaryFieldRef()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            patch[faceI] = cellMask[faceCells[faceI]];
        }
    }
}


void validateBoundaryMaskOwnerValues
(
    const Foam::volScalarField& mask,
    const Foam::word& name
)
{
    const Foam::fvMesh& mesh = mask.mesh();

    forAll(mask.boundaryField(), patchI)
    {
        const Foam::scalarField& patch = mask.boundaryField()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            const Foam::scalar ownerValue = mask[faceCells[faceI]];

            if (Foam::mag(patch[faceI] - ownerValue) > Foam::SMALL)
            {
                FatalErrorInFunction
                    << name << " boundary value does not match owner cell"
                    << nl
                    << "    patch = "
                    << mesh.boundary()[patchI].name() << nl
                    << "    face = " << faceI << nl
                    << "    owner cell = " << faceCells[faceI] << nl
                    << "    boundary value = " << patch[faceI] << nl
                    << "    owner value = " << ownerValue
                    << abort(Foam::FatalError);
            }
        }
    }
}


void validateBoundaryMaskOwnerValues
(
    const Foam::surfaceScalarField& faceMask,
    const Foam::volScalarField& cellMask,
    const Foam::word& name
)
{
    const Foam::fvMesh& mesh = cellMask.mesh();

    forAll(faceMask.boundaryField(), patchI)
    {
        const Foam::scalarField& patch = faceMask.boundaryField()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            const Foam::scalar ownerValue = cellMask[faceCells[faceI]];

            if (Foam::mag(patch[faceI] - ownerValue) > Foam::SMALL)
            {
                FatalErrorInFunction
                    << name << " boundary value does not match owner cell"
                    << nl
                    << "    patch = "
                    << mesh.boundary()[patchI].name() << nl
                    << "    face = " << faceI << nl
                    << "    owner cell = " << faceCells[faceI] << nl
                    << "    boundary value = " << patch[faceI] << nl
                    << "    owner value = " << ownerValue
                    << abort(Foam::FatalError);
            }
        }
    }
}


void reportBoundaryMaskSamples
(
    const Foam::volScalarField& mask,
    const Foam::word& name
)
{
    const Foam::fvMesh& mesh = mask.mesh();
    bool foundSelected = false;
    bool foundUnselected = false;

    forAll(mask.boundaryField(), patchI)
    {
        const Foam::scalarField& patch = mask.boundaryField()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            const Foam::scalar ownerValue = mask[faceCells[faceI]];

            if (!foundSelected && ownerValue > 0.5)
            {
                Info<< "    " << name << " selected boundary sample: patch "
                    << mesh.boundary()[patchI].name()
                    << ", face " << faceI
                    << ", owner cell " << faceCells[faceI]
                    << ", owner/boundary = "
                    << ownerValue << " / " << patch[faceI] << nl;
                foundSelected = true;
            }
            else if (!foundUnselected && ownerValue < 0.5)
            {
                Info<< "    " << name
                    << " unselected boundary sample: patch "
                    << mesh.boundary()[patchI].name()
                    << ", face " << faceI
                    << ", owner cell " << faceCells[faceI]
                    << ", owner/boundary = "
                    << ownerValue << " / " << patch[faceI] << nl;
                foundUnselected = true;
            }

            if (foundSelected && foundUnselected)
            {
                return;
            }
        }
    }
}


void reportBoundaryMaskSamples
(
    const Foam::surfaceScalarField& faceMask,
    const Foam::volScalarField& cellMask,
    const Foam::word& name
)
{
    const Foam::fvMesh& mesh = cellMask.mesh();
    bool foundSelected = false;
    bool foundUnselected = false;

    forAll(faceMask.boundaryField(), patchI)
    {
        const Foam::scalarField& patch = faceMask.boundaryField()[patchI];
        const Foam::labelUList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(patch, faceI)
        {
            const Foam::scalar ownerValue = cellMask[faceCells[faceI]];

            if (!foundSelected && ownerValue > 0.5)
            {
                Info<< "    " << name << " selected boundary sample: patch "
                    << mesh.boundary()[patchI].name()
                    << ", face " << faceI
                    << ", owner cell " << faceCells[faceI]
                    << ", owner/face = "
                    << ownerValue << " / " << patch[faceI] << nl;
                foundSelected = true;
            }
            else if (!foundUnselected && ownerValue < 0.5)
            {
                Info<< "    " << name
                    << " unselected boundary sample: patch "
                    << mesh.boundary()[patchI].name()
                    << ", face " << faceI
                    << ", owner cell " << faceCells[faceI]
                    << ", owner/face = "
                    << ownerValue << " / " << patch[faceI] << nl;
                foundUnselected = true;
            }

            if (foundSelected && foundUnselected)
            {
                return;
            }
        }
    }
}


struct MagnitudeStats
{
    Foam::label n;
    Foam::scalar minBefore;
    Foam::scalar maxBefore;
    Foam::scalar minAfter;
    Foam::scalar maxAfter;

    MagnitudeStats()
    :
        n(0),
        minBefore(Foam::VGREAT),
        maxBefore(-Foam::VGREAT),
        minAfter(Foam::VGREAT),
        maxAfter(-Foam::VGREAT)
    {}

    void addBefore(const Foam::scalar value)
    {
        ++n;
        minBefore = Foam::min(minBefore, value);
        maxBefore = Foam::max(maxBefore, value);
    }

    void addAfter(const Foam::scalar value)
    {
        minAfter = Foam::min(minAfter, value);
        maxAfter = Foam::max(maxAfter, value);
    }
};


void printMagnitudeStats
(
    const Foam::word& name,
    const MagnitudeStats& stats
)
{
    Info<< "    " << name
        << " magnitude before normalisation min/max = "
        << stats.minBefore << " / " << stats.maxBefore << nl
        << "    " << name
        << " magnitude after normalisation min/max = "
        << stats.minAfter << " / " << stats.maxAfter << nl;
}


void normaliseVectorField
(
    Foam::volVectorField& field,
    const Foam::word& name,
    const Foam::scalar tolerance,
    MagnitudeStats& stats
)
{
#ifdef OPENFOAM_NOT_EXTEND
    Foam::vectorField& internal = field.primitiveFieldRef();
#else
    Foam::vectorField& internal = field.internalField();
#endif

    forAll(internal, cellI)
    {
        const Foam::scalar m = Foam::mag(internal[cellI]);
        stats.addBefore(m);
        if (!finiteVector(internal[cellI]) || m <= tolerance)
        {
            FatalErrorInFunction
                << "Invalid " << name << " cell fibre" << nl
                << "    cell = " << cellI << nl
                << "    value = " << internal[cellI] << nl
                << "    magnitude = " << m << abort(Foam::FatalError);
        }
        internal[cellI] /= m;
        stats.addAfter(Foam::mag(internal[cellI]));
    }

    forAll(field.boundaryField(), patchI)
    {
        Foam::vectorField& patch = field.boundaryFieldRef()[patchI];
        forAll(patch, faceI)
        {
            const Foam::scalar m = Foam::mag(patch[faceI]);
            stats.addBefore(m);
            if (!finiteVector(patch[faceI]) || m <= tolerance)
            {
                FatalErrorInFunction
                    << "Invalid " << name << " boundary fibre" << nl
                    << "    patch = " << field.mesh().boundary()[patchI].name() << nl
                    << "    face = " << faceI << nl
                    << "    value = " << patch[faceI] << nl
                    << "    magnitude = " << m << abort(Foam::FatalError);
            }
            patch[faceI] /= m;
            stats.addAfter(Foam::mag(patch[faceI]));
        }
    }
}


void normaliseVectorField
(
    Foam::surfaceVectorField& field,
    const Foam::word& name,
    const Foam::scalar tolerance,
    MagnitudeStats& stats
)
{
#ifdef OPENFOAM_NOT_EXTEND
    Foam::vectorField& internal = field.primitiveFieldRef();
#else
    Foam::vectorField& internal = field.internalField();
#endif

    forAll(internal, faceI)
    {
        const Foam::scalar m = Foam::mag(internal[faceI]);
        stats.addBefore(m);
        if (!finiteVector(internal[faceI]) || m <= tolerance)
        {
            FatalErrorInFunction
                << "Invalid " << name << " internal face fibre" << nl
                << "    face = " << faceI << nl
                << "    value = " << internal[faceI] << nl
                << "    magnitude = " << m << abort(Foam::FatalError);
        }
        internal[faceI] /= m;
        stats.addAfter(Foam::mag(internal[faceI]));
    }

    forAll(field.boundaryField(), patchI)
    {
        Foam::vectorField& patch = field.boundaryFieldRef()[patchI];
        forAll(patch, faceI)
        {
            const Foam::scalar m = Foam::mag(patch[faceI]);
            stats.addBefore(m);
            if (!finiteVector(patch[faceI]) || m <= tolerance)
            {
                FatalErrorInFunction
                    << "Invalid " << name << " boundary face fibre" << nl
                    << "    patch = " << field.mesh().boundary()[patchI].name() << nl
                    << "    face = " << faceI << nl
                    << "    value = " << patch[faceI] << nl
                    << "    magnitude = " << m << abort(Foam::FatalError);
            }
            patch[faceI] /= m;
            stats.addAfter(Foam::mag(patch[faceI]));
        }
    }
}


struct DifferenceStats
{
    Foam::label n;
    Foam::scalar maxValue;
    Foam::scalar sumSqr;

    DifferenceStats()
    :
        n(0),
        maxValue(0),
        sumSqr(0)
    {}

    void add(const Foam::scalar value)
    {
        ++n;
        maxValue = Foam::max(maxValue, value);
        sumSqr += Foam::sqr(value);
    }

    Foam::scalar rms() const
    {
        return n ? Foam::sqrt(sumSqr/Foam::scalar(n)) : 0;
    }
};


DifferenceStats compareDyads
(
    const Foam::volSymmTensorField& a,
    const Foam::volSymmTensorField& b
)
{
    DifferenceStats stats;

#ifdef OPENFOAM_NOT_EXTEND
    const Foam::symmTensorField& aI = a.primitiveField();
    const Foam::symmTensorField& bI = b.primitiveField();
#else
    const Foam::symmTensorField& aI = a.internalField();
    const Foam::symmTensorField& bI = b.internalField();
#endif

    forAll(aI, cellI)
    {
        stats.add(Foam::mag(aI[cellI] - bI[cellI]));
    }

    return stats;
}


DifferenceStats compareDyads
(
    const Foam::surfaceSymmTensorField& a,
    const Foam::surfaceSymmTensorField& b
)
{
    DifferenceStats stats;

#ifdef OPENFOAM_NOT_EXTEND
    const Foam::symmTensorField& aI = a.primitiveField();
    const Foam::symmTensorField& bI = b.primitiveField();
#else
    const Foam::symmTensorField& aI = a.internalField();
    const Foam::symmTensorField& bI = b.internalField();
#endif

    forAll(aI, faceI)
    {
        stats.add(Foam::mag(aI[faceI] - bI[faceI]));
    }

    forAll(a.boundaryField(), patchI)
    {
        const Foam::symmTensorField& aPatch = a.boundaryField()[patchI];
        const Foam::symmTensorField& bPatch = b.boundaryField()[patchI];

        forAll(aPatch, faceI)
        {
            stats.add(Foam::mag(aPatch[faceI] - bPatch[faceI]));
        }
    }

    return stats;
}


bool hasApexRegularisationDict(const Foam::dictionary& dict)
{
    return dict.found("apexRegularisation");
}


const Foam::dictionary& apexRegularisationDict
(
    const Foam::dictionary& dict
)
{
    return dict.subDict("apexRegularisation");
}


Foam::Switch readApexRegularisationSwitch
(
    const Foam::dictionary& dict,
    const Foam::word& keyword,
    const bool defaultValue
)
{
    if (!hasApexRegularisationDict(dict))
    {
        return Foam::Switch(defaultValue);
    }

    return
        apexRegularisationDict(dict).lookupOrDefault<Foam::Switch>
        (
            keyword,
            Foam::Switch(defaultValue)
        );
}


Foam::word readApexRegularisationZoneName(const Foam::dictionary& dict)
{
    if (!hasApexRegularisationDict(dict))
    {
        return Foam::word::null;
    }

    return
        apexRegularisationDict(dict).lookupOrDefault<Foam::word>
        (
            "cellZone",
            Foam::word::null
        );
}


Foam::word passiveMechanicalLawType(const Foam::dictionary& dict)
{
    return Foam::word(dict.subDict("passiveMechanicalLaw").lookup("type"));
}


Foam::dictionary passiveMechanicalLawDict(const Foam::dictionary& dict)
{
    Foam::dictionary passiveDict(dict.subDict("passiveMechanicalLaw"));

    if
    (
        hasApexRegularisationDict(dict)
     && !passiveDict.found("apexRegularisation")
    )
    {
        passiveDict.add
        (
            "apexRegularisation",
            apexRegularisationDict(dict)
        );
    }

    return passiveDict;
}

}

void Foam::electroMechanicalLaw::checkFieldTa() const
{
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else if (mag(Ta_.value()) > SMALL)
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }
}


Foam::dimensionedScalar Foam::electroMechanicalLaw::activeTension() const
{
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    return currentTa;
}


void Foam::electroMechanicalLaw::makeApexRegularisationMasks() const
{
    if (apexRegularisationMaskPtr_.valid())
    {
        return;
    }

    if (!apexRegularisationEnabled_)
    {
        return;
    }

    if
    (
        apexRegularisationZoneName_ != "apexLayer0"
     && apexRegularisationZoneName_ != "apexLayers0to1"
    )
    {
        FatalErrorInFunction
            << "Invalid apexRegularisation cellZone "
            << apexRegularisationZoneName_ << nl
            << "Allowed zones are apexLayer0 and apexLayers0to1"
            << abort(FatalError);
    }

    const cellZoneID zoneID
    (
        apexRegularisationZoneName_,
        mesh().cellZones()
    );

    if (!zoneID.active())
    {
        FatalErrorInFunction
            << "apexRegularisation cellZone not found: "
            << apexRegularisationZoneName_ << abort(FatalError);
    }

    const labelList& zoneCells = mesh().cellZones()[zoneID.index()];
    const label globalZoneSize =
        returnReduce(zoneCells.size(), sumOp<label>());

    if (globalZoneSize == 0)
    {
        FatalErrorInFunction
            << "apexRegularisation cellZone is empty: "
            << apexRegularisationZoneName_ << abort(FatalError);
    }

    apexRegularisationMaskPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "apexRegularisationMask",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                writeApexDiagnosticFields_
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& mask = apexRegularisationMaskPtr_();

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& maskI = mask.primitiveFieldRef();
#else
    scalarField& maskI = mask.internalField();
#endif

    label selectedCells = 0;
    scalar selectedVolume = 0.0;
    label selectedHexCells = 0;
    label selectedWedgeCells = 0;

    const cellShapeList& cellShapes = mesh().cellShapes();

    forAll(zoneCells, i)
    {
        const label cellI = zoneCells[i];

        if (cellI < 0 || cellI >= mesh().nCells())
        {
            FatalErrorInFunction
                << "Invalid cell label " << cellI
                << " in apexRegularisation cellZone "
                << apexRegularisationZoneName_ << abort(FatalError);
        }

        if (maskI[cellI] < 0.5)
        {
            ++selectedCells;
            selectedVolume += mesh().V()[cellI];

            const word cellModelName = cellShapes[cellI].model().name();
            if (cellModelName == "hex")
            {
                ++selectedHexCells;
            }
            else if (cellModelName == "wedge")
            {
                ++selectedWedgeCells;
            }
        }

        maskI[cellI] = 1.0;
    }

    mask.correctBoundaryConditions();
    setBoundaryMaskToOwner(mask);
    validateBoundaryMaskOwnerValues(mask, "apexRegularisationMask");

    const label globalSelectedCells =
        returnReduce(selectedCells, sumOp<label>());

    if (globalSelectedCells != globalZoneSize)
    {
        FatalErrorInFunction
            << "apexRegularisation mask count " << globalSelectedCells
            << " does not match cellZone size " << globalZoneSize
            << abort(FatalError);
    }

    const scalar globalSelectedVolume =
        returnReduce(selectedVolume, sumOp<scalar>());
    const scalar totalVolume = gSum(mesh().V());
    const label globalHexCells =
        returnReduce(selectedHexCells, sumOp<label>());
    const label globalWedgeCells =
        returnReduce(selectedWedgeCells, sumOp<label>());

    apexActiveStressMaskPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "apexActiveStressMask",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                writeApexDiagnosticFields_
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (apexDisableActiveStress_)
    {
        apexActiveStressMaskPtr_() =
            dimensionedScalar("one", dimless, 1.0) - mask;
    }

    apexActiveStressMaskPtr_().correctBoundaryConditions();
    setBoundaryMaskToOwner(apexActiveStressMaskPtr_());
    validateBoundaryMaskOwnerValues
    (
        apexActiveStressMaskPtr_(),
        "apexActiveStressMask"
    );

    const tmp<surfaceScalarField> tMaskf =
        fvc::interpolate(apexActiveStressMaskPtr_());
    apexActiveStressMaskfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "apexActiveStressMaskf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tMaskf()
        )
    );
    setBoundaryMaskToOwner
    (
        apexActiveStressMaskfPtr_(),
        apexActiveStressMaskPtr_()
    );
    validateBoundaryMaskOwnerValues
    (
        apexActiveStressMaskfPtr_(),
        apexActiveStressMaskPtr_(),
        "apexActiveStressMaskf"
    );

    Info<< "electroMechanicalLaw apexRegularisation:" << nl
        << "    enabled = " << apexRegularisationEnabled_ << nl
        << "    cellZone = " << apexRegularisationZoneName_ << nl
        << "    selected cell count = " << globalSelectedCells << nl
        << "    selected reference volume = " << globalSelectedVolume << nl
        << "    selected volume fraction = "
        << globalSelectedVolume/(totalVolume + VSMALL) << nl
        << "    selected hex cells = " << globalHexCells << nl
        << "    selected wedge cells = " << globalWedgeCells << nl
        << "    disableActiveStress = " << apexDisableActiveStress_ << nl
        << "    writeDiagnosticFields = "
        << writeApexDiagnosticFields_ << endl;

    if (writeApexDiagnosticFields_)
    {
        Info<< "electroMechanicalLaw apexRegularisation boundary mask samples"
            << nl;
        reportBoundaryMaskSamples(mask, "apexRegularisationMask");
        reportBoundaryMaskSamples
        (
            apexActiveStressMaskPtr_(),
            "apexActiveStressMask"
        );
        reportBoundaryMaskSamples
        (
            apexActiveStressMaskfPtr_(),
            apexActiveStressMaskPtr_(),
            "apexActiveStressMaskf"
        );
    }

    writeApexDiagnosticFields();
}


const Foam::volScalarField&
Foam::electroMechanicalLaw::apexRegularisationMask() const
{
    if (!apexRegularisationMaskPtr_.valid())
    {
        makeApexRegularisationMasks();
    }

    return apexRegularisationMaskPtr_();
}


const Foam::volScalarField&
Foam::electroMechanicalLaw::apexActiveStressMask() const
{
    if (!apexActiveStressMaskPtr_.valid())
    {
        makeApexRegularisationMasks();
    }

    return apexActiveStressMaskPtr_();
}


const Foam::surfaceScalarField&
Foam::electroMechanicalLaw::apexActiveStressMaskf() const
{
    if (!apexActiveStressMaskfPtr_.valid())
    {
        makeApexRegularisationMasks();
    }

    return apexActiveStressMaskfPtr_();
}


void Foam::electroMechanicalLaw::writeApexDiagnosticFields() const
{
    if (writeApexDiagnosticFields_ && apexRegularisationMaskPtr_.valid())
    {
        apexRegularisationMaskPtr_().write();
        apexActiveStressMaskPtr_().write();
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
            passiveMechanicalLawType(dict),
            mesh,
            passiveMechanicalLawDict(dict),
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
	    f0f0_
	    (
	        IOobject
	        (
	            "f0f0",
	            mesh.time().timeName(),
	            mesh,
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        mesh,
	        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
	    ),
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
	    f0f0f_
	    (
	        IOobject
	        (
	            "f0f0f",
	            mesh.time().timeName(),
	            mesh,
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
	        mesh,
	        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
	    ),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    useFieldTa_(false),
    fieldTaChecked_(false),
    apexRegularisationEnabled_
    (
        readApexRegularisationSwitch(dict, "enabled", false)
    ),
    apexRegularisationZoneName_(readApexRegularisationZoneName(dict)),
    apexDisableActiveStress_
    (
        apexRegularisationEnabled_
      ? readApexRegularisationSwitch(dict, "disableActiveStress", false)
      : Switch(false)
    ),
    writeApexDiagnosticFields_
    (
        apexRegularisationEnabled_
      ? readApexRegularisationSwitch(dict, "writeDiagnosticFields", false)
      : Switch(false)
    ),
    apexRegularisationMaskPtr_(),
    apexActiveStressMaskPtr_(),
    apexActiveStressMaskfPtr_()
{
    if (rampTime_ < 0.0)
    {
        FatalErrorIn("electroMechanicalLaw::electroMechanicalLaw(...)")
            << "rampTime should be greater than or equal to zero"
            << abort(FatalError);
    }

    const scalar fibreMagnitudeTolerance
    (
        dict.lookupOrDefault<scalar>("fibreMagnitudeTolerance", SMALL)
    );

    MagnitudeStats f0Stats;
    MagnitudeStats f0fStats;

    normaliseVectorField
    (
        f0_,
        "electroMechanicalLaw f0",
        fibreMagnitudeTolerance,
        f0Stats
    );
    normaliseVectorField
    (
        f0f_,
        "electroMechanicalLaw f0f",
        fibreMagnitudeTolerance,
        f0fStats
    );

    Info<< "electroMechanicalLaw fibre normalisation" << nl;
    printMagnitudeStats("f0", f0Stats);
    printMagnitudeStats("f0f", f0fStats);

    f0f0_ = sqr(f0_);
    f0f0f_ = sqr(f0f_);

    if (isA<GuccioneElastic>(passiveMechLawPtr_()))
    {
        const GuccioneElastic& guccione =
            refCast<const GuccioneElastic>(passiveMechLawPtr_());

        const DifferenceStats cellStats =
            compareDyads(f0f0_, guccione.fibreDyad());
        const DifferenceStats faceStats =
            compareDyads(f0f0f_, guccione.faceFibreDyad());

        Info<< "electroMechanicalLaw/GuccioneElastic fibre dyad consistency"
            << nl
            << "    cell max Frobenius difference = "
            << cellStats.maxValue << nl
            << "    cell RMS Frobenius difference = "
            << cellStats.rms() << nl
            << "    face max Frobenius difference = "
            << faceStats.maxValue << nl
            << "    face RMS Frobenius difference = "
            << faceStats.rms() << nl;

        const Switch strictFibreConsistency
        (
            dict.lookupOrDefault<Switch>
            (
                "strictBenchmarkFibreConsistency",
                dict.lookupOrDefault<Switch>("strictFibreConsistency", false)
            )
        );
        const scalar fibreDyadTolerance
        (
            dict.lookupOrDefault<scalar>("fibreDyadTolerance", 1e-10)
        );

        if
        (
            strictFibreConsistency
         && (
                cellStats.maxValue > fibreDyadTolerance
             || faceStats.maxValue > fibreDyadTolerance
            )
        )
        {
            FatalErrorInFunction
                << "Active/passive fibre dyads differ above tolerance" << nl
                << "    fibreDyadTolerance = " << fibreDyadTolerance << nl
                << "    cell max = " << cellStats.maxValue << nl
                << "    face max = " << faceStats.maxValue
                << abort(FatalError);
        }
    }

    if (apexRegularisationEnabled_)
    {
        makeApexRegularisationMasks();
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


const Foam::volTensorField&
Foam::electroMechanicalLaw::deformationGradient() const
{
    return passiveMechLawPtr_->deformationGradient();
}


const Foam::surfaceTensorField&
Foam::electroMechanicalLaw::faceDeformationGradient() const
{
    return passiveMechLawPtr_->faceDeformationGradient();
}


bool Foam::electroMechanicalLaw::hasActiveStress() const
{
    checkFieldTa();

    if (useFieldTa_)
    {
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        return max(mag(Ta)).value() > SMALL;
    }

    return mag(activeTension().value()) > SMALL;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::electroMechanicalLaw::activeCauchyStress() const
{
    checkFieldTa();

    // Take a reference to the deformation gradient
    const volTensorField& F = deformationGradient();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    if (useFieldTa_)
    {
        // Field-based active tension from the coupling model
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
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
                apexDisableActiveStress_
              ? symm(F & ((Ta*apexActiveStressMask())*f0f0_) & F.T())/J
              : symm(F & (Ta*f0f0_) & F.T())/J
            )
        );
    }

    // Constant active tension with optional ramp
    const dimensionedScalar currentTa = activeTension();

    // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
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
            apexDisableActiveStress_
          ? symm
            (
                F
              & ((currentTa*apexActiveStressMask())*f0f0_)
              & F.T()
            )/J
          : symm(F & (currentTa*f0f0_) & F.T())/J
        )
    );
}


Foam::tmp<Foam::surfaceSymmTensorField>
Foam::electroMechanicalLaw::activeCauchyStressf() const
{
    checkFieldTa();

    // Take a reference to the deformation gradient
    const surfaceTensorField& F = faceDeformationGradient();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    if (useFieldTa_)
    {
        // Interpolate field-based active tension to faces
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        const surfaceScalarField Taf(fvc::interpolate(Ta));

        // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
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
                apexDisableActiveStress_
              ? (1.0/J)
               *symm
                (
                    F
                  & ((Taf*apexActiveStressMaskf())*f0f0f_)
                  & F.T()
                )
              : (1.0/J)*symm(F & (Taf*f0f0f_) & F.T())
            )
        );
    }

    // Constant active tension with optional ramp
    const dimensionedScalar currentTa = activeTension();

    // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
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
            apexDisableActiveStress_
          ? (1.0/J)
           *symm
            (
                F
              & ((currentTa*apexActiveStressMaskf())*f0f0f_)
              & F.T()
            )
          : (1.0/J)*symm(F & (currentTa*f0f0f_) & F.T())
        )
    );
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if
    (
        dict().lookupOrDefault<Switch>
        (
            "deformationGradientOwnershipAudit",
            false
        )
    )
    {
        const volTensorField& FActive = deformationGradient();
        const volTensorField& FPassive =
            passiveMechLawPtr_->deformationGradient();

        if (&FActive != &FPassive)
        {
            FatalErrorInFunction
                << "Active and passive laws do not reference the same "
                << "deformation-gradient object"
                << abort(FatalError);
        }

        if
        (
            nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN
         && !incremental()
        )
        {
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");
            scalar maxError = 0;

            forAll(FActive, cellI)
            {
                maxError = max
                (
                    maxError,
                    mag(FActive[cellI] - (I + gradD[cellI].T()))
                );
            }

            forAll(FActive.boundaryField(), patchI)
            {
                const tensorField& FPatch =
                    FActive.boundaryField()[patchI];
                const tensorField& gradDPatch =
                    gradD.boundaryField()[patchI];

                forAll(FPatch, faceI)
                {
                    maxError = max
                    (
                        maxError,
                        mag
                        (
                            FPatch[faceI]
                          - (I + gradDPatch[faceI].T())
                        )
                    );
                }
            }

            reduce(maxError, maxOp<scalar>());

            if (maxError > 100*SMALL)
            {
                FatalErrorInFunction
                    << "Active/passive deformation gradient differs from "
                    << "I + grad(D).T(): max error = " << maxError
                    << abort(FatalError);
            }

            static label auditEvaluations = 0;
            static label reportedTimeIndex = -1;
            ++auditEvaluations;

            if (reportedTimeIndex != mesh().time().timeIndex())
            {
                reportedTimeIndex = mesh().time().timeIndex();
                Info<< "electroMechanicalLaw deformation-gradient invariant: "
                    << "time = " << mesh().time().value()
                    << ", evaluation = " << auditEvaluations
                    << ", active/passive same object = true"
                    << ", max |F_active - (I + grad(D).T())| = "
                    << maxError << nl;
            }
        }
    }

    if (hasActiveStress())
    {
        const tmp<volSymmTensorField> tSigmaActive = activeCauchyStress();
        sigma += tSigmaActive();
    }
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if (hasActiveStress())
    {
        const tmp<surfaceSymmTensorField> tSigmaActive = activeCauchyStressf();
        sigma += tSigmaActive();
    }
}


void Foam::electroMechanicalLaw::updateTotalFields()
{
    passiveMechLawPtr_->updateTotalFields();
}


void Foam::electroMechanicalLaw::setRestart()
{
    passiveMechLawPtr_->setRestart();
}


// ************************************************************************* //
