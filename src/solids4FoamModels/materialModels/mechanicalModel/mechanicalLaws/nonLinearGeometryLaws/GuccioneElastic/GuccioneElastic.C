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

#include "GuccioneElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "eig3.H"
#include "ZoneIDs.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GuccioneElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GuccioneElastic, nonLinGeomMechLaw
    );
}

namespace
{

using namespace Foam;

Foam::dimensionedScalar readBulkModulus
(
    const Foam::dictionary& dict
)
{
    if
    (
        dict.lookupOrDefault<Foam::Switch>
        (
            "pressureDisplacement",
            Foam::Switch(false)
        )
     && !dict.found("bulkModulus")
    )
    {
        return Foam::dimensionedScalar("K", Foam::dimPressure, Foam::GREAT);
    }

    return Foam::dimensionedScalar(dict.lookup("bulkModulus"));
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


Foam::scalar readApexIsotropicB(const Foam::dictionary& dict)
{
    if (!hasApexRegularisationDict(dict))
    {
        return 2.0;
    }

    return
        apexRegularisationDict(dict).lookupOrDefault<Foam::scalar>
        (
            "isotropicB",
            2.0
        );
}


Foam::IOobject findFibreFieldIOobject
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
                << exit(Foam::FatalError);
        }
    }

    io.readOpt() = Foam::IOobject::MUST_READ;

    return io;
}


Foam::IOobject findSurfaceFibreFieldIOobject
(
    const Foam::word& fieldName,
    const Foam::fvMesh& mesh,
    const bool required
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
    bool ok = io.headerOk();
#elif defined(OPENFOAM_ORG)
    bool ok = io.typeHeaderOk<Foam::surfaceVectorField>(true);
#else
    bool ok = io.typeHeaderOk<Foam::surfaceVectorField>(true, false, false);
#endif

    if (!ok)
    {
        io.instance() = "0";

#ifdef FOAMEXTEND
        ok = io.headerOk();
#elif defined(OPENFOAM_ORG)
        ok = io.typeHeaderOk<Foam::surfaceVectorField>(true);
#else
        ok = io.typeHeaderOk<Foam::surfaceVectorField>(true, false, false);
#endif
    }

    if (!ok && required)
    {
        FatalErrorInFunction
            << "Cannot find required surface fibre field " << fieldName
            << " in either " << mesh.time().timeName() << " or 0"
            << exit(Foam::FatalError);
    }

    io.readOpt() = ok ? Foam::IOobject::MUST_READ : Foam::IOobject::NO_READ;

    return io;
}


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


Foam::tmp<Foam::surfaceVectorField> signAwareInterpolateF0
(
    const Foam::volVectorField& f0,
    const Foam::fvMesh& mesh
)
{
    Foam::tmp<Foam::surfaceVectorField> tResult
    (
        new Foam::surfaceVectorField
        (
            Foam::IOobject
            (
                "f0fSignAwareInterpolated",
                mesh.time().timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE
            ),
            mesh,
            Foam::dimensionedVector("zero", Foam::dimless, Foam::vector::zero)
        )
    );

    Foam::surfaceVectorField& result = tResult.ref();

#ifdef OPENFOAM_NOT_EXTEND
    Foam::vectorField& resultI = result.primitiveFieldRef();
#else
    Foam::vectorField& resultI = result.internalField();
#endif

    const Foam::labelUList& owner = mesh.owner();
    const Foam::labelUList& neighbour = mesh.neighbour();

    forAll(resultI, faceI)
    {
        Foam::vector own = f0[owner[faceI]];
        Foam::vector nei = f0[neighbour[faceI]];
        if ((own & nei) < 0)
        {
            nei = -nei;
        }
        resultI[faceI] = own + nei;
    }

    forAll(result.boundaryField(), patchI)
    {
        Foam::vectorField& patch = result.boundaryFieldRef()[patchI];
        const Foam::labelUList& faceCells = mesh.boundary()[patchI].faceCells();
        forAll(patch, faceI)
        {
            Foam::vector boundaryValue = f0.boundaryField()[patchI][faceI];
            const Foam::vector own = f0[faceCells[faceI]];
            if ((own & boundaryValue) < 0)
            {
                boundaryValue = -boundaryValue;
            }
            patch[faceI] = own + boundaryValue;
        }
    }

    return tResult;
}


Foam::vector seedAxis(const Foam::vector& f)
{
    Foam::vector seed(1, 0, 0);
    Foam::scalar minDot = Foam::mag(f.x());

    if (Foam::mag(f.y()) < minDot)
    {
        seed = Foam::vector(0, 1, 0);
        minDot = Foam::mag(f.y());
    }

    if (Foam::mag(f.z()) < minDot)
    {
        seed = Foam::vector(0, 0, 1);
    }

    return seed;
}


Foam::tensor basisTensor
(
    const Foam::vector& f,
    const Foam::vector& s,
    const Foam::vector& n
)
{
    return Foam::tensor
    (
        f.x(), s.x(), n.x(),
        f.y(), s.y(), n.y(),
        f.z(), s.z(), n.z()
    );
}


void constructBasis
(
    const Foam::vector& f,
    const Foam::scalar tolerance,
    Foam::vector& s,
    Foam::vector& n,
    Foam::tensor& R
)
{
    const Foam::vector seed = seedAxis(f);
    s = seed - (seed & f)*f;

    if (Foam::mag(s) <= tolerance)
    {
        FatalErrorInFunction
            << "Failed to construct Guccione sheet direction" << nl
            << "    f0 = " << f << nl
            << "    seed = " << seed << abort(Foam::FatalError);
    }

    s /= Foam::mag(s);
    n = f ^ s;

    if (Foam::mag(n) <= tolerance)
    {
        FatalErrorInFunction
            << "Failed to construct Guccione sheet-normal direction" << nl
            << "    f0 = " << f << nl
            << "    s0 = " << s << abort(Foam::FatalError);
    }

    n /= Foam::mag(n);
    s = n ^ f;
    s /= Foam::mag(s);
    R = basisTensor(f, s, n);
}


struct BasisStats
{
    Foam::scalar magF;
    Foam::scalar magS;
    Foam::scalar magN;
    Foam::scalar fdots;
    Foam::scalar fdotn;
    Foam::scalar sdotn;
    Foam::scalar RtR;
    Foam::scalar detR;

    BasisStats()
    :
        magF(0),
        magS(0),
        magN(0),
        fdots(0),
        fdotn(0),
        sdotn(0),
        RtR(0),
        detR(0)
    {}

    void add
    (
        const Foam::vector& f,
        const Foam::vector& s,
        const Foam::vector& n,
        const Foam::tensor& R
    )
    {
        magF = Foam::max(magF, Foam::mag(Foam::mag(f) - 1));
        magS = Foam::max(magS, Foam::mag(Foam::mag(s) - 1));
        magN = Foam::max(magN, Foam::mag(Foam::mag(n) - 1));
        fdots = Foam::max(fdots, Foam::mag(f & s));
        fdotn = Foam::max(fdotn, Foam::mag(f & n));
        sdotn = Foam::max(sdotn, Foam::mag(s & n));
        RtR = Foam::max(RtR, Foam::mag((R.T() & R) - Foam::tensor::I));
        detR = Foam::max(detR, Foam::mag(Foam::det(R) - 1));
    }
};


void printBasisStats(const Foam::word& name, const BasisStats& stats)
{
    Info<< "    " << name << " max |f0|-1 = " << stats.magF << nl
        << "    " << name << " max |s0|-1 = " << stats.magS << nl
        << "    " << name << " max |n0|-1 = " << stats.magN << nl
        << "    " << name << " max |f0.s0| = " << stats.fdots << nl
        << "    " << name << " max |f0.n0| = " << stats.fdotn << nl
        << "    " << name << " max |s0.n0| = " << stats.sdotn << nl
        << "    " << name << " max ||R^T R-I|| = " << stats.RtR << nl
        << "    " << name << " max |det(R)-1| = " << stats.detR << nl;
}

}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::GuccioneElastic::makeF0
(
    const Switch& uniformFibreField,
    const fvMesh& mesh,
    const dictionary& dict
) const
{
    if (uniformFibreField)
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "f0",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("f0", dimless, dict.lookup("f0"))
            )
        );
    }

    return tmp<volVectorField>
    (
        new volVectorField
        (
            findFibreFieldIOobject("f0", mesh),
            mesh
        )
    );
}


Foam::tmp<Foam::surfaceVectorField> Foam::GuccioneElastic::makeF0f
(
    const Switch& uniformFibreField,
    const fvMesh& mesh,
    const dictionary& dict,
    const volVectorField& f0
) const
{
    if (uniformFibreField)
    {
        return tmp<surfaceVectorField>
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "f0f",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("f0", dimless, dict.lookup("f0"))
            )
        );
    }

    IOobject f0fIO(findSurfaceFibreFieldIOobject("f0f", mesh, false));

    if (f0fIO.readOpt() == IOobject::MUST_READ)
    {
        return tmp<surfaceVectorField>
        (
            new surfaceVectorField
            (
                f0fIO,
                mesh
            )
        );
    }

    const Switch allowInterpolatedFaceFibres
    (
        dict.lookupOrDefault<Switch>
        (
            "allowInterpolatedFaceFibres",
            Switch(false)
        )
    );

    if (!allowInterpolatedFaceFibres)
    {
        FatalErrorInFunction
            << "Cannot find required surfaceVectorField f0f. "
            << "GuccioneElastic no longer interpolates f0 to faces unless "
            << "allowInterpolatedFaceFibres true is set explicitly."
            << abort(FatalError);
    }

    WarningInFunction
        << "Interpolating f0 to create f0f because "
        << "allowInterpolatedFaceFibres true. This is an explicit fallback; "
        << "use the analytic supplied f0f field for Land2015 benchmark runs."
        << endl;

    return signAwareInterpolateF0(f0, mesh);
}


bool Foam::GuccioneElastic::useApexPassiveIsotropisation() const
{
    return
        apexRegularisationEnabled_
     && apexIsotropisePassive_;
}


bool Foam::GuccioneElastic::useApexImplicitStiffness() const
{
    return
        useApexPassiveIsotropisation()
     && !preserveBaselineImpK_;
}


void Foam::GuccioneElastic::makeApexRegularisationMasks() const
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
                "apexGuccioneRegularisationMask",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
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
        }

        maskI[cellI] = 1.0;
    }

    mask.correctBoundaryConditions();
    setBoundaryMaskToOwner(mask);
    validateBoundaryMaskOwnerValues(mask, "apexGuccioneRegularisationMask");

    const label globalSelectedCells =
        returnReduce(selectedCells, sumOp<label>());

    if (globalSelectedCells != globalZoneSize)
    {
        FatalErrorInFunction
            << "apexRegularisation mask count " << globalSelectedCells
            << " does not match cellZone size " << globalZoneSize
            << abort(FatalError);
    }

    const tmp<surfaceScalarField> tMaskf = fvc::interpolate(mask);
    apexRegularisationMaskfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "apexGuccioneRegularisationMaskf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tMaskf()
        )
    );
    setBoundaryMaskToOwner(apexRegularisationMaskfPtr_(), mask);
    validateBoundaryMaskOwnerValues
    (
        apexRegularisationMaskfPtr_(),
        mask,
        "apexGuccioneRegularisationMaskf"
    );

    apexPassiveIsotropicMaskPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "apexPassiveIsotropicMask",
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

    if (useApexPassiveIsotropisation())
    {
        apexPassiveIsotropicMaskPtr_() = mask;
    }
    apexPassiveIsotropicMaskPtr_().correctBoundaryConditions();
    setBoundaryMaskToOwner(apexPassiveIsotropicMaskPtr_());
    validateBoundaryMaskOwnerValues
    (
        apexPassiveIsotropicMaskPtr_(),
        "apexPassiveIsotropicMask"
    );

    apexBaselineImpKMaskPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "apexBaselineImpKMask",
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

    if (useApexPassiveIsotropisation() && preserveBaselineImpK_)
    {
        apexBaselineImpKMaskPtr_() = mask;
    }
    apexBaselineImpKMaskPtr_().correctBoundaryConditions();
    setBoundaryMaskToOwner(apexBaselineImpKMaskPtr_());
    validateBoundaryMaskOwnerValues
    (
        apexBaselineImpKMaskPtr_(),
        "apexBaselineImpKMask"
    );

    const scalar globalSelectedVolume =
        returnReduce(selectedVolume, sumOp<scalar>());
    const dimensionedScalar muIso
    (
        "muIsoApex",
        dimPressure,
        0.5*k_.value()*apexIsotropicB_
    );

    Info<< "GuccioneElastic apexRegularisation:" << nl
        << "    enabled = " << apexRegularisationEnabled_ << nl
        << "    cellZone = " << apexRegularisationZoneName_ << nl
        << "    selected cell count = " << globalSelectedCells << nl
        << "    selected reference volume = " << globalSelectedVolume << nl
        << "    isotropisePassive = " << apexIsotropisePassive_ << nl
        << "    isotropicB = " << apexIsotropicB_ << nl
        << "    preserveBaselineImpK = " << preserveBaselineImpK_ << nl
        << "    original cf/ct/cfs = "
        << cf_ << " / " << ct_ << " / " << cfs_ << nl
        << "    apex effective cf/ct/cfs = "
        << (
               useApexPassiveIsotropisation() ? apexIsotropicB_ : cf_
           )
        << " / "
        << (
               useApexPassiveIsotropisation() ? apexIsotropicB_ : ct_
           )
        << " / "
        << (
               useApexPassiveIsotropisation() ? apexIsotropicB_ : cfs_
           )
        << nl
        << "    baseline small-strain shear modulus = " << mu_.value() << nl
        << "    apex isotropic small-strain shear modulus = "
        << muIso.value() << endl;

    if (useApexPassiveIsotropisation() && preserveBaselineImpK_)
    {
        Info<< "    baseline impK/shearModulus will be preserved for "
            << "momentum stabilisation and approximate Jacobian" << endl;
    }

    if (useApexImplicitStiffness())
    {
        WarningInFunction
            << "preserveBaselineImpK is false: this experiment changes "
            << "the effective shear modulus, impK, momentum-stabilisation "
            << "coefficient and approximate displacement Jacobian"
            << endl;
    }

    if (writeApexDiagnosticFields_)
    {
        Info<< "GuccioneElastic apexRegularisation boundary mask samples"
            << nl;
        reportBoundaryMaskSamples(mask, "apexGuccioneRegularisationMask");
        reportBoundaryMaskSamples
        (
            apexRegularisationMaskfPtr_(),
            mask,
            "apexGuccioneRegularisationMaskf"
        );
        reportBoundaryMaskSamples
        (
            apexPassiveIsotropicMaskPtr_(),
            "apexPassiveIsotropicMask"
        );
        reportBoundaryMaskSamples
        (
            apexBaselineImpKMaskPtr_(),
            "apexBaselineImpKMask"
        );
    }

    writeApexDiagnosticFields();
}


const Foam::volScalarField&
Foam::GuccioneElastic::apexRegularisationMask() const
{
    if (!apexRegularisationMaskPtr_.valid())
    {
        makeApexRegularisationMasks();
    }

    return apexRegularisationMaskPtr_();
}


const Foam::surfaceScalarField&
Foam::GuccioneElastic::apexRegularisationMaskf() const
{
    if (!apexRegularisationMaskfPtr_.valid())
    {
        makeApexRegularisationMasks();
    }

    return apexRegularisationMaskfPtr_();
}


void Foam::GuccioneElastic::writeApexDiagnosticFields() const
{
    if (writeApexDiagnosticFields_ && apexPassiveIsotropicMaskPtr_.valid())
    {
        apexPassiveIsotropicMaskPtr_().write();
        apexBaselineImpKMaskPtr_().write();
    }
}


void Foam::GuccioneElastic::reportApexRegularisationCoefficients() const
{
    if (!apexRegularisationEnabled_)
    {
        return;
    }

    const volScalarField& mask = apexRegularisationMask();

#ifdef OPENFOAM_NOT_EXTEND
    const scalarField& maskI = mask.primitiveField();
    const scalarField& muEffI = muEff_.primitiveField();
#else
    const scalarField& maskI = mask.internalField();
    const scalarField& muEffI = muEff_.internalField();
#endif

    scalar insideVolume = 0.0;
    scalar outsideVolume = 0.0;
    scalar insideMuEff = 0.0;
    scalar outsideMuEff = 0.0;
    scalar insideImpK = 0.0;
    scalar outsideImpK = 0.0;

    forAll(maskI, cellI)
    {
        const scalar m = maskI[cellI];
        const scalar V = mesh().V()[cellI];
        const scalar muValue =
            pressureDisplacement_
          ? muEffI[cellI]
          : 0.5*k_.value()*implicitCoefficient
            (
                (cf_ + cfs_ + ct_)/3.0,
                cellI
            );
        const scalar impKValue =
            pressureDisplacement_
          ? impKcoeff_*muEffI[cellI]
          : (4.0/3.0)*muValue + bulkModulus_.value();

        insideVolume += V*m;
        outsideVolume += V*(1.0 - m);
        insideMuEff += V*m*muValue;
        outsideMuEff += V*(1.0 - m)*muValue;
        insideImpK += V*m*impKValue;
        outsideImpK += V*(1.0 - m)*impKValue;
    }

    insideVolume = returnReduce(insideVolume, sumOp<scalar>());
    outsideVolume = returnReduce(outsideVolume, sumOp<scalar>());
    insideMuEff = returnReduce(insideMuEff, sumOp<scalar>());
    outsideMuEff = returnReduce(outsideMuEff, sumOp<scalar>());
    insideImpK = returnReduce(insideImpK, sumOp<scalar>());
    outsideImpK = returnReduce(outsideImpK, sumOp<scalar>());

    const scalar physicalCfInside =
        useApexPassiveIsotropisation() ? apexIsotropicB_ : cf_;
    const scalar physicalCtInside =
        useApexPassiveIsotropisation() ? apexIsotropicB_ : ct_;
    const scalar physicalCfsInside =
        useApexPassiveIsotropisation() ? apexIsotropicB_ : cfs_;

    const scalar implicitCfInside =
        useApexImplicitStiffness() ? apexIsotropicB_ : cf_;
    const scalar implicitCtInside =
        useApexImplicitStiffness() ? apexIsotropicB_ : ct_;
    const scalar implicitCfsInside =
        useApexImplicitStiffness() ? apexIsotropicB_ : cfs_;

    Info<< "GuccioneElastic apexRegularisation coefficient summary:" << nl
        << "    outside physical cf/ct/cfs = "
        << cf_ << " / " << ct_ << " / " << cfs_ << nl
        << "    inside physical cf/ct/cfs = "
        << physicalCfInside << " / "
        << physicalCtInside << " / "
        << physicalCfsInside << nl
        << "    outside implicit cf/ct/cfs = "
        << cf_ << " / " << ct_ << " / " << cfs_ << nl
        << "    inside implicit cf/ct/cfs = "
        << implicitCfInside << " / "
        << implicitCtInside << " / "
        << implicitCfsInside << nl
        << "    inside average effective shear modulus = "
        << insideMuEff/(insideVolume + VSMALL) << nl
        << "    outside average effective shear modulus = "
        << outsideMuEff/(outsideVolume + VSMALL) << nl
        << "    inside average impK = "
        << insideImpK/(insideVolume + VSMALL) << nl
        << "    outside average impK = "
        << outsideImpK/(outsideVolume + VSMALL) << endl;
}


Foam::scalar Foam::GuccioneElastic::implicitCoefficient
(
    const scalar baseline,
    const label cellI
) const
{
    if (useApexImplicitStiffness() && apexRegularisationMask()[cellI] > 0.5)
    {
        return apexIsotropicB_;
    }

    return baseline;
}


void Foam::GuccioneElastic::calculateStress
(
    surfaceSymmTensorField& sigma,
    const surfaceTensorField& gradD
)
{
    const surfaceTensorField Ff("FfFromGradD", I + gradD.T());
    const surfaceScalarField Jf("Jf", det(Ff));

    surfaceTensorField FfWork("FfWork", Ff);
    if (useIsochoricSplit_)
    {
        FfWork = pow(Jf, -1.0/3.0)*Ff;
    }
    const surfaceTensorField FfWorkT("FfWorkT", FfWork.T());

    // Calculate the right Cauchy-Green deformation tensor
    const surfaceSymmTensorField C("C", symm(FfWorkT & FfWork));

    // Calculate the Green-Lagrange strain
    const surfaceSymmTensorField E("E", 0.5*(C - I));

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const surfaceTensorField RfT("RfT", Rf_.T());
        const surfaceSymmTensorField EStar("EStar", symm(RfT & E & Rf_));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const surfaceScalarField E11
        (
            "E11", EStar.component(symmTensor::XX)
        );
        const surfaceScalarField E12
        (
            "E12", EStar.component(symmTensor::XY)
        );
        const surfaceScalarField E13
        (
            "E13", EStar.component(symmTensor::XZ)
        );
        const surfaceScalarField E22
        (
            "E22", EStar.component(symmTensor::YY)
        );
        const surfaceScalarField E23
        (
            "E23", EStar.component(symmTensor::YZ)
        );
        const surfaceScalarField E33
        (
            "E33", EStar.component(symmTensor::ZZ)
        );

        if (useApexPassiveIsotropisation())
        {
            surfaceScalarField cf
            (
                IOobject
                (
                    "cfApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cf", dimless, cf_)
            );
            surfaceScalarField ct
            (
                IOobject
                (
                    "ctApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("ct", dimless, ct_)
            );
            surfaceScalarField cfs
            (
                IOobject
                (
                    "cfsApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cfs", dimless, cfs_)
            );

            cf += apexRegularisationMaskf()*(apexIsotropicB_ - cf_);
            ct += apexRegularisationMaskf()*(apexIsotropicB_ - ct_);
            cfs += apexRegularisationMaskf()*(apexIsotropicB_ - cfs_);

            // Calculate Q
            const surfaceScalarField Q
            (
                "Q",
                cf*sqr(E11)
              + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
              + cfs*(2*sqr(E12) + 2*sqr(E13))
            );

            // Calculate the derivative of Q wrt to EStar
            surfaceSymmTensorField dQdEStar
            (
                IOobject
                (
                    "dQdEStar",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedSymmTensor("0", dimless, symmTensor::zero)
            );

            dQdEStar.replace(symmTensor::XX, 2*cf*E11);
            dQdEStar.replace(symmTensor::XY, 2*cfs*E12);
            dQdEStar.replace(symmTensor::XZ, 2*cfs*E13);
            dQdEStar.replace(symmTensor::YY, 2*ct*E22);
            dQdEStar.replace(symmTensor::YZ, 2*ct*E23);
            dQdEStar.replace(symmTensor::ZZ, 2*ct*E33);

            const surfaceScalarField expQ("expQ", exp(Q));

            // Calculate the local 2nd Piola-Kirchhoff stress (without the
            // hydrostatic term)
            Sf_ = dQdEStar*0.5*k_*expQ;
        }
        else
        {
            // Calculate Q
            const surfaceScalarField Q
            (
                "Q",
                cf_*sqr(E11)
              + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
              + cfs_*(2*sqr(E12) + 2*sqr(E13))
            );

            // Calculate the derivative of Q wrt to EStar
            surfaceSymmTensorField dQdEStar
            (
                IOobject
                (
                    "dQdEStar",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedSymmTensor("0", dimless, symmTensor::zero)
            );

            dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
            dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
            dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
            dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
            dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
            dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

            const surfaceScalarField expQ("expQ", exp(Q));

            // Calculate the local 2nd Piola-Kirchhoff stress (without the
            // hydrostatic term)
            Sf_ = dQdEStar*0.5*k_*expQ;
        }

        // Rotate S from the local fibre coordinate system to the global
        // coordinate system
        Sf_ = symm(Rf_ & Sf_ & RfT);
    }
    else
    {
        // Calculate E . E
        const surfaceSymmTensorField sqrE("sqrE", symm(E & E));

        // Calculate the invariants of E
        const surfaceScalarField I1("I1", tr(E));
        const surfaceScalarField I2
        (
            "I2",
            0.5*(sqr(tr(E)) - tr(sqrE))
        );
        const surfaceScalarField I4("I4", E && f0f0f_);
        const surfaceScalarField I5("I5", sqrE && f0f0f_);

        if (useApexPassiveIsotropisation())
        {
            surfaceScalarField cf
            (
                IOobject
                (
                    "cfApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cf", dimless, cf_)
            );
            surfaceScalarField ct
            (
                IOobject
                (
                    "ctApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("ct", dimless, ct_)
            );
            surfaceScalarField cfs
            (
                IOobject
                (
                    "cfsApexf",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cfs", dimless, cfs_)
            );

            cf += apexRegularisationMaskf()*(apexIsotropicB_ - cf_);
            ct += apexRegularisationMaskf()*(apexIsotropicB_ - ct_);
            cfs += apexRegularisationMaskf()*(apexIsotropicB_ - cfs_);

            // Calculate Q
            const surfaceScalarField Q
            (
                "Q",
                ct*sqr(I1)
              - 2.0*ct*I2
              + (cf - 2.0*cfs + ct)*sqr(I4)
              + 2.0*(cfs - ct)*I5
            );

            // Calculate the derivative of Q wrt to E
            const surfaceSymmTensorField dQdE
            (
                2.0*ct*E
              + 2.0*(cf - 2.0*cfs + ct)*I4*f0f0f_
              + 2.0*(cfs - ct)*symm((E & f0f0f_) + (f0f0f_ & E))
            );

            const surfaceScalarField expQ("expQ", exp(Q));

            // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic
            // term)
            Sf_ = dQdE*0.5*k_*expQ;
        }
        else
        {
            // Calculate Q
            const surfaceScalarField Q
            (
                "Q",
                ct_*sqr(I1)
              - 2.0*ct_*I2
              + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
              + 2.0*(cfs_ - ct_)*I5
            );

            // Calculate the derivative of Q wrt to E
            const surfaceSymmTensorField dQdE
            (
                2.0*ct_*E
              + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0f_
              + 2.0*(cfs_ - ct_)*symm((E & f0f0f_) + (f0f0f_ & E))
            );

            const surfaceScalarField expQ("expQ", exp(Q));

            // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic
            // term)
            Sf_ = dQdE*0.5*k_*expQ;
        }
    }

    // Convert the second Piola-Kirchhoff stress to the deviatoric Cauchy
    // stress
    const surfaceSymmTensorField sf
    (
        "sf",
        dev(symm(FfWork & Sf_ & FfWorkT))/Jf
    );

    if (pressureDisplacement_)
    {
        // During solid-model construction pf may not be registered yet.
        // The mixed solid model applies the pressure split once it is available.
        if (!mesh().foundObject<surfaceScalarField>("pf"))
        {
            sigma = sf;
            return;
        }

        // Add the pressure-displacement hydrostatic term
        const surfaceScalarField& pf =
            mesh().lookupObject<surfaceScalarField>("pf");
        sigma = sf - pf*I;
        return;
    }

    sigma =
        sf
      + (
            0.5*bulkModulus_*(pow(Jf, 2.0) - 1.0)/Jf
        )*I;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GuccioneElastic::GuccioneElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    bulkModulus_(readBulkModulus(dict)),
    k_(dict.lookup("k")),
    cf_(readScalar(dict.lookup("cf"))),
    ct_(readScalar(dict.lookup("ct"))),
    cfs_(readScalar(dict.lookup("cfs"))),
    // Linearised (small-strain) shear modulus: average of the three
    // material constants scaled by k. Same form is used in both the
    // pressure-displacement and standard branches.
    mu_(0.5*k_*(cf_ + cfs_ + ct_)/3.0),
    pressureDisplacement_
    (
        dict.lookupOrDefault<Switch>("pressureDisplacement", false)
    ),
    useIsochoricSplit_
    (
        dict.lookupOrDefault<Switch>
        (
            "useIsochoricSplit",
            false
        )
    ),
    apexRegularisationEnabled_
    (
        readApexRegularisationSwitch(dict, "enabled", false)
    ),
    apexRegularisationZoneName_(readApexRegularisationZoneName(dict)),
    apexIsotropisePassive_
    (
        apexRegularisationEnabled_
      ? readApexRegularisationSwitch(dict, "isotropisePassive", false)
      : Switch(false)
    ),
    apexIsotropicB_(readApexIsotropicB(dict)),
    preserveBaselineImpK_
    (
        apexRegularisationEnabled_
      ? readApexRegularisationSwitch(dict, "preserveBaselineImpK", true)
      : Switch(true)
    ),
    writeApexDiagnosticFields_
    (
        apexRegularisationEnabled_
      ? readApexRegularisationSwitch(dict, "writeDiagnosticFields", false)
      : Switch(false)
    ),
    apexRegularisationMaskPtr_(),
    apexRegularisationMaskfPtr_(),
    apexPassiveIsotropicMaskPtr_(),
    apexBaselineImpKMaskPtr_(),
    muEff_
    (
        IOobject
        (
            "muEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        mu_,
        zeroGradientFvPatchScalarField::typeName
    ),
    uniformFibreField_
    (
        dict.lookupOrDefault<Switch>("uniformFibreField", false)
    ),
    f0_(makeF0(uniformFibreField_, mesh, dict)),
    f0f_(makeF0f(uniformFibreField_, mesh, dict, f0_)),
    s0_
    (
        IOobject
        (
            "s0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("i", dimless, vector(1, 0, 0))
    ),
    s0f_
    (
        IOobject
        (
            "s0f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("i", dimless, vector(1, 0, 0))
    ),
    n0_
    (
        IOobject
        (
            "n0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    ),
    n0f_
    (
        IOobject
        (
            "n0f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    ),
    R_
    (
        IOobject
        (
            "R",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    Rf_
    (
        IOobject
        (
            "Rf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
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
    S_
    (
        IOobject
        (
            "S2PK",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
    ),
    Sf_
    (
        IOobject
        (
            "S2PKf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
    ),
    expQf_
    (
        IOobject
        (
            "expQf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("1", dimless, 1)
    ),
    impKcoeff_
    (
        dict.lookupOrDefault<scalar>("impKcoeff", 1.0)
    )
{
    Info<< "GuccioneElastic:" << nl
        << "    useIsochoricSplit = " << useIsochoricSplit_ << endl;

    if (pressureDisplacement_)
    {
        muEff_ = mu_;
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
        "GuccioneElastic f0",
        fibreMagnitudeTolerance,
        f0Stats
    );
    normaliseVectorField
    (
        f0f_,
        "GuccioneElastic f0f",
        fibreMagnitudeTolerance,
        f0fStats
    );

    Info<< "GuccioneElastic fibre normalisation" << nl;
    printMagnitudeStats("f0", f0Stats);
    printMagnitudeStats("f0f", f0fStats);

    f0f0_ = sqr(f0_);
    f0f0f_ = sqr(f0f_);

    surfaceVectorField f0fInterpolated
    (
        signAwareInterpolateF0(f0_, mesh)
    );
    MagnitudeStats f0fInterpolatedStats;
    normaliseVectorField
    (
        f0fInterpolated,
        "GuccioneElastic sign-aware interpolated f0f diagnostic",
        fibreMagnitudeTolerance,
        f0fInterpolatedStats
    );

    surfaceScalarField suppliedVsInterpolatedFaceFibreAngle
    (
        IOobject
        (
            "suppliedVsInterpolatedFaceFibreAngle",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    scalar angleMin = VGREAT;
    scalar angleMax = -VGREAT;
    scalar angleSum = 0;
    scalar angleSumSqr = 0;
    label angleN = 0;

    const scalar radToDeg = 180.0/constant::mathematical::pi;

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& angleI =
        suppliedVsInterpolatedFaceFibreAngle.primitiveFieldRef();
    const vectorField& f0fI = f0f_.primitiveField();
    const vectorField& interpI = f0fInterpolated.primitiveField();
#else
    scalarField& angleI =
        suppliedVsInterpolatedFaceFibreAngle.internalField();
    const vectorField& f0fI = f0f_.internalField();
    const vectorField& interpI = f0fInterpolated.internalField();
#endif

    forAll(angleI, faceI)
    {
        const scalar c =
            min(scalar(1), max(scalar(-1), mag(f0fI[faceI] & interpI[faceI])));
        const scalar angle = acos(c)*radToDeg;
        angleI[faceI] = angle;
        angleMin = min(angleMin, angle);
        angleMax = max(angleMax, angle);
        angleSum += angle;
        angleSumSqr += sqr(angle);
        ++angleN;
    }

    forAll(suppliedVsInterpolatedFaceFibreAngle.boundaryField(), patchI)
    {
        scalarField& patchAngle =
            suppliedVsInterpolatedFaceFibreAngle.boundaryFieldRef()[patchI];
        const vectorField& patchF0f = f0f_.boundaryField()[patchI];
        const vectorField& patchInterp =
            f0fInterpolated.boundaryField()[patchI];

        forAll(patchAngle, faceI)
        {
            const scalar c =
                min
                (
                    scalar(1),
                    max(scalar(-1), mag(patchF0f[faceI] & patchInterp[faceI]))
                );
            const scalar angle = acos(c)*radToDeg;
            patchAngle[faceI] = angle;
            angleMin = min(angleMin, angle);
            angleMax = max(angleMax, angle);
            angleSum += angle;
            angleSumSqr += sqr(angle);
            ++angleN;
        }
    }

    Info<< "GuccioneElastic supplied-vs-interpolated face fibre angle "
        << "statistics [deg]: n = " << angleN
        << ", min = " << angleMin
        << ", max = " << angleMax
        << ", mean = " << (angleN ? angleSum/scalar(angleN) : 0)
        << ", RMS = " << (angleN ? sqrt(angleSumSqr/scalar(angleN)) : 0)
        << nl;
    suppliedVsInterpolatedFaceFibreAngle.write();

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();

    const scalar basisTolerance
    (
        dict.lookupOrDefault<scalar>("fibreBasisTolerance", 1e-12)
    );
    BasisStats cellBasisStats;
    BasisStats faceBasisStats;

#ifdef OPENFOAM_NOT_EXTEND
    vectorField& f0Internal = f0_.primitiveFieldRef();
    vectorField& s0Internal = s0_.primitiveFieldRef();
    vectorField& n0Internal = n0_.primitiveFieldRef();
    tensorField& RInternal = R_.primitiveFieldRef();
#else
    vectorField& f0Internal = f0_.internalField();
    vectorField& s0Internal = s0_.internalField();
    vectorField& n0Internal = n0_.internalField();
    tensorField& RInternal = R_.internalField();
#endif

    forAll(f0Internal, cellI)
    {
        constructBasis
        (
            f0Internal[cellI],
            basisTolerance,
            s0Internal[cellI],
            n0Internal[cellI],
            RInternal[cellI]
        );
        cellBasisStats.add
        (
            f0Internal[cellI],
            s0Internal[cellI],
            n0Internal[cellI],
            RInternal[cellI]
        );
    }

    forAll(f0_.boundaryField(), patchI)
    {
        const vectorField& patchF0 = f0_.boundaryField()[patchI];
        vectorField& patchS0 = s0_.boundaryFieldRef()[patchI];
        vectorField& patchN0 = n0_.boundaryFieldRef()[patchI];
        tensorField& patchR = R_.boundaryFieldRef()[patchI];

        forAll(patchF0, faceI)
        {
            constructBasis
            (
                patchF0[faceI],
                basisTolerance,
                patchS0[faceI],
                patchN0[faceI],
                patchR[faceI]
            );
            cellBasisStats.add
            (
                patchF0[faceI],
                patchS0[faceI],
                patchN0[faceI],
                patchR[faceI]
            );
        }
    }

#ifdef OPENFOAM_NOT_EXTEND
    vectorField& f0fInternal = f0f_.primitiveFieldRef();
    vectorField& s0fInternal = s0f_.primitiveFieldRef();
    vectorField& n0fInternal = n0f_.primitiveFieldRef();
    tensorField& RfInternal = Rf_.primitiveFieldRef();
#else
    vectorField& f0fInternal = f0f_.internalField();
    vectorField& s0fInternal = s0f_.internalField();
    vectorField& n0fInternal = n0f_.internalField();
    tensorField& RfInternal = Rf_.internalField();
#endif

    forAll(f0fInternal, faceI)
    {
        constructBasis
        (
            f0fInternal[faceI],
            basisTolerance,
            s0fInternal[faceI],
            n0fInternal[faceI],
            RfInternal[faceI]
        );
        faceBasisStats.add
        (
            f0fInternal[faceI],
            s0fInternal[faceI],
            n0fInternal[faceI],
            RfInternal[faceI]
        );
    }

    forAll(f0f_.boundaryField(), patchI)
    {
        const vectorField& patchF0f = f0f_.boundaryField()[patchI];
        vectorField& patchS0f = s0f_.boundaryFieldRef()[patchI];
        vectorField& patchN0f = n0f_.boundaryFieldRef()[patchI];
        tensorField& patchRf = Rf_.boundaryFieldRef()[patchI];

        forAll(patchF0f, faceI)
        {
            constructBasis
            (
                patchF0f[faceI],
                basisTolerance,
                patchS0f[faceI],
                patchN0f[faceI],
                patchRf[faceI]
            );
            faceBasisStats.add
            (
                patchF0f[faceI],
                patchS0f[faceI],
                patchN0f[faceI],
                patchRf[faceI]
            );
        }
    }

    Info<< "GuccioneElastic basis diagnostics" << nl
        << "    calculateStressInLocalCoordinateSystem = "
        << dict.lookupOrDefault<Switch>
           (
               "calculateStressInLocalCoordinateSystem",
               Switch(false)
           )
        << nl;
    printBasisStats("cell", cellBasisStats);
    printBasisStats("face", faceBasisStats);

    if (dict.lookupOrDefault<Switch>("writeS0N0R", Switch(false)))
    {
        Info<< "Writing s0, n0 and R" << endl;
        s0_.write();
        n0_.write();
        R_.write();
        s0f_.write();
        n0f_.write();
        Rf_.write();
    }

    if (apexRegularisationEnabled_)
    {
        makeApexRegularisationMasks();
    }

    if (pressureDisplacement_)
    {
        calcInitialShearModulus();
        calcEffectiveShearModulus();
    }

    if (apexRegularisationEnabled_)
    {
        reportApexRegularisationCoefficients();
    }

    if (dict.lookupOrDefault<Switch>("validateIsochoricSplit", Switch(false)))
    {
        validateIsochoricSplitImplementation();
    }

    if
    (
        dict.lookupOrDefault<Switch>
        (
            "validateApexRegularisationTangent",
            Switch(false)
        )
    )
    {
        validateApexRegularisationTangent();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GuccioneElastic::~GuccioneElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::impK() const
{
    if (pressureDisplacement_)
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
                impKcoeff_*muEff_
            )
        );
    }

    if (useApexImplicitStiffness())
    {
        const tmp<volScalarField> tMu = shearModulus();

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
                (4.0/3.0)*tMu() + bulkModulus_
            )
        );
    }

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
            (4.0/3.0)*mu_ + bulkModulus_
        )
    );
}


void Foam::GuccioneElastic::validateMaterialTangentColumns
(
    const List<mat66>& matTan,
    const word& context
) const
{
    label zeroColumnFaces[symmTensor::nComponents];
    scalar minColumnMag[symmTensor::nComponents];
    scalar maxColumnMag[symmTensor::nComponents];

    for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
    {
        zeroColumnFaces[cmptI] = 0;
        minColumnMag[cmptI] = VGREAT;
        maxColumnMag[cmptI] = 0.0;
    }

    label nonFiniteEntries = 0;

    forAll(matTan, faceI)
    {
        for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
        {
            scalar columnMag = 0.0;

            for (label rowI = 0; rowI < symmTensor::nComponents; rowI++)
            {
                const scalar value = matTan[faceI](rowI, cmptI);
                columnMag += mag(value);

                if (!std::isfinite(value))
                {
                    ++nonFiniteEntries;
                }
            }

            minColumnMag[cmptI] = min(minColumnMag[cmptI], columnMag);
            maxColumnMag[cmptI] = max(maxColumnMag[cmptI], columnMag);

            if (columnMag <= SMALL)
            {
                ++zeroColumnFaces[cmptI];
            }
        }
    }

    if (nonFiniteEntries)
    {
        FatalErrorInFunction
            << "Non-finite entries in GuccioneElastic material tangent "
            << "validation for " << context << nl
            << "    count = " << nonFiniteEntries
            << abort(FatalError);
    }

    for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
    {
        if (zeroColumnFaces[cmptI] == matTan.size())
        {
            FatalErrorInFunction
                << "GuccioneElastic material tangent column " << cmptI
                << " is zero on every face during " << context << nl
                << "This indicates that a previously assembled column may "
                << "have been cleared."
                << abort(FatalError);
        }
    }

    Info<< "GuccioneElastic material tangent validation (" << context << ")"
        << nl
        << "    faces checked = " << matTan.size() << nl
        << "    zero-column face counts = ";

    for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
    {
        Info<< (cmptI ? " " : "")
            << cmptI << ":" << zeroColumnFaces[cmptI];
    }

    Info<< nl
        << "    column magnitude min/max = ";

    for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
    {
        Info<< (cmptI ? " " : "")
            << cmptI << ":" << minColumnMag[cmptI]
            << "/" << maxColumnMag[cmptI];
    }

    Info<< nl;

    if (apexRegularisationMaskfPtr_.valid())
    {
        const surfaceScalarField& maskf = apexRegularisationMaskf();
        label baselineFaces = 0;
        label apexFaces = 0;
        label transitionFaces = 0;

#ifdef OPENFOAM_NOT_EXTEND
        const scalarField& maskfI = maskf.primitiveField();
#else
        const scalarField& maskfI = maskf.internalField();
#endif

        forAll(maskfI, faceI)
        {
            if (maskfI[faceI] <= SMALL)
            {
                ++baselineFaces;
            }
            else if (maskfI[faceI] >= 1.0 - SMALL)
            {
                ++apexFaces;
            }
            else
            {
                ++transitionFaces;
            }
        }

        forAll(maskf.boundaryField(), patchI)
        {
            const scalarField& patch = maskf.boundaryField()[patchI];

            forAll(patch, faceI)
            {
                if (patch[faceI] <= SMALL)
                {
                    ++baselineFaces;
                }
                else if (patch[faceI] >= 1.0 - SMALL)
                {
                    ++apexFaces;
                }
                else
                {
                    ++transitionFaces;
                }
            }
        }

        Info<< "    apex passive tangent face-mask counts "
            << "baseline/apex/transition = "
            << baselineFaces << " / " << apexFaces
            << " / " << transitionFaces << nl
            << "    calculateStress() supplied the passive coefficient "
            << "mask used by the finite-difference tangent" << nl;
    }
}


void Foam::GuccioneElastic::validateApexRegularisationTangent() const
{
    if (!dict().found("tangentEps"))
    {
        FatalErrorInFunction
            << "validateApexRegularisationTangent requires tangentEps "
            << "in the GuccioneElastic dictionary"
            << abort(FatalError);
    }

    const dimensionedTensor validationGradD
    (
        "validationGradD",
        dimless,
        tensor
        (
             0.071,  0.017, -0.011,
            -0.023,  0.043,  0.019,
             0.029, -0.013, -0.037
        )
    );

    surfaceTensorField gradDRef
    (
        IOobject
        (
            "apexRegularisationTangentValidationGradD",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        validationGradD
    );

    surfaceSymmTensorField sigmaRef
    (
        IOobject
        (
            "apexRegularisationTangentValidationSigma",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    const_cast<GuccioneElastic&>(*this).calculateStress(sigmaRef, gradDRef);

    surfaceSymmTensorField sigmaPerturb("sigmaPerturb", sigmaRef);
    surfaceTensorField gradDPerturb("gradDPerturb", gradDRef);
    const scalar eps(readScalar(dict().lookup("tangentEps")));

    List<mat66> matTan(mesh().nFaces());
    forAll(matTan, faceI)
    {
        matTan[faceI].clear();
    }

    const label XX = symmTensor::XX;
    const label YY = symmTensor::YY;
    const label ZZ = symmTensor::ZZ;
    const label XY = symmTensor::XY;
    const label YZ = symmTensor::YZ;
    const label XZ = symmTensor::XZ;

    for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
    {
        label tensorCmptI = -1;
        if (cmptI == symmTensor::XX)
        {
            tensorCmptI = tensor::XX;
        }
        else if (cmptI == symmTensor::XY)
        {
            tensorCmptI = tensor::XY;
        }
        else if (cmptI == symmTensor::XZ)
        {
            tensorCmptI = tensor::XZ;
        }
        else if (cmptI == symmTensor::YY)
        {
            tensorCmptI = tensor::YY;
        }
        else if (cmptI == symmTensor::YZ)
        {
            tensorCmptI = tensor::YZ;
        }
        else
        {
            tensorCmptI = tensor::ZZ;
        }

        gradDPerturb = 1.0*gradDRef;
        gradDPerturb.replace
        (
            tensorCmptI,
            gradDRef.component(tensorCmptI) + eps
        );

        const_cast<GuccioneElastic&>(*this).calculateStress
        (
            sigmaPerturb,
            gradDPerturb
        );

        const surfaceSymmTensorField tangCmpt((sigmaPerturb - sigmaRef)/eps);
        const symmTensorField& tangCmptI = tangCmpt.internalField();

        forAll(tangCmptI, faceI)
        {
            mat66& curMatTan = matTan[faceI];

            curMatTan(XX, cmptI) = tangCmptI[faceI][XX];
            curMatTan(YY, cmptI) = tangCmptI[faceI][YY];
            curMatTan(ZZ, cmptI) = tangCmptI[faceI][ZZ];
            curMatTan(XY, cmptI) = tangCmptI[faceI][XY];
            curMatTan(YZ, cmptI) = tangCmptI[faceI][YZ];
            curMatTan(XZ, cmptI) = tangCmptI[faceI][XZ];
        }

        forAll(tangCmpt.boundaryField(), patchI)
        {
            const symmTensorField& tangCmptP =
                tangCmpt.boundaryField()[patchI];
            const label start = mesh().boundaryMesh()[patchI].start();

            forAll(tangCmptP, fI)
            {
                const label faceID = start + fI;
                mat66& curMatTan = matTan[faceID];

                curMatTan(XX, cmptI) = tangCmptP[fI][XX];
                curMatTan(YY, cmptI) = tangCmptP[fI][YY];
                curMatTan(ZZ, cmptI) = tangCmptP[fI][ZZ];
                curMatTan(XY, cmptI) = tangCmptP[fI][XY];
                curMatTan(YZ, cmptI) = tangCmptP[fI][YZ];
                curMatTan(XZ, cmptI) = tangCmptP[fI][XZ];
            }
        }
    }

    validateMaterialTangentColumns
    (
        matTan,
        "synthetic apexRegularisation perturbation"
    );
}


void Foam::GuccioneElastic::materialTangentField(List<mat66>& matTan) const
{
    // Set the list size
    matTan.resize(mesh().nFaces());

    forAll(matTan, faceI)
    {
        matTan[faceI].clear();
    }

    // Calculate tangent field
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradDRef =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Lookup current stress and store it as the reference
        // const surfaceSymmTensorField& sigmaRef =
        //     mesh().lookupObject<surfaceSymmTensorField>("sigmaf")
        // Calculate sigmaRef to be consistent with gradDRef;
        surfaceSymmTensorField sigmaRef
        (
            "sigmaRef", 1.0*mesh().lookupObject<surfaceSymmTensorField>("sigmaf")
        );
        const_cast<GuccioneElastic&>(*this).calculateStress(sigmaRef, gradDRef);

        // Create fields to be used for perturbations
        surfaceSymmTensorField sigmaPerturb("sigmaPerturb", sigmaRef);
        surfaceTensorField gradDPerturb("gradDPerturb", gradDRef);

        // Small number used for perturbations
        const scalar eps(readScalar(dict().lookup("tangentEps")));

        // Define matrix indices for readability
        const label XX = symmTensor::XX;
        const label YY = symmTensor::YY;
        const label ZZ = symmTensor::ZZ;
        const label XY = symmTensor::XY;
        const label YZ = symmTensor::YZ;
        const label XZ = symmTensor::XZ;

        // For each component of gradD, sequentially apply a perturbation and
        // then calculate the resulting sigma
        for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
        {
            // Map tensor component to symmTensor
            // We can avoid this is we perturb epsilon directly
            label tensorCmptI = -1;
            if (cmptI == symmTensor::XX)
            {
                tensorCmptI = tensor::XX;
            }
            else if (cmptI == symmTensor::XY)
            {
                tensorCmptI = tensor::XY;
            }
            else if (cmptI == symmTensor::XZ)
            {
                tensorCmptI = tensor::XZ;
            }
            else if (cmptI == symmTensor::YY)
            {
                tensorCmptI = tensor::YY;
            }
            else if (cmptI == symmTensor::YZ)
            {
                tensorCmptI = tensor::YZ;
            }
            else // if (cmptI == symmTensor::ZZ)
            {
                tensorCmptI = tensor::ZZ;
            }

            // Reset gradDPerturb and multiply by 1.0 to avoid it being removed
            // from the object registry
            gradDPerturb = 1.0*gradDRef;

            // Perturb this component of gradD
            gradDPerturb.replace
            (
                tensorCmptI, gradDRef.component(tensorCmptI) + eps
            );

            // Calculate perturbed stress
            const_cast<GuccioneElastic&>(*this).calculateStress(sigmaPerturb, gradDPerturb);

            // Calculate tangent component
            const surfaceSymmTensorField tangCmpt((sigmaPerturb - sigmaRef)/eps);
            const symmTensorField& tangCmptI = tangCmpt.internalField();

            // Insert tangent component
            forAll(tangCmptI, faceI)
            {
                // Take a reference to the current tangent
                mat66& curMatTan = matTan[faceI];

                curMatTan(XX, cmptI) = tangCmptI[faceI][XX];
                curMatTan(YY, cmptI) = tangCmptI[faceI][YY];
                curMatTan(ZZ, cmptI) = tangCmptI[faceI][ZZ];
                curMatTan(XY, cmptI) = tangCmptI[faceI][XY];
                curMatTan(YZ, cmptI) = tangCmptI[faceI][YZ];
                curMatTan(XZ, cmptI) = tangCmptI[faceI][XZ];
            }

            forAll(tangCmpt.boundaryField(), patchI)
            {
                const symmTensorField& tangCmptP =
                    tangCmpt.boundaryField()[patchI];
                const label start = mesh().boundaryMesh()[patchI].start();

                forAll(tangCmptP, fI)
                {
                    const label faceID = start + fI;

                    // Take a reference to the current tangent
                    mat66& curMatTan = matTan[faceID];

                    curMatTan(XX, cmptI) = tangCmptP[fI][XX];
                    curMatTan(YY, cmptI) = tangCmptP[fI][YY];
                    curMatTan(ZZ, cmptI) = tangCmptP[fI][ZZ];
                    curMatTan(XY, cmptI) = tangCmptP[fI][XY];
                    curMatTan(YZ, cmptI) = tangCmptP[fI][YZ];
                    curMatTan(XZ, cmptI) = tangCmptP[fI][XZ];
                }
            }
        }
    }

    if
    (
        dict().lookupOrDefault<Switch>
        (
            "validateApexRegularisationTangent",
            Switch(false)
        )
    )
    {
        validateMaterialTangentColumns(matTan, "materialTangentField");
    }
}


Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulus",
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


Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::shearModulus() const
{
    tmp<volScalarField> tMu
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
            mu_
        )
    );

    if (useApexImplicitStiffness())
    {
        const dimensionedScalar muIso
        (
            "muIsoApex",
            mu_.dimensions(),
            0.5*k_.value()*apexIsotropicB_
        );

        tMu.ref() += apexRegularisationMask()*(muIso - mu_);
        tMu.ref().correctBoundaryConditions();
    }

    return tMu;
}


void Foam::GuccioneElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, bulkModulus_))
    {
        return;
    }

    if (pressureDisplacement_)
    {
        // Take a reference to the deformation gradient to make the code easier
        // to read
        const volTensorField& F = this->F();
        const volScalarField J("J", det(F));

        volTensorField Fwork("Fwork", F);
        if (useIsochoricSplit_)
        {
            Fwork = pow(J, -1.0/3.0)*F;
        }
        const volTensorField FworkT("FworkT", Fwork.T());

        // Calculate the right Cauchy-Green deformation tensor
        const volSymmTensorField C("C", symm(FworkT & Fwork));

        // Calculate the Green-Lagrange strain
        const volSymmTensorField E("E", 0.5*(C - I));

        const Switch useLocalCoordSys
        (
            dict().lookupOrDefault<Switch>
            (
                "calculateStressInLocalCoordinateSystem",
                Switch(false)
            )
        );

        if (useLocalCoordSys)
        {
            // Calculate the Green strain in the local coordinate system
            const volTensorField RT("RT", R_.T());
            const volSymmTensorField EStar("EStar", symm(RT & E & R_));

            // Extract the components of EStar
            // Note: EStar is symmetric
            const volScalarField E11("E11", EStar.component(symmTensor::XX));
            const volScalarField E12("E12", EStar.component(symmTensor::XY));
            const volScalarField E13("E13", EStar.component(symmTensor::XZ));
            const volScalarField E22("E22", EStar.component(symmTensor::YY));
            const volScalarField E23("E23", EStar.component(symmTensor::YZ));
            const volScalarField E33("E33", EStar.component(symmTensor::ZZ));

            if (useApexPassiveIsotropisation())
            {
                volScalarField cf
                (
                    IOobject
                    (
                        "cfApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cf", dimless, cf_)
                );
                volScalarField ct
                (
                    IOobject
                    (
                        "ctApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("ct", dimless, ct_)
                );
                volScalarField cfs
                (
                    IOobject
                    (
                        "cfsApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cfs", dimless, cfs_)
                );

                cf += apexRegularisationMask()*(apexIsotropicB_ - cf_);
                ct += apexRegularisationMask()*(apexIsotropicB_ - ct_);
                cfs += apexRegularisationMask()*(apexIsotropicB_ - cfs_);

                // Calculate Q
                const volScalarField Q
                (
                    "Q",
                    cf*sqr(E11)
                  + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
                  + cfs*(2*sqr(E12) + 2*sqr(E13))
                );

                // Calculate the derivative of Q wrt to EStar
                volSymmTensorField dQdEStar
                (
                    IOobject
                    (
                        "dQdEStar",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedSymmTensor("0", dimless, symmTensor::zero)
                );

                dQdEStar.replace(symmTensor::XX, 2*cf*E11);
                dQdEStar.replace(symmTensor::XY, 2*cfs*E12);
                dQdEStar.replace(symmTensor::XZ, 2*cfs*E13);
                dQdEStar.replace(symmTensor::YY, 2*ct*E22);
                dQdEStar.replace(symmTensor::YZ, 2*ct*E23);
                dQdEStar.replace(symmTensor::ZZ, 2*ct*E33);

                // Calculate the local 2nd Piola-Kirchhoff stress (without the
                // hydrostatic term)
                S_ = dQdEStar*0.5*k_*exp(Q);
            }
            else
            {
                // Calculate Q
                const volScalarField Q
                (
                    "Q",
                    cf_*sqr(E11)
                  + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
                  + cfs_*(2*sqr(E12) + 2*sqr(E13))
                );

                // Calculate the derivative of Q wrt to EStar
                volSymmTensorField dQdEStar
                (
                    IOobject
                    (
                        "dQdEStar",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedSymmTensor("0", dimless, symmTensor::zero)
                );

                dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
                dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
                dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
                dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
                dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
                dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

                // Calculate the local 2nd Piola-Kirchhoff stress (without the
                // hydrostatic term)
                S_ = dQdEStar*0.5*k_*exp(Q);
            }

            // Rotate S from the local fibre coordinate system to the global
            // coordinate system
            S_ = symm(R_ & S_ & RT);
        }
        else
        {
            // Calculate E . E
            const volSymmTensorField sqrE("sqrE", symm(E & E));

            // Calculate the invariants of E
            const volScalarField I1("I1", tr(E));
            const volScalarField I2("I2", 0.5*(sqr(tr(E)) - tr(sqrE)));
            const volScalarField I4("I4", E && f0f0_);
            const volScalarField I5("I5", sqrE && f0f0_);

            if (useApexPassiveIsotropisation())
            {
                volScalarField cf
                (
                    IOobject
                    (
                        "cfApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cf", dimless, cf_)
                );
                volScalarField ct
                (
                    IOobject
                    (
                        "ctApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("ct", dimless, ct_)
                );
                volScalarField cfs
                (
                    IOobject
                    (
                        "cfsApex",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cfs", dimless, cfs_)
                );

                cf += apexRegularisationMask()*(apexIsotropicB_ - cf_);
                ct += apexRegularisationMask()*(apexIsotropicB_ - ct_);
                cfs += apexRegularisationMask()*(apexIsotropicB_ - cfs_);

                // Calculate Q
                const volScalarField Q
                (
                    "Q",
                    ct*sqr(I1)
                  - 2.0*ct*I2
                 + (cf - 2.0*cfs + ct)*sqr(I4)
                 + 2.0*(cfs - ct)*I5
                );

                // Calculate the derivative of Q wrt to E
                const volSymmTensorField dQdE
                (
                    2.0*ct*E
                  + 2.0*(cf - 2.0*cfs + ct)*I4*f0f0_
                  + 2.0*(cfs - ct)*symm((E & f0f0_) + (f0f0_ & E))
                );

                // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic
                // term)
                S_ = dQdE*0.5*k_*exp(Q);
            }
            else
            {
                // Calculate Q
                const volScalarField Q
                (
                    "Q",
                    ct_*sqr(I1)
                  - 2.0*ct_*I2
                 + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
                 + 2.0*(cfs_ - ct_)*I5
                );

                // Calculate the derivative of Q wrt to E
                const volSymmTensorField dQdE
                (
                    2.0*ct_*E
                  + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0_
                  + 2.0*(cfs_ - ct_)*symm((E & f0f0_) + (f0f0_ & E))
                );

                // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic
                // term)
                S_ = dQdE*0.5*k_*exp(Q);
            }
        }

        // Convert the second Piola-Kirchhoff stress to the deviatoric Cauchy
        // stress
        const volSymmTensorField s("s", dev(symm(Fwork & S_ & FworkT))/J);

        // Lookup pressure field
        // During solid-model construction p may not be registered yet.
        // The mixed solid model applies the pressure split once it is available.
        if (!mesh().foundObject<volScalarField>("p"))
        {
            sigma = s;
            return;
        }

        // Add the pressure-displacement hydrostatic term
        const volScalarField& p =
            mesh().lookupObject<volScalarField>("p");
        sigma = s - p*I;
        return;
    }

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    volTensorField Fwork(F);
    if (useIsochoricSplit_)
    {
        Fwork = pow(J, -1.0/3.0)*F;
    }
    const volTensorField FworkT(Fwork.T());

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField C(symm(FworkT & Fwork));

    // Calculate the Green-Lagrange strain
    const volSymmTensorField E(0.5*(C - I));

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const volTensorField RT(R_.T());
        const volSymmTensorField EStar("EStar", symm(RT & E & R_));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const volScalarField E11("E11", EStar.component(symmTensor::XX));
        const volScalarField E12("E12", EStar.component(symmTensor::XY));
        const volScalarField E13("E13", EStar.component(symmTensor::XZ));
        const volScalarField E22("E22", EStar.component(symmTensor::YY));
        const volScalarField E23("E23", EStar.component(symmTensor::YZ));
        const volScalarField E33("E33", EStar.component(symmTensor::ZZ));

        if (useApexPassiveIsotropisation())
        {
            volScalarField cf
            (
                IOobject
                (
                    "cfApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cf", dimless, cf_)
            );
            volScalarField ct
            (
                IOobject
                (
                    "ctApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("ct", dimless, ct_)
            );
            volScalarField cfs
            (
                IOobject
                (
                    "cfsApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cfs", dimless, cfs_)
            );

            cf += apexRegularisationMask()*(apexIsotropicB_ - cf_);
            ct += apexRegularisationMask()*(apexIsotropicB_ - ct_);
            cfs += apexRegularisationMask()*(apexIsotropicB_ - cfs_);

            // Calculate Q
            const volScalarField Q
            (
                "Q",
                cf*sqr(E11)
              + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
              + cfs*(2*sqr(E12) + 2*sqr(E13))
            );

            // Calculate the derivative of Q wrt to EStar
            volSymmTensorField dQdEStar
            (
                IOobject
                (
                    "dQdEStar",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedSymmTensor("0", dimless, symmTensor::zero)
            );

            dQdEStar.replace(symmTensor::XX, 2*cf*E11);
            dQdEStar.replace(symmTensor::XY, 2*cfs*E12);
            dQdEStar.replace(symmTensor::XZ, 2*cfs*E13);
            dQdEStar.replace(symmTensor::YY, 2*ct*E22);
            dQdEStar.replace(symmTensor::YZ, 2*ct*E23);
            dQdEStar.replace(symmTensor::ZZ, 2*ct*E33);

            // Calculate the local 2nd Piola-Kirchhoff stress (without the
            // hydrostatic term)
            S_ = dQdEStar*0.5*k_*exp(Q);
        }
        else
        {
            // Calculate Q
            const volScalarField Q
            (
                "Q",
                cf_*sqr(E11)
              + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
              + cfs_*(2*sqr(E12) + 2*sqr(E13))
            );

            // Calculate the derivative of Q wrt to EStar
            volSymmTensorField dQdEStar
            (
                IOobject
                (
                    "dQdEStar",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedSymmTensor("0", dimless, symmTensor::zero)
            );

            dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
            dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
            dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
            dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
            dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
            dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

            // Calculate the local 2nd Piola-Kirchhoff stress (without the
            // hydrostatic term)
            S_ = dQdEStar*0.5*k_*exp(Q);
        }

        // Rotate S from the local fibre coordinate system to the global
        // coordinate system
        S_ = symm(R_ & S_ & RT);
    }
    else
    {
        // Calculate E . E
        const volSymmTensorField sqrE(symm(E & E));

        // Calculate the invariants of E
        const volScalarField I1(tr(E));
        const volScalarField I2(0.5*(sqr(tr(E)) - tr(sqrE)));
        const volScalarField I4(E && f0f0_);
        const volScalarField I5(sqrE && f0f0_);

        if (useApexPassiveIsotropisation())
        {
            volScalarField cf
            (
                IOobject
                (
                    "cfApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cf", dimless, cf_)
            );
            volScalarField ct
            (
                IOobject
                (
                    "ctApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("ct", dimless, ct_)
            );
            volScalarField cfs
            (
                IOobject
                (
                    "cfsApex",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("cfs", dimless, cfs_)
            );

            cf += apexRegularisationMask()*(apexIsotropicB_ - cf_);
            ct += apexRegularisationMask()*(apexIsotropicB_ - ct_);
            cfs += apexRegularisationMask()*(apexIsotropicB_ - cfs_);

            // Calculate Q
            const volScalarField Q
            (
                ct*sqr(I1)
              - 2.0*ct*I2
             + (cf - 2.0*cfs + ct)*sqr(I4)
             + 2.0*(cfs - ct)*I5
            );

            // Calculate the derivative of Q wrt to E
            const volSymmTensorField dQdE
            (
                2.0*ct*E
              + 2.0*(cf - 2.0*cfs + ct)*I4*f0f0_
              + 2.0*(cfs - ct)*symm((E & f0f0_) + (f0f0_ & E))
            );

            // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic term)
            S_ = dQdE*0.5*k_*exp(Q);
        }
        else
        {
            // Calculate Q
            const volScalarField Q
            (
                ct_*sqr(I1)
              - 2.0*ct_*I2
             + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
             + 2.0*(cfs_ - ct_)*I5
            );

            // Calculate the derivative of Q wrt to E
            const volSymmTensorField dQdE
            (
                2.0*ct_*E
              + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0_
              + 2.0*(cfs_ - ct_)*symm((E & f0f0_) + (f0f0_ & E))
            );

            // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic term)
            S_ = dQdE*0.5*k_*exp(Q);
        }
    }

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and take
    // the deviatoric component
    // s = dev(Fwork & S & Fwork.T)/J, where Fwork is either F or Fbar
    const volSymmTensorField s(dev(symm(Fwork & S_ & FworkT))/J);

    // Calculate the hydrostatic stress
    updateSigmaHyd
    (
        0.5*bulkModulus_*(pow(J, 2.0) - 1.0)/J,
        (4.0/3.0)*mu_ + bulkModulus_
    );

    // Convert the second Piola-Kirchhoff deviatoric stress to the Cauchy stress
    // and add hydrostatic stress term
    sigma = s + sigmaHyd()*I;
}


void Foam::GuccioneElastic::correct(surfaceSymmTensorField& sigma)
{
    if (pressureDisplacement_)
    {
        // Update the deformation gradient field
        // Note: if true is returned, it means that linearised elasticity was
        // enforced by the solver via the enforceLinear switch
        if (updateF(sigma, mu_, bulkModulus_))
        {
            return;
        }

        expQf_.storePrevIter();

        // Take a reference to the deformation gradient to make the code easier
        // to read
        const surfaceTensorField& Ff = this->Ff();
        const surfaceScalarField Jf("Jf", det(Ff));

        surfaceTensorField FfWork("FfWork", Ff);
        if (useIsochoricSplit_)
        {
            FfWork = pow(Jf, -1.0/3.0)*Ff;
        }
        const surfaceTensorField FfWorkT("FfWorkT", FfWork.T());

        // Calculate the right Cauchy-Green deformation tensor
        const surfaceSymmTensorField C("C", symm(FfWorkT & FfWork));

        // Calculate the Green-Lagrange strain
        const surfaceSymmTensorField E("E", 0.5*(C - I));

        const Switch useLocalCoordSys
        (
            dict().lookupOrDefault<Switch>
            (
                "calculateStressInLocalCoordinateSystem",
                Switch(false)
            )
        );

        if (useLocalCoordSys)
        {
            // Calculate the Green strain in the local coordinate system
            const surfaceTensorField RfT("RfT", Rf_.T());
            const surfaceSymmTensorField EStar("EStar", symm(RfT & E & Rf_));

            // Extract the components of EStar
            // Note: EStar is symmetric
            const surfaceScalarField E11
            (
                "E11", EStar.component(symmTensor::XX)
            );
            const surfaceScalarField E12
            (
                "E12", EStar.component(symmTensor::XY)
            );
            const surfaceScalarField E13
            (
                "E13", EStar.component(symmTensor::XZ)
            );
            const surfaceScalarField E22
            (
                "E22", EStar.component(symmTensor::YY)
            );
            const surfaceScalarField E23
            (
                "E23", EStar.component(symmTensor::YZ)
            );
            const surfaceScalarField E33
            (
                "E33", EStar.component(symmTensor::ZZ)
            );

            if (useApexPassiveIsotropisation())
            {
                surfaceScalarField cf
                (
                    IOobject
                    (
                        "cfApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cf", dimless, cf_)
                );
                surfaceScalarField ct
                (
                    IOobject
                    (
                        "ctApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("ct", dimless, ct_)
                );
                surfaceScalarField cfs
                (
                    IOobject
                    (
                        "cfsApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cfs", dimless, cfs_)
                );

                cf += apexRegularisationMaskf()*(apexIsotropicB_ - cf_);
                ct += apexRegularisationMaskf()*(apexIsotropicB_ - ct_);
                cfs += apexRegularisationMaskf()*(apexIsotropicB_ - cfs_);

                // Calculate Q
                const surfaceScalarField Q
                (
                    "Q",
                    cf*sqr(E11)
                  + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
                  + cfs*(2*sqr(E12) + 2*sqr(E13))
                );

                // Calculate the derivative of Q wrt to EStar
                surfaceSymmTensorField dQdEStar
                (
                    IOobject
                    (
                        "dQdEStar",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedSymmTensor("0", dimless, symmTensor::zero)
                );

                dQdEStar.replace(symmTensor::XX, 2*cf*E11);
                dQdEStar.replace(symmTensor::XY, 2*cfs*E12);
                dQdEStar.replace(symmTensor::XZ, 2*cfs*E13);
                dQdEStar.replace(symmTensor::YY, 2*ct*E22);
                dQdEStar.replace(symmTensor::YZ, 2*ct*E23);
                dQdEStar.replace(symmTensor::ZZ, 2*ct*E33);

                expQf_ = exp(Q);
                expQf_.relax();

                // Calculate the local 2nd Piola-Kirchhoff stress (without the
                // hydrostatic term)
                Sf_ = dQdEStar*0.5*k_*expQf_;
            }
            else
            {
                // Calculate Q
                const surfaceScalarField Q
                (
                    "Q",
                    cf_*sqr(E11)
                  + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
                  + cfs_*(2*sqr(E12) + 2*sqr(E13))
                );

                // Calculate the derivative of Q wrt to EStar
                surfaceSymmTensorField dQdEStar
                (
                    IOobject
                    (
                        "dQdEStar",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedSymmTensor("0", dimless, symmTensor::zero)
                );

                dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
                dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
                dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
                dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
                dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
                dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

                expQf_ = exp(Q);
                expQf_.relax();

                // Calculate the local 2nd Piola-Kirchhoff stress (without the
                // hydrostatic term)
                Sf_ = dQdEStar*0.5*k_*expQf_;
            }

            // Rotate S from the local fibre coordinate system to the global
            // coordinate system
            Sf_ = symm(Rf_ & Sf_ & RfT);
        }
        else
        {
            // Calculate E . E
            const surfaceSymmTensorField sqrE("sqrE", symm(E & E));

            // Calculate the invariants of E
            const surfaceScalarField I1("I1", tr(E));
            const surfaceScalarField I2
            (
                "I2",
                0.5*(sqr(tr(E)) - tr(sqrE))
            );
            const surfaceScalarField I4("I4", E && f0f0f_);
            const surfaceScalarField I5("I5", sqrE && f0f0f_);

            if (useApexPassiveIsotropisation())
            {
                surfaceScalarField cf
                (
                    IOobject
                    (
                        "cfApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cf", dimless, cf_)
                );
                surfaceScalarField ct
                (
                    IOobject
                    (
                        "ctApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("ct", dimless, ct_)
                );
                surfaceScalarField cfs
                (
                    IOobject
                    (
                        "cfsApexf",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("cfs", dimless, cfs_)
                );

                cf += apexRegularisationMaskf()*(apexIsotropicB_ - cf_);
                ct += apexRegularisationMaskf()*(apexIsotropicB_ - ct_);
                cfs += apexRegularisationMaskf()*(apexIsotropicB_ - cfs_);

                // Calculate Q
                const surfaceScalarField Q
                (
                    "Q",
                    ct*sqr(I1)
                  - 2.0*ct*I2
                  + (cf - 2.0*cfs + ct)*sqr(I4)
                  + 2.0*(cfs - ct)*I5
                );

                // Calculate the derivative of Q wrt to E
                const surfaceSymmTensorField dQdE
                (
                    2.0*ct*E
                  + 2.0*(cf - 2.0*cfs + ct)*I4*f0f0f_
                  + 2.0*(cfs - ct)*symm((E & f0f0f_) + (f0f0f_ & E))
                );

                expQf_ = exp(Q);
                expQf_.relax();

                // Update the 2nd Piola-Kirchhoff stress (without the
                // hydrostatic term)
                Sf_ = dQdE*0.5*k_*expQf_;
            }
            else
            {
                // Calculate Q
                const surfaceScalarField Q
                (
                    "Q",
                    ct_*sqr(I1)
                  - 2.0*ct_*I2
                  + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
                  + 2.0*(cfs_ - ct_)*I5
                );

                // Calculate the derivative of Q wrt to E
                const surfaceSymmTensorField dQdE
                (
                    2.0*ct_*E
                  + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0f_
                  + 2.0*(cfs_ - ct_)*symm((E & f0f0f_) + (f0f0f_ & E))
                );

                expQf_ = exp(Q);
                expQf_.relax();

                // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic
                // term)
                Sf_ = dQdE*0.5*k_*expQf_;
            }
        }

        // Convert the second Piola-Kirchhoff stress to the deviatoric Cauchy
        // stress
        const surfaceSymmTensorField sf
        (
            "sf",
            dev(symm(FfWork & Sf_ & FfWorkT))/Jf
        );

        // Lookup pressure field
        // During solid-model construction pf may not be registered yet.
        // The mixed solid model applies the pressure split once it is available.
        if (!mesh().foundObject<surfaceScalarField>("pf"))
        {
            sigma = sf;
            return;
        }

        // Add the pressure-displacement hydrostatic term
        const surfaceScalarField& pf =
            mesh().lookupObject<surfaceScalarField>("pf");
        sigma = sf - pf*I;
        return;
    }

    const surfaceTensorField gradD("gradDFromFf", this->Ff().T() - I);
    calculateStress(sigma, gradD);
}


void Foam::GuccioneElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}


void Foam::GuccioneElastic::calcDevCauchy
(
    const tensor& F,
    const symmTensor& f0f0,
    const tensor& R,
    symmTensor& devSigma
) const
{
    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    calcDevCauchy
    (
        F,
        f0f0,
        R,
        useIsochoricSplit_,
        useLocalCoordSys,
        cf_,
        ct_,
        cfs_,
        devSigma
    );
}


void Foam::GuccioneElastic::calcDevCauchy
(
    const tensor& F,
    const symmTensor& f0f0,
    const tensor& R,
    const Switch& useIsochoricSplit,
    const Switch& useLocalCoordSys,
    const scalar cf,
    const scalar ct,
    const scalar cfs,
    symmTensor& devSigma
) const
{
    const scalar J = det(F);
    const tensor Fwork(useIsochoricSplit ? Foam::pow(J, -1.0/3.0)*F : F);

    // Calculate the right Cauchy-Green deformation tensor
    const tensor FT(Fwork.T());
    const symmTensor C(symm(FT & Fwork));

    // Calculate the Green-Lagrange strain
    const symmTensor E(0.5*(C - I));

    symmTensor S = symmTensor::zero;

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const tensor RT(R.T());
        const symmTensor EStar(symm(RT & E & R));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const scalar E11(EStar.xx());
        const scalar E12(EStar.xy());
        const scalar E13(EStar.xz());
        const scalar E22(EStar.yy());
        const scalar E23(EStar.yz());
        const scalar E33(EStar.zz());

        // Calculate Q
        const scalar Q
        (
            cf*sqr(E11)
          + ct*(sqr(E22) + sqr(E33) + 2*sqr(E23))
          + cfs*(2*sqr(E12) + 2*sqr(E13))
        );

        // Calculate the derivative of Q wrt to EStar
        symmTensor dQdEStar = symmTensor::zero;

        dQdEStar.xx() = 2*cf*E11;
        dQdEStar.xy() = 2*cfs*E12;
        dQdEStar.xz() = 2*cfs*E13;
        dQdEStar.yy() = 2*ct*E22;
        dQdEStar.yz() = 2*ct*E23;
        dQdEStar.zz() = 2*ct*E33;

        // Calculate the local 2nd Piola-Kirchhoff stress
        // (without the hydrostatic term)
        S = dQdEStar*0.5*k_.value()*exp(Q);

        // Rotate S from the local fibre coordinate system
        // to the global coordinate system
        S = symm(R & S & RT);
    }
    else
    {
        // Calculate E . E
        const symmTensor sqrE(symm(E & E));

        // Calculate the invariants of E
        const scalar I1(tr(E));
        const scalar I2(0.5*(sqr(tr(E)) - tr(sqrE)));
        const scalar I4(E && f0f0);
        const scalar I5(sqrE && f0f0);

        // Calculate Q
        const scalar Q
        (
            ct*sqr(I1)
          - 2.0*ct*I2
          + (cf - 2.0*cfs + ct)*sqr(I4)
          + 2.0*(cfs - ct)*I5
        );

        // Calculate the derivative of Q wrt to E
        const symmTensor dQdE
        (
            2.0*ct*E
          + 2.0*(cf - 2.0*cfs + ct)*I4*f0f0
          + 2.0*(cfs - ct)*symm((E & f0f0) + (f0f0 & E))
        );

        // Calculate the 2nd Piola-Kirchhoff stress
        // (without the hydrostatic term)
        S = dQdE*0.5*k_.value()*::exp(Q);
    }

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and take
    // the deviatoric component
    devSigma = dev(symm(Fwork & S & FT))/J;
}


void Foam::GuccioneElastic::validateIsochoricSplitImplementation() const
{
    const symmTensor f0f0(f0f0_[0]);
    const tensor R(R_[0]);

    const scalar tolerance = 1e-10;

    symmTensor sigmaLegacy = symmTensor::zero;
    symmTensor sigmaIso = symmTensor::zero;

    calcDevCauchy(tensor::I, f0f0, R, false, false, cf_, ct_, cfs_, sigmaLegacy);
    calcDevCauchy(tensor::I, f0f0, R, true, false, cf_, ct_, cfs_, sigmaIso);

    const scalar identityError =
        max(mag(sigmaLegacy), mag(sigmaIso));

    tensor Fiso
    (
        1.2, 0.1, 0.0,
        0.0, 0.9, 0.2,
        0.0, 0.0, 1.0/(1.2*0.9)
    );

    calcDevCauchy(Fiso, f0f0, R, false, false, cf_, ct_, cfs_, sigmaLegacy);
    calcDevCauchy(Fiso, f0f0, R, true, false, cf_, ct_, cfs_, sigmaIso);

    const scalar isoAgreement = mag(sigmaLegacy - sigmaIso);

    const tensor Fdilation(1.15*tensor::I);
    calcDevCauchy(Fdilation, f0f0, R, true, false, cf_, ct_, cfs_, sigmaIso);
    const scalar dilationError = mag(sigmaIso);

    const tensor Fgeneral
    (
        1.12, 0.08, 0.03,
        0.04, 0.95, 0.11,
        0.02, 0.05, 1.18
    );
    const scalar Jgeneral = det(Fgeneral);
    const tensor Fbar(Foam::pow(Jgeneral, -1.0/3.0)*Fgeneral);
    const scalar FbarDetError = mag(det(Fbar) - 1.0);

    symmTensor sigmaIsoInvariant = symmTensor::zero;
    symmTensor sigmaIsoLocal = symmTensor::zero;
    calcDevCauchy
    (
        Fgeneral, f0f0, R, true, false, cf_, ct_, cfs_, sigmaIsoInvariant
    );
    calcDevCauchy
    (
        Fgeneral, f0f0, R, true, true, cf_, ct_, cfs_, sigmaIsoLocal
    );
    const scalar localInvariantDifference =
        mag(sigmaIsoInvariant - sigmaIsoLocal);

    Info<< "GuccioneElastic isochoric split validation" << nl
        << "    identity passive stress error = " << identityError << nl
        << "    det(F)=1 legacy/isochoric stress difference = "
        << isoAgreement << nl
        << "    pure dilation isochoric stress error = "
        << dilationError << nl
        << "    det(Fbar)-1 error = " << FbarDetError << nl
        << "    local/invariant isochoric stress difference = "
        << localInvariantDifference << endl;

    if
    (
        identityError > tolerance
     || isoAgreement > tolerance
     || dilationError > tolerance
     || FbarDetError > tolerance
    )
    {
        FatalErrorInFunction
            << "GuccioneElastic isochoric split validation failed" << nl
            << "    tolerance = " << tolerance << abort(FatalError);
    }
}


void Foam::GuccioneElastic::calcInitialShearModulus()
{
    scalarField mu(3, 0);
    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xy() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0, 1, 0) & devSigma);

        mu[0] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0, 0, 1) & devSigma);

        mu[1] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.yz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0, 0, 1) & devSigma);

        mu[2] = mag(dTraction.y()/gamma);
    }

    Info<< "m12 = " << mu[0] << endl;
    Info<< "m13 = " << mu[1] << endl;
    Info<< "m23 = " << mu[2] << endl;
    Info<< "Current mu = " << mu_.value() << endl;

    mu_.value() = max(mu);
}


void Foam::GuccioneElastic::calcEffectiveShearModulus()
{
    // This function updates shear modulus muEff_ at the end of time step,
    // using current state of deformation

    const volTensorField& F = this->F();
    const tensorField& FI = F.internalField();
    const tensorField FIT(FI.T()());

    const symmTensorField& f0f0I = f0f0_.internalField();
    const tensorField& RI = R_.internalField();

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& muEffI = muEff_.primitiveFieldRef();
#else
    scalarField& muEffI = muEff_.internalField();
#endif

    forAll(muEffI, cellI)
    {
        scalarField mu(3, 0);

        const scalar cf = implicitCoefficient(cf_, cellI);
        const scalar ct = implicitCoefficient(ct_, cellI);
        const scalar cfs = implicitCoefficient(cfs_, cellI);

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy
        (
            FI[cellI],
            f0f0I[cellI],
            RI[cellI],
            useIsochoricSplit_,
            useLocalCoordSys,
            cf,
            ct,
            cfs,
            devSigma
        );

        tensor pVectors = tensor::zero;
        vector pValues = vector::zero;
        eig3().eigen_decomposition(devSigma, pVectors, pValues);

        tensor pVectorsT = pVectors.T();
        devSigma = symm(pVectors & devSigma & pVectorsT);

        // mu12
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xy() = 0.001;
            pertGradDD.yx() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & F.oldTime()[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 1, 0) & (devSigmaNew - devSigma));

            mu[0] = mag(dTraction.x()/pertGradDD.xy());
        }

        // mu13
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xz() = 0.001;
            pertGradDD.zx() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & F.oldTime()[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 0, 1) & (devSigmaNew - devSigma));

            mu[1] = mag(dTraction.x()/pertGradDD.xz());
        }

        // mu23
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.yz() = 0.001;
            pertGradDD.zy() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & F.oldTime()[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 0, 1) & (devSigmaNew - devSigma));

            mu[2] = mag(dTraction.y()/pertGradDD.yz());
        }

        // mu12+
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xy() = 0.001;
            pertGradDD.yx() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & FI[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 1, 0) & (devSigmaNew - devSigma));

            mu[0] = max(mu[0], mag(dTraction.x()/pertGradDD.xy()));
        }

        // mu13+
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xz() = 0.001;
            pertGradDD.zx() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & FI[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 0, 1) & (devSigmaNew - devSigma));

            mu[1] = max(mu[1], mag(dTraction.x()/pertGradDD.xz()));
        }

        // mu23+
        {
            // Displacement gradient perturbation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.yz() = 0.001;
            pertGradDD.zy() = 0.001;

            // Transform from local to global coordinate system
            pertGradDD = pVectorsT & pertGradDD & pVectors;

            // Calculate new F
            tensor Fnew = (I + pertGradDD.T()) & FI[cellI];

            symmTensor devSigmaNew = symmTensor::zero;
            calcDevCauchy
            (
                Fnew,
                f0f0I[cellI],
                RI[cellI],
                useIsochoricSplit_,
                useLocalCoordSys,
                cf,
                ct,
                cfs,
                devSigmaNew
            );

            // Transform to principal coordinate system
            devSigmaNew = symm(pVectors & devSigmaNew & pVectorsT);

            vector dTraction =
                (vector(0, 0, 1) & (devSigmaNew - devSigma));

            mu[2] = max(mu[2], mag(dTraction.y()/pertGradDD.yz()));
        }

        muEffI[cellI] = max(mu);
    }

    muEff_.correctBoundaryConditions();
}


void Foam::GuccioneElastic::updateTotalFields()
{
    if (pressureDisplacement_)
    {
        calcEffectiveShearModulus();
    }
}

// ************************************************************************* //
