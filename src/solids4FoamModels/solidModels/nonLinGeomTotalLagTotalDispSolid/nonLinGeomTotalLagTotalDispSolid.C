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

#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "arosticaSpringDashpotTractionFvPatchVectorField.H"
#include "arosticaNormalSpringDashpotTractionFvPatchVectorField.H"
#include "arosticaVectorSpringDashpotTractionFvPatchVectorField.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "processorFvPatch.H"
#include "calculatedFvPatchFields.H"
#include "compatibilityFunctions.H"
#include "electroMechanicalLaw.H"
#include "primitiveMeshTools.H"
#include "globalIndex.H"
#include "pointFields.H"
#include "PstreamBuffers.H"
#ifdef USE_PETSC
    #include "petscErrorHandling.H"
#endif

#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagTotalDispSolid, dictionary
);


namespace
{

struct ScalarDiagnosticStats
{
    scalar min;
    scalar max;
    scalar mean;
    scalar rms;
    scalar maxAbs;

    ScalarDiagnosticStats()
    :
        min(VGREAT),
        max(-VGREAT),
        mean(0),
        rms(0),
        maxAbs(0)
    {}
};


ScalarDiagnosticStats scalarStats
(
    const volScalarField& field,
    const scalarField& V
)
{
    ScalarDiagnosticStats stats;

    scalar volume = 0;
    scalar volumeValue = 0;
    scalar volumeValueSqr = 0;

    forAll(field, cellI)
    {
        const scalar value = field[cellI];
        const scalar Vc = V[cellI];

        stats.min = min(stats.min, value);
        stats.max = max(stats.max, value);
        stats.maxAbs = max(stats.maxAbs, mag(value));
        volume += Vc;
        volumeValue += Vc*value;
        volumeValueSqr += Vc*sqr(value);
    }

    reduce(stats.min, minOp<scalar>());
    reduce(stats.max, maxOp<scalar>());
    reduce(stats.maxAbs, maxOp<scalar>());
    reduce(volume, sumOp<scalar>());
    reduce(volumeValue, sumOp<scalar>());
    reduce(volumeValueSqr, sumOp<scalar>());

    stats.mean = volumeValue/(volume + VSMALL);
    stats.rms = sqrt(volumeValueSqr/(volume + VSMALL));

    return stats;
}


scalar vectorRms
(
    const volVectorField& field,
    const scalarField& V
)
{
    scalar volume = 0;
    scalar volumeValueSqr = 0;

    forAll(field, cellI)
    {
        const scalar Vc = V[cellI];
        volume += Vc;
        volumeValueSqr += Vc*magSqr(field[cellI]);
    }

    reduce(volume, sumOp<scalar>());
    reduce(volumeValueSqr, sumOp<scalar>());

    return sqrt(volumeValueSqr/(volume + VSMALL));
}


scalar vectorSumRms
(
    const volVectorField& a,
    const volVectorField& b,
    const scalarField& V
)
{
    scalar volume = 0;
    scalar volumeValueSqr = 0;

    forAll(a, cellI)
    {
        const scalar Vc = V[cellI];
        volume += Vc;
        volumeValueSqr += Vc*magSqr(a[cellI] + b[cellI]);
    }

    reduce(volume, sumOp<scalar>());
    reduce(volumeValueSqr, sumOp<scalar>());

    return sqrt(volumeValueSqr/(volume + VSMALL));
}


void printScalarStats
(
    const word& name,
    const ScalarDiagnosticStats& stats
)
{
    Info<< "    " << name << ":" << nl
        << "        min                  = " << stats.min << nl
        << "        max                  = " << stats.max << nl
        << "        volume-weighted mean = " << stats.mean << nl
        << "        volume-weighted RMS  = " << stats.rms << nl
        << "        max(abs)             = " << stats.maxAbs << nl;
}


bool hasOptionalSubDict(const dictionary& dict, const word& subDictName)
{
    return dict.found(subDictName);
}


Switch optionalSubDictSwitch
(
    const dictionary& dict,
    const word& subDictName,
    const word& keyword,
    const bool defaultValue
)
{
    if (!hasOptionalSubDict(dict, subDictName))
    {
        return Switch(defaultValue);
    }

    return
        dict.subDict(subDictName).lookupOrDefault<Switch>
        (
            keyword,
            Switch(defaultValue)
        );
}


bool isFixedDisplacementBoundary(const fvPatchVectorField& patchField)
{
    return
        isA<fixedDisplacementFvPatchVectorField>(patchField)
     || isA<fixedDisplacementZeroShearFvPatchVectorField>(patchField);
}


bool isTractionBoundary(const fvPatchVectorField& patchField)
{
    return isA<solidTractionFvPatchVectorField>(patchField);
}


wordList diagnosticVectorPatchTypes
(
    const fvMesh& mesh,
    const volVectorField& referenceField
)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    forAll(patchTypes, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (patch.coupled() || fvPatch::constraintType(patch.type()))
        {
            patchTypes[patchI] = referenceField.boundaryField()[patchI].type();
        }
    }

    return patchTypes;
}


wordList diagnosticTensorPatchTypes
(
    const fvMesh& mesh,
    const volTensorField& referenceField
)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        calculatedFvPatchTensorField::typeName
    );

    forAll(patchTypes, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (patch.coupled() || fvPatch::constraintType(patch.type()))
        {
            patchTypes[patchI] = referenceField.boundaryField()[patchI].type();
        }
    }

    return patchTypes;
}


bool hasLandProblem3DiagnosticsDict(const dictionary& dict)
{
    return dict.found("landProblem3Diagnostics");
}


Switch landProblem3DiagnosticsSwitch
(
    const dictionary& dict,
    const word& keyword,
    const bool defaultValue
)
{
    if (!hasLandProblem3DiagnosticsDict(dict))
    {
        return Switch(defaultValue);
    }

    return
        dict.subDict("landProblem3Diagnostics").lookupOrDefault<Switch>
        (
            keyword,
            Switch(defaultValue)
        );
}


scalar landProblem3DiagnosticsScalar
(
    const dictionary& dict,
    const word& keyword,
    const scalar defaultValue
)
{
    if (!hasLandProblem3DiagnosticsDict(dict))
    {
        return defaultValue;
    }

    return
        dict.subDict("landProblem3Diagnostics").lookupOrDefault<scalar>
        (
            keyword,
            defaultValue
        );
}


vector affineDisplacement
(
    const tensor& displacementMap,
    const vector& translation,
    const vector& X
)
{
    return (displacementMap & X) + translation;
}


tensor gradDFromDisplacementMap(const tensor& displacementMap)
{
    return displacementMap.T();
}


tensor skewDisplacementMap(const vector& omega)
{
    return tensor
    (
        0.0,       -omega.z(),  omega.y(),
        omega.z(),  0.0,       -omega.x(),
       -omega.y(),  omega.x(),  0.0
    );
}


tensor zRotation(const scalar angle)
{
    const scalar c = cos(angle);
    const scalar s = sin(angle);

    return tensor
    (
        c,   -s,  0.0,
        s,    c,  0.0,
        0.0,  0.0, 1.0
    );
}


vector quadraticDisplacement(const vector& X)
{
    const scalar x = X.x();
    const scalar y = X.y();
    const scalar z = X.z();

    return vector
    (
        20.0*sqr(x) - 12.0*y*z,
       -15.0*sqr(y) + 10.0*x*z,
        18.0*sqr(z) - 8.0*x*y
    );
}


tensor quadraticGradD(const vector& X)
{
    const scalar x = X.x();
    const scalar y = X.y();
    const scalar z = X.z();

    // OpenFOAM stores grad(D) with rows as spatial derivative directions and
    // columns as displacement components; for D = B & X, grad(D) = B.T().
    return tensor
    (
        40.0*x,  10.0*z,  -8.0*y,
       -12.0*z, -30.0*y,  -8.0*x,
       -12.0*y,  10.0*x,  36.0*z
    );
}

}


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


tmp<volScalarField>
nonLinGeomTotalLagTotalDispSolid::mixedVolumetricConstraint
(
    const volScalarField& J
) const
{
    // Existing generic total-Lagrangian mixed constraint:
    // g(J) = 0.5*(J - 1/J).
    return 0.5*(pow(J, 2.0) - 1.0)/J;
}


scalar
nonLinGeomTotalLagTotalDispSolid::mixedVolumetricConstraintValue
(
    const scalar J
) const
{
    return 0.5*(J - 1.0/J);
}


scalar
nonLinGeomTotalLagTotalDispSolid::mixedVolumetricConstraintDerivative
(
    const scalar J
) const
{
    return 0.5*(1.0 + 1.0/sqr(J));
}


string
nonLinGeomTotalLagTotalDispSolid::mixedVolumetricConstraintDescription() const
{
    return "0.5*(J^2 - 1)/J";
}


void nonLinGeomTotalLagTotalDispSolid::applyMixedPressureStressSplit()
{
    static bool printedActiveStressPreservingSplit = false;

    if (retainFullPassiveStressInMixedSplit())
    {
        // Isolated opt-in path: the material stress is already the complete
        // non-volumetric Cauchy stress and pressure is inserted exactly once.
        sigma() = sigma() - p()*I;
        return;
    }

    const PtrList<mechanicalLaw>& laws = mechanical();

    if (laws.size() == 1)
    {
        if (isA<electroMechanicalLaw>(laws[0]))
        {
            const electroMechanicalLaw& activeLaw =
                refCast<const electroMechanicalLaw>(laws[0]);

            if (activeLaw.hasActiveStress())
            {
                if (!printedActiveStressPreservingSplit)
                {
                    Info<< "Using active-stress-preserving mixed pressure split"
                        << endl;
                    printedActiveStressPreservingSplit = true;
                }

                const tmp<volSymmTensorField> tSigmaActive =
                    activeLaw.activeCauchyStress();
                const volSymmTensorField& sigmaActive = tSigmaActive();

                sigma() =
                    dev(sigma() - sigmaActive) + sigmaActive - p()*I;
                return;
            }
        }

        sigma() = dev(sigma());
        sigma() = sigma() - p()*I;
        return;
    }

    bool hasActiveStress = false;
    forAll(laws, lawI)
    {
        if (isA<electroMechanicalLaw>(laws[lawI]))
        {
            const electroMechanicalLaw& activeLaw =
                refCast<const electroMechanicalLaw>(laws[lawI]);

            if (activeLaw.hasActiveStress())
            {
                hasActiveStress = true;
                break;
            }
        }
    }

    if (!hasActiveStress)
    {
        sigma() = dev(sigma());
        sigma() = sigma() - p()*I;
        return;
    }

    if (!printedActiveStressPreservingSplit)
    {
        Info<< "Using active-stress-preserving mixed pressure split" << endl;
        printedActiveStressPreservingSplit = true;
    }

    PtrList<volSymmTensorField> subMeshSigmaActive(laws.size());

    forAll(laws, lawI)
    {
        const fvMesh& lawMesh =
            mechanical().solSubMeshes().subMeshes()[lawI].subMesh();

        if (isA<electroMechanicalLaw>(laws[lawI]))
        {
            const electroMechanicalLaw& activeLaw =
                refCast<const electroMechanicalLaw>(laws[lawI]);

            if (activeLaw.hasActiveStress())
            {
                const tmp<volSymmTensorField> tSigmaActive =
                    activeLaw.activeCauchyStress();

                subMeshSigmaActive.set
                (
                    lawI,
                    new volSymmTensorField(tSigmaActive())
                );
                continue;
            }
        }

        subMeshSigmaActive.set
        (
            lawI,
            new volSymmTensorField
            (
                IOobject
                (
                    "sigmaActive",
                    lawMesh.time().timeName(),
                    lawMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawMesh,
                dimensionedSymmTensor
                (
                    "zero",
                    dimPressure,
                    symmTensor::zero
                )
            )
        );
    }

    volSymmTensorField sigmaActive
    (
        IOobject
        (
            "sigmaActive",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    mechanical().solSubMeshes().mapSubMeshVolFields<symmTensor>
    (
        subMeshSigmaActive,
        sigmaActive
    );

    sigma() = dev(sigma() - sigmaActive) + sigmaActive - p()*I;
}


void nonLinGeomTotalLagTotalDispSolid::completeDirectFaceCauchyStress
(
    surfaceSymmTensorField&
) const
{}


void nonLinGeomTotalLagTotalDispSolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    // D() = D().oldTime() + U()*runTime().deltaT();
    D() = D().oldTime() + U()*runTime().deltaT()
        + 0.5*sqr(runTime().deltaT())*A_;

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);
    updateJacobianDiagnostics();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        applyMixedPressureStressSplit();
    }
}


void nonLinGeomTotalLagTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& force,
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        vectorField& forceP = boundaryFieldRef(force)[patchI];

        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // traction.boundaryFieldRef()[patchI] =
            //     tracPatch.traction() - nPatch*tracPatch.pressure();
            if (tracPatch.useUndeformedArea())
            {
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfCurrentPatch;
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<slipFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Set shear traction to zero
            // traction.boundaryFieldRef()[patchI] =
                // sqr(nPatch) & traction.boundaryField()[patchI];
            forceP = sqr(nPatch) & force.boundaryField()[patchI];
        }
    }
}


Foam::tmp<Foam::volVectorField>
Foam::solidModels::nonLinGeomTotalLagTotalDispSolid::
momentumSurfaceForceDensity
(
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    const fvMesh& mesh = this->mesh();

    surfaceSymmTensorField sigmaFace
    (
        IOobject
        (
            "faceSigmaDirectConstitutive",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );

    if (faceStressTreatment_ == "interpolatedCell")
    {
        sigmaFace = fvc::interpolate(sigma());
    }
    else if (faceStressTreatment_ == "directConstitutive")
    {
        surfaceTensorField gradDf
        (
            IOobject
            (
                "grad(D)f",
                mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimless, tensor::zero),
            calculatedFvPatchTensorField::typeName
        );
        const_cast<mechanicalModel&>(mechanical()).grad(D, pointD(), gradDf);
        const_cast<mechanicalModel&>(mechanical()).correct(sigmaFace);
        completeDirectFaceCauchyStress(sigmaFace);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown faceStressTreatment " << faceStressTreatment_
            << ". Valid options are interpolatedCell and directConstitutive"
            << exit(FatalError);
    }

    const surfaceVectorField tractionPhysical(nCurrent & sigmaFace);

    surfaceVectorField forcePhysical(magSfCurrent*tractionPhysical);
    enforceTractionBoundaries(forcePhysical, D, nCurrent, magSfCurrent);

    momentumStabilisation().updateVector(D, &gradD());

    const surfaceVectorField tractionWithStab
    (
        tractionPhysical + impKf_*momentumStabilisation().faceVector()
    );

    surfaceVectorField forceWithStab(magSfCurrent*tractionWithStab);
    enforceTractionBoundaries(forceWithStab, D, nCurrent, magSfCurrent);

    // Isolate the stabilisation that remains after solidTraction patches have
    // overwritten face forces and symmetry/slip patches have projected shear.
    // The optional boundary treatment is then applied to this isolated
    // numerical contribution only; physical stress forces are unchanged.
    surfaceVectorField forceStab(forceWithStab - forcePhysical);
    applyMomentumStabilisationBoundaryTreatment(forceStab);

    vectorField stabCellForce(fvc::div(forceStab));
    stabCellForce *= mesh.V();

    momentumStabilisation().projectExtensiveVectorForce(stabCellForce);

    tmp<volVectorField> tMomentumSurfaceForceDensity
    (
        fvc::div(forcePhysical)
    );

    volVectorField& momentumSurfaceForceDensity =
        tmpRef(tMomentumSurfaceForceDensity);

    if (reportMomentumStabDiagnostics_ || writeMomentumStabFields_)
    {
        reportMomentumStabilisationDiagnostics
        (
            forceStab,
            stabCellForce,
            momentumSurfaceForceDensity
        );
    }

    primitiveFieldRef(momentumSurfaceForceDensity) +=
        stabCellForce/mesh.V();

    return tMomentumSurfaceForceDensity;
}


void nonLinGeomTotalLagTotalDispSolid::momentumResidualDecomposition
(
    momentumResidualDecompositionData& decomposition
) const
{
    // This routine is deliberately an audit path.  It mirrors the production
    // face-force route, but does not write fields, assemble a matrix, or alter
    // any accepted state.  The returned face arrays use polyMesh global face
    // numbering so an external audit can inspect owner/neighbour and patch
    // contributions without reconstructing OpenFOAM surface-field storage.
    const fvMesh& mesh = this->mesh();

    decomposition.names = wordList
    (
        makeList<word>
        ({
            "passive",
            "viscous",
            "pressure",
            "faceInterpolation",
            "epicardialSpring",
            "epicardialDashpot",
            "basalSpring",
            "basalDashpot",
            "otherTraction",
            "momentumStabilisation",
            "inertia",
            "damping"
        })
    );

    const label nTerms = decomposition.names.size();
    decomposition.cellActions.setSize(nTerms);
    decomposition.faceForces.setSize(nTerms);
    for (label termI = 0; termI < nTerms; ++termI)
    {
        decomposition.cellActions[termI].setSize(mesh.nCells(), vector::zero);
        decomposition.faceForces[termI].setSize(mesh.nFaces(), vector::zero);
    }

    decomposition.passiveFaceStress.setSize(mesh.nFaces(), symmTensor::zero);
    decomposition.viscousFaceStress.setSize(mesh.nFaces(), symmTensor::zero);
    decomposition.pressureFaceStress.setSize(mesh.nFaces(), symmTensor::zero);
    decomposition.directPassiveFaceStress.setSize
    (
        mesh.nFaces(),
        symmTensor::zero
    );
    decomposition.directViscousFaceStress.setSize
    (
        mesh.nFaces(),
        symmTensor::zero
    );
    decomposition.faceAreaVectors.setSize(mesh.nFaces(), vector::zero);

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    surfaceSymmTensorField directMaterialStress
    (
        IOobject
        (
            "arosticaAuditDirectMaterialStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    surfaceTensorField directGradDf
    (
        IOobject
        (
            "grad(D)f",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero),
        calculatedFvPatchTensorField::typeName
    );
    const_cast<mechanicalModel&>(mechanical()).grad
    (
        D(),
        pointD(),
        directGradDf
    );
    const_cast<mechanicalModel&>(mechanical()).correct(directMaterialStress);

    surfaceSymmTensorField directViscousStress
    (
        IOobject
        (
            "arosticaAuditDirectViscousStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    volSymmTensorField cellViscousStress
    (
        IOobject
        (
            "arosticaAuditCellViscousStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );

    // The Aróstica viscoelastic law registers these diagnostic fields even
    // when writing is disabled.  Looking them up by their documented name
    // keeps this generic audit method free of an Aróstica-only dynamic cast.
    const PtrList<mechanicalLaw>& laws = mechanical();
    forAll(laws, lawI)
    {
        const word cellFieldName
        (
            word("ArosticasigmaViscous_") + laws[lawI].name()
        );
        const word faceFieldName
        (
            word("ArosticasigmaViscousf_") + laws[lawI].name()
        );

        if (mesh.foundObject<volSymmTensorField>(cellFieldName))
        {
            cellViscousStress +=
                mesh.lookupObject<volSymmTensorField>(cellFieldName);
        }

        if (mesh.foundObject<surfaceSymmTensorField>(faceFieldName))
        {
            directViscousStress +=
                mesh.lookupObject<surfaceSymmTensorField>(faceFieldName);
        }
    }

    surfaceSymmTensorField selectedMaterialStress
    (
        IOobject
        (
            "arosticaAuditSelectedMaterialStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    surfaceSymmTensorField selectedViscousStress
    (
        IOobject
        (
            "arosticaAuditSelectedViscousStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    surfaceSymmTensorField selectedPressureStress
    (
        IOobject
        (
            "arosticaAuditSelectedPressureStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    surfaceSymmTensorField directPressureStress
    (
        IOobject
        (
            "arosticaAuditDirectPressureStress",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );

    if (faceStressTreatment_ == "interpolatedCell")
    {
        selectedMaterialStress = fvc::interpolate(sigma());
        selectedViscousStress = fvc::interpolate(cellViscousStress);

        if (solvePressure())
        {
            selectedPressureStress =
                -fvc::interpolate(p())
               *dimensionedSymmTensor("I", dimless, I);
        }
    }
    else if (faceStressTreatment_ == "directConstitutive")
    {
        const surfaceSymmTensorField directStressBeforeCompletion
        (
            directMaterialStress
        );
        completeDirectFaceCauchyStress(directMaterialStress);
        directPressureStress =
            directMaterialStress - directStressBeforeCompletion;
        selectedMaterialStress = directMaterialStress;
        selectedViscousStress = directViscousStress;
        selectedPressureStress = directPressureStress;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown faceStressTreatment " << faceStressTreatment_
            << exit(FatalError);
    }

    const surfaceSymmTensorField directPassiveStress
    (
        directMaterialStress - directViscousStress - directPressureStress
    );
    const surfaceSymmTensorField selectedPassiveStress
    (
        selectedMaterialStress
      - selectedViscousStress
      - selectedPressureStress
    );

    surfaceVectorField passiveForce
    (
        IOobject
        (
            "arosticaAuditPassiveForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        SfCurrent & directPassiveStress
    );
    surfaceVectorField viscousForce
    (
        IOobject
        (
            "arosticaAuditViscousForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        SfCurrent & directViscousStress
    );
    surfaceVectorField pressureForce
    (
        IOobject
        (
            "arosticaAuditPressureForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        SfCurrent & selectedPressureStress
    );
    surfaceVectorField interpolationForce
    (
        IOobject
        (
            "arosticaAuditFaceInterpolationForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        SfCurrent & (selectedPassiveStress + selectedViscousStress
                   - directPassiveStress - directViscousStress)
    );
    surfaceVectorField epicardialSpringForce
    (
        IOobject
        (
            "arosticaAuditEpicardialSpringForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );
    surfaceVectorField epicardialDashpotForce
    (
        IOobject
        (
            "arosticaAuditEpicardialDashpotForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );
    surfaceVectorField basalSpringForce
    (
        IOobject
        (
            "arosticaAuditBasalSpringForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );
    surfaceVectorField basalDashpotForce
    (
        IOobject
        (
            "arosticaAuditBasalDashpotForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );
    surfaceVectorField otherTractionForce
    (
        IOobject
        (
            "arosticaAuditOtherTractionForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );

    // Replace physical forces on traction boundaries with the exact
    // stateless boundary contributions, and apply the same projection used
    // by the production residual on fixed/shear-free boundaries.
    forAll(mesh.boundary(), patchI)
    {
        const fvPatchVectorField& DPatch = D().boundaryField()[patchI];
        const vectorField& nPatch = nCurrent.boundaryField()[patchI];

        if (isA<solidTractionFvPatchVectorField>(DPatch))
        {
            passiveForce.boundaryFieldRef()[patchI] = vector::zero;
            viscousForce.boundaryFieldRef()[patchI] = vector::zero;
            pressureForce.boundaryFieldRef()[patchI] = vector::zero;
            interpolationForce.boundaryFieldRef()[patchI] = vector::zero;

            const solidTractionFvPatchVectorField& tractionPatch =
                refCast<const solidTractionFvPatchVectorField>(DPatch);
            scalarField patchArea
            (
                tractionPatch.useUndeformedArea()
              ? mesh.boundary()[patchI].magSf()
              : magSfCurrent.boundaryField()[patchI]
            );

            if
            (
                isA<arosticaSpringDashpotTractionFvPatchVectorField>(DPatch)
            )
            {
                const arosticaSpringDashpotTractionFvPatchVectorField& bc =
                    refCast
                    <const arosticaSpringDashpotTractionFvPatchVectorField>
                    (DPatch);
                const tmp<volVectorField> tDdot(fvc::ddt(D()));
                vectorField spring;
                vectorField dashpot;
                bc.diagnosticContributions
                (
                    spring,
                    dashpot,
                    DPatch,
                    tDdot().boundaryField()[patchI]
                );

                if (isA<arosticaNormalSpringDashpotTractionFvPatchVectorField>
                    (DPatch))
                {
                    epicardialSpringForce.boundaryFieldRef()[patchI] =
                        spring*patchArea;
                    epicardialDashpotForce.boundaryFieldRef()[patchI] =
                        dashpot*patchArea;
                }
                else
                {
                    basalSpringForce.boundaryFieldRef()[patchI] =
                        spring*patchArea;
                    basalDashpotForce.boundaryFieldRef()[patchI] =
                        dashpot*patchArea;
                }
            }
            else
            {
                otherTractionForce.boundaryFieldRef()[patchI] =
                    (
                        tractionPatch.traction()
                      - nPatch*tractionPatch.pressure()
                    )*patchArea;
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>(DPatch)
         || isA<symmetryFvPatchVectorField>(DPatch)
         || isA<slipFvPatchVectorField>(DPatch)
        )
        {
            passiveForce.boundaryFieldRef()[patchI] =
                sqr(nPatch) & passiveForce.boundaryField()[patchI];
            viscousForce.boundaryFieldRef()[patchI] =
                sqr(nPatch) & viscousForce.boundaryField()[patchI];
            pressureForce.boundaryFieldRef()[patchI] =
                sqr(nPatch) & pressureForce.boundaryField()[patchI];
            interpolationForce.boundaryFieldRef()[patchI] =
                sqr(nPatch) & interpolationForce.boundaryField()[patchI];
        }
    }

    surfaceVectorField physicalForce
    (
        IOobject
        (
            "arosticaAuditPhysicalForce",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        passiveForce
      + viscousForce
      + pressureForce
      + interpolationForce
      + epicardialSpringForce
      + epicardialDashpotForce
      + basalSpringForce
      + basalDashpotForce
      + otherTractionForce
    );

    momentumStabilisation().updateVector(D(), &gradD());

    const surfaceVectorField forceWithStabRaw
    (
        physicalForce
      + magSfCurrent*impKf_*momentumStabilisation().faceVector()
    );
    surfaceVectorField forceWithStab(forceWithStabRaw);
    enforceTractionBoundaries(forceWithStab, D(), nCurrent, magSfCurrent);
    surfaceVectorField stabilisationForce(forceWithStab - physicalForce);
    applyMomentumStabilisationBoundaryTreatment(stabilisationForce);

    vectorField stabilisationCellForce(fvc::div(stabilisationForce));
    stabilisationCellForce *= mesh.V();
    momentumStabilisation().projectExtensiveVectorForce
    (
        stabilisationCellForce
    );

    const surfaceVectorField* faceFields[] =
    {
        &passiveForce,
        &viscousForce,
        &pressureForce,
        &interpolationForce,
        &epicardialSpringForce,
        &epicardialDashpotForce,
        &basalSpringForce,
        &basalDashpotForce,
        &otherTractionForce,
        &stabilisationForce
    };

    for (label termI = 0; termI < 10; ++termI)
    {
        const surfaceVectorField& faceField = *faceFields[termI];
        const tmp<volVectorField> tCellDiv(fvc::div(faceField));
        decomposition.cellActions[termI] = tCellDiv().internalField();

        forAll(faceField, faceI)
        {
            decomposition.faceForces[termI][faceI] = faceField[faceI];
        }
        forAll(faceField.boundaryField(), patchI)
        {
            const vectorField& patchField =
                faceField.boundaryField()[patchI];
            forAll(patchField, faceI)
            {
                decomposition.faceForces[termI]
                    [mesh.boundary()[patchI].start() + faceI] =
                    patchField[faceI];
            }
        }
    }

    const tmp<volVectorField> tInertia
    (
        -rho()*fvc::d2dt2(D())
    );
    const tmp<volVectorField> tDamping
    (
        -rho()*dampingCoeff()*fvc::ddt(D())
    );
    decomposition.cellActions[10] = tInertia().internalField();
    decomposition.cellActions[11] = tDamping().internalField();

    forAll(SfCurrent, faceI)
    {
        decomposition.faceAreaVectors[faceI] = SfCurrent[faceI];
    }
    forAll(SfCurrent.boundaryField(), patchI)
    {
        const vectorField& patchField = SfCurrent.boundaryField()[patchI];
        forAll(patchField, faceI)
        {
            decomposition.faceAreaVectors
                [mesh.boundary()[patchI].start() + faceI] = patchField[faceI];
        }
    }

    auto flattenStress =
    [&mesh](const surfaceSymmTensorField& field, symmTensorField& flat)
    {
        forAll(field, faceI)
        {
            flat[faceI] = field[faceI];
        }
        forAll(field.boundaryField(), patchI)
        {
            const symmTensorField& patchField =
                field.boundaryField()[patchI];
            forAll(patchField, faceI)
            {
                flat[mesh.boundary()[patchI].start() + faceI] =
                    patchField[faceI];
            }
        }
    };

    flattenStress(selectedPassiveStress, decomposition.passiveFaceStress);
    flattenStress(selectedViscousStress, decomposition.viscousFaceStress);
    flattenStress(selectedPressureStress, decomposition.pressureFaceStress);
    flattenStress(directPassiveStress, decomposition.directPassiveFaceStress);
    flattenStress(directViscousStress, decomposition.directViscousFaceStress);
}


void nonLinGeomTotalLagTotalDispSolid::reportMomentumStabilisationDiagnostics
(
    const surfaceVectorField& momentumStabFaceForce,
    const vectorField& momentumStabCellForce,
    const volVectorField& physicalSurfaceForceDensity
) const
{
    momentumStabDiagnosticCall_++;

    if
    (
        momentumStabDiagnosticCall_ != 1
     && (momentumStabDiagnosticCall_ % momentumStabReportInterval_) != 0
    )
    {
        return;
    }

    const fvMesh& mesh = this->mesh();

    const vectorField momentumStabForceDensity
    (
        momentumStabCellForce/mesh.V()
    );

    const vectorField physicalCellForce
    (
        primitiveField(physicalSurfaceForceDensity)*mesh.V()
    );

    label nCells = mesh.nCells();
    reduce(nCells, sumOp<label>());

    const vector netStabForce = gSum(momentumStabCellForce);
    const scalar sumMagStabForce = gSum(mag(momentumStabCellForce));
    const scalar maxStabForce = gMax(mag(momentumStabCellForce));
    const scalar rmsStabForceDensity =
        Foam::sqrt(gSum(magSqr(momentumStabForceDensity))/scalar(nCells));

    const scalar sumMagPhysicalForce = gSum(mag(physicalCellForce));
    const scalar rmsPhysicalForceDensity =
        Foam::sqrt
        (
            gSum(magSqr(primitiveField(physicalSurfaceForceDensity)))
           /scalar(nCells)
        );

    if (reportMomentumStabDiagnostics_)
    {
        Info<< "Momentum stabilisation diagnostics:"
            << " call=" << momentumStabDiagnosticCall_ << nl
            << "    net extensive stabilisation force = "
            << netStabForce << nl
            << "    sum of extensive cell-force magnitudes = "
            << sumMagStabForce << nl
            << "    maximum extensive cell force = "
            << maxStabForce << nl
            << "    RMS stabilisation force density = "
            << rmsStabForceDensity << nl
            << "    sumMag(stabilisation force)/sumMag(physical force) = "
            << sumMagStabForce/(sumMagPhysicalForce + VSMALL) << nl
            << "    RMS(stabilisation force density)"
            << "/RMS(physical force density) = "
            << rmsStabForceDensity/(rmsPhysicalForceDensity + VSMALL)
            << endl;
    }

    if (writeMomentumStabFields_)
    {
        surfaceVectorField momentumStabFaceForceWrite(momentumStabFaceForce);
        momentumStabFaceForceWrite.rename("momentumStabFaceForce");
        momentumStabFaceForceWrite.write();

        volVectorField momentumStabForceDensityWrite
        (
            IOobject
            (
                "momentumStabForceDensity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, vector::zero)
        );
        primitiveFieldRef(momentumStabForceDensityWrite) =
            momentumStabForceDensity;
        momentumStabForceDensityWrite.correctBoundaryConditions();
        momentumStabForceDensityWrite.write();

        volVectorField physicalSurfaceForceDensityWrite
        (
            physicalSurfaceForceDensity
        );
        physicalSurfaceForceDensityWrite.rename("physicalSurfaceForceDensity");
        physicalSurfaceForceDensityWrite.write();

        volScalarField momentumStabToPhysicalRatio
        (
            IOobject
            (
                "momentumStabToPhysicalRatio",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        scalarField& ratioI = primitiveFieldRef(momentumStabToPhysicalRatio);
        const vectorField& physicalDensityI =
            primitiveField(physicalSurfaceForceDensity);

        forAll(ratioI, cellI)
        {
            ratioI[cellI] =
                mag(momentumStabForceDensity[cellI])
               /(mag(physicalDensityI[cellI]) + VSMALL);
        }

        momentumStabToPhysicalRatio.correctBoundaryConditions();
        momentumStabToPhysicalRatio.write();
    }
}


bool nonLinGeomTotalLagTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    if (predictor_ && newTimeStep())
    {
        predict();
    }

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    int iCorr = 0;
    scalar currentResidualNorm = 0;
    scalar initialResidualNorm = 0;
    scalar deltaXNorm = 0;
    scalar xNorm = 0;
    const convergenceParameters convParam =
        readConvergenceParameters(solidModelDict());

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Calculate deformed area vectors and normals
        const surfaceVectorField SfCurrent
        (
            fvc::interpolate(J_*Finv_.T()) & mesh().Sf()
        );
        const surfaceScalarField magSfCurrent(mag(SfCurrent));
        const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

        const tmp<volVectorField> tMomentumSurfaceForceDensity
        (
            momentumSurfaceForceDensity(D(), nCurrent, magSfCurrent)
        );

        // Momentum equation total displacement total Lagrangian form
#ifndef OPENFOAM_COM
        // Assemble the RHS in stages.
        // The equivalent chained tmp fvMatrix expression is stable on OpenFOAM.com.
        tmp<fvVectorMatrix> tRhsEqn
        (
            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
        );
        tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
        tmpRef(tRhsEqn) += tMomentumSurfaceForceDensity();
        tmpRef(tRhsEqn) += rho()*g();

        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == tRhsEqn
        );
#else
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + tMomentumSurfaceForceDensity()
          + rho()*g()
        );
#endif

        // Add damping
        if (dampingCoeff().value() > SMALL)
        {
            DEqn += dampingCoeff()*rho()*fvm::ddt(D());
        }

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system and store the residual
        currentResidualNorm = mag(DEqn.solve().initialResidual());

        // Norm of the solution correction
        deltaXNorm =
            sqrt
            (
                gSum
                (
                    magSqr
                    (
#ifdef OPENFOAM_NOT_EXTEND
                        D().primitiveField() - D().prevIter().primitiveField()
#else
                        D().internalField() - D().prevIter().internalField()
#endif
                    )
                )
            );

        // Norm of the solution
#ifdef OPENFOAM_NOT_EXTEND
        xNorm = sqrt(gSum(magSqr(D().primitiveField())));
#else
        xNorm = sqrt(gSum(magSqr(D().internalField())));
#endif

        // Store the initial residual
        if (iCorr == 0)
        {
            initialResidualNorm = currentResidualNorm;
            Info<< "Initial Residual Norm = " << initialResidualNorm << nl
                << "Initial Solution Norm = " << xNorm << endl;
        }

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);
        updateJacobianDiagnostics();

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
        !checkConvergence
        (
            currentResidualNorm,
            initialResidualNorm,
            deltaXNorm,
            xNorm,
            ++iCorr,
            convParam
        )
    );

    updateFinalStateDiagnostics();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool nonLinGeomTotalLagTotalDispSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Seed the PETSc solution vector from the predicted fields
        packSolution(foamPetscSnesHelper::solution());
    }

    if
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "exactAffineResidualTest", false
        )
    )
    {
        evaluateExactAffineResidual();
        return true;
    }

    if
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "exactMmsResidualTest", false
        )
    )
    {
        evaluateExactMmsResidual();
        return true;
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Map the PETSc solution back into the D field (and p when active),
    // refreshing dependent kinematic fields and boundary conditions
    unpackSolution(foamPetscSnesHelper::solution(), false);

    if (solvePressure())
    {
        // Update dpdt
        autoPtrRef(dpdtPtr_) = fvc::ddt(p());
    }

    updateFinalStateDiagnostics();

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


void nonLinGeomTotalLagTotalDispSolid::evaluateExactMmsResidual()
{
    const fvMesh& mesh = this->mesh();
    Info<< "Evaluating the production residual once at DExact/pExact MMS state" << endl;

    volVectorField exactD
    (
        IOobject("DExact", mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField exactP
    (
        IOobject("pExact", mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volVectorField& D = const_cast<volVectorField&>(this->D());
    D = exactD;
    D.correctBoundaryConditions();
    if (solvePressure())
    {
        volScalarField& p = const_cast<volScalarField&>(this->p());
        p = exactP;
        p.correctBoundaryConditions();
    }
    mechanical().grad(D, gradD());
    mechanical().interpolate(D, gradD(), pointD());
    pointD().correctBoundaryConditions();
    F_ = I + gradD().T();
    Finv_ = inv(F_);
    J_ = det(F_);
    mechanical().correct(sigma());
    packSolution(foamPetscSnesHelper::solution());

    Vec residual = nullptr;
    AssertPETSc(VecDuplicate(foamPetscSnesHelper::solution(), &residual));
    AssertPETSc(formResidual(residual, foamPetscSnesHelper::solution()));

    PetscInt blockSize = 0;
    PetscInt localSize = 0;
    AssertPETSc(VecGetBlockSize(residual, &blockSize));
    AssertPETSc(VecGetLocalSize(residual, &localSize));
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetArrayRead(residual, &values));
    vectorField momentum(mesh.nCells(), vector::zero);
    scalarField pressure(mesh.nCells(), 0.0);
    scalar momSqr = 0.0;
    scalar pSqr = 0.0;
    scalar momInf = 0.0;
    scalar pInf = 0.0;
    for (label cellI = 0; cellI < mesh.nCells(); ++cellI)
    {
        const PetscInt offset = cellI*blockSize;
        for (label cmpt = 0; cmpt < 3; ++cmpt)
        {
            momentum[cellI][cmpt] = PetscRealPart(values[offset + cmpt])/mesh.V()[cellI];
        }
        if (blockSize > 3)
        {
            pressure[cellI] = PetscRealPart(values[offset + 3])/mesh.V()[cellI];
        }
        momSqr += mesh.V()[cellI]*magSqr(momentum[cellI]);
        pSqr += mesh.V()[cellI]*sqr(pressure[cellI]);
        momInf = max(momInf, mag(momentum[cellI]));
        pInf = max(pInf, mag(pressure[cellI]));
    }
    AssertPETSc(VecRestoreArrayRead(residual, &values));
    AssertPETSc(VecDestroy(&residual));
    const scalar volume = gSum(mesh.V());
    volVectorField exactMomentumResidual
    (
        IOobject("exactMmsMomentumResidual", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );
    exactMomentumResidual.primitiveFieldRef() = momentum;
    exactMomentumResidual.write();
    volScalarField exactPressureResidual
    (
        IOobject("exactMmsPressureResidual", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("zero", dimless, 0)
    );
    exactPressureResidual.primitiveFieldRef() = pressure;
    exactPressureResidual.write();
    Info<< "S4F_EXACT_MMS_RESIDUAL,momentumRMS="
        << sqrt(momSqr/(volume + VSMALL))
        << ",momentumLinf=" << momInf
        << ",pressureRMS=" << sqrt(pSqr/(volume + VSMALL))
        << ",pressureLinf=" << pInf << endl;
}


void nonLinGeomTotalLagTotalDispSolid::evaluateExactAffineResidual()
{
    const fvMesh& mesh = this->mesh();
    const scalar gamma =
        solidModelDict().lookupOrDefault<scalar>("exactAffineGamma", 0.2);

    Info<< "Evaluating the production residual once at exact affine state"
        << " (gamma=" << gamma << ")" << endl;

    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
    forAll(DI, cellI)
    {
        DI[cellI] = vector(gamma*mesh.C()[cellI].z(), 0, 0);
    }
    D.correctBoundaryConditions();

    if (solvePressure())
    {
        volScalarField& p = const_cast<volScalarField&>(this->p());
        p = dimensionedScalar("zero", p.dimensions(), 0);
        p.correctBoundaryConditions();
    }

    // Refresh exactly the same dependent fields used by the production
    // residual before packing the state into PETSc.
    mechanical().grad(D, gradD());
    mechanical().interpolate(D, gradD(), pointD());
    pointD().correctBoundaryConditions();
    F_ = I + gradD().T();
    Finv_ = inv(F_);
    J_ = det(F_);
    updateJacobianDiagnostics();
    mechanical().correct(sigma());
    packSolution(foamPetscSnesHelper::solution());

    Vec residual = nullptr;
    AssertPETSc
    (
        VecDuplicate(foamPetscSnesHelper::solution(), &residual)
    );
    AssertPETSc
    (
        formResidual(residual, foamPetscSnesHelper::solution())
    );

    // These are the same decomposition fields used by the completed-state
    // diagnostics, refreshed after the production residual call.
    updateQuasiStaticDiagnostics();
    if (solvePressure())
    {
        updatePressureConstraintDiagnostics();
    }

    // Preserve the two legitimate face-stress constructions at this exact
    // state.  This is a diagnostic only; the production residual above has
    // already been evaluated and is not replaced by these fields.
    surfaceSymmTensorField exactFaceStress
    (
        IOobject
        (
            "exactAffineFaceStress",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        calculatedFvPatchSymmTensorField::typeName
    );
    if (faceStressTreatment_ == "interpolatedCell")
    {
        exactFaceStress = fvc::interpolate(sigma());
    }
    else
    {
        surfaceTensorField gradDf
        (
            IOobject
            (
                // The mechanical-law surface correct() overload obtains the
                // current surface gradient by this production registry name.
                // Keep the diagnostic field name identical so direct
                // constitutive evaluation exercises the same face-kernel
                // path as the selectable production treatment.
                "grad(D)f",
                mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimless, tensor::zero),
            calculatedFvPatchTensorField::typeName
        );
        mechanical().grad(D, pointD(), gradDf);
        mechanical().correct(exactFaceStress);
    }
    surfaceSymmTensorField interpolatedFaceStress
    (
        IOobject
        (
            "exactAffineInterpolatedFaceStress",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        fvc::interpolate(sigma())
    );
    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);
    surfaceVectorField exactFaceTraction
    (
        IOobject
        (
            "exactAffineFaceTraction",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        nCurrent & exactFaceStress
    );
    surfaceVectorField interpolatedFaceTraction
    (
        IOobject
        (
            "exactAffineInterpolatedFaceTraction",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        nCurrent & interpolatedFaceStress
    );
    scalar maxFaceStressDifference = 0;
    scalar maxFaceTractionDifference = 0;
    forAll(exactFaceStress, faceI)
    {
        maxFaceStressDifference = max
        (
            maxFaceStressDifference,
            mag(exactFaceStress[faceI] - interpolatedFaceStress[faceI])
        );
        maxFaceTractionDifference = max
        (
            maxFaceTractionDifference,
            mag(exactFaceTraction[faceI] - interpolatedFaceTraction[faceI])
        );
    }
    forAll(exactFaceStress.boundaryField(), patchI)
    {
        forAll(exactFaceStress.boundaryField()[patchI], faceI)
        {
            maxFaceStressDifference = max
            (
                maxFaceStressDifference,
                mag
                (
                    exactFaceStress.boundaryField()[patchI][faceI]
                  - interpolatedFaceStress.boundaryField()[patchI][faceI]
                )
            );
            maxFaceTractionDifference = max
            (
                maxFaceTractionDifference,
                mag
                (
                    exactFaceTraction.boundaryField()[patchI][faceI]
                  - interpolatedFaceTraction.boundaryField()[patchI][faceI]
                )
            );
        }
    }
    reduce(maxFaceStressDifference, maxOp<scalar>());
    reduce(maxFaceTractionDifference, maxOp<scalar>());
    exactFaceStress.write();
    interpolatedFaceStress.write();
    exactFaceTraction.write();
    interpolatedFaceTraction.write();
    Info<< "S4F_EXACT_AFFINE_FACE,"
        << "treatment=" << faceStressTreatment_ << ','
        << "maxStressDifference=" << maxFaceStressDifference << ','
        << "maxTractionDifference=" << maxFaceTractionDifference
        << endl;

    PetscInt blockSize = 0;
    PetscInt localSize = 0;
    AssertPETSc(VecGetBlockSize(residual, &blockSize));
    AssertPETSc(VecGetLocalSize(residual, &localSize));
    const PetscScalar* values = nullptr;
    AssertPETSc(VecGetArrayRead(residual, &values));

    volVectorField exactMomentumResidual
    (
        IOobject
        (
            "exactAffineMomentumResidual",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );
    volScalarField exactPressureResidual
    (
        IOobject
        (
            "exactAffinePressureResidual",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    vectorField& exactMomentumI = exactMomentumResidual.primitiveFieldRef();
    scalarField& exactPressureI = exactPressureResidual.primitiveFieldRef();
    const scalarField& V = mesh.V();

    scalar rawMomentumSqr = 0;
    scalar rawPressureSqr = 0;
    scalar momentumSqr = 0;
    scalar pressureSqr = 0;
    scalar momentumLinf = 0;
    scalar pressureLinf = 0;
    forAll(exactMomentumI, cellI)
    {
        const PetscInt i = blockSize*cellI;
        vector mom
        (
            PetscRealPart(values[i]),
            PetscRealPart(values[i + 1]),
            PetscRealPart(values[i + 2])
        );
        const scalar pResidual =
            solvePressure() && blockSize > 3
          ? PetscRealPart(values[i + 3])
          : 0;

        exactMomentumI[cellI] = mom/V[cellI];
        exactPressureI[cellI] =
            solvePressure()
          ? pResidual/(V[cellI]*pressureEqnScale_)
          : 0;

        rawMomentumSqr += magSqr(mom);
        rawPressureSqr += sqr(pResidual);
        momentumSqr += V[cellI]*magSqr(exactMomentumI[cellI]);
        pressureSqr += V[cellI]*sqr(exactPressureI[cellI]);
        momentumLinf = max(momentumLinf, mag(exactMomentumI[cellI]));
        pressureLinf = max(pressureLinf, mag(exactPressureI[cellI]));
    }
    AssertPETSc(VecRestoreArrayRead(residual, &values));

    reduce(rawMomentumSqr, sumOp<scalar>());
    reduce(rawPressureSqr, sumOp<scalar>());
    reduce(momentumSqr, sumOp<scalar>());
    reduce(pressureSqr, sumOp<scalar>());
    reduce(momentumLinf, maxOp<scalar>());
    reduce(pressureLinf, maxOp<scalar>());
    exactMomentumResidual.correctBoundaryConditions();
    exactPressureResidual.correctBoundaryConditions();
    exactMomentumResidual.write();
    exactPressureResidual.write();

    volVectorField exactPhysicalForceDensity
    (
        IOobject
        (
            "exactAffinePhysicalStressForceDensity",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        physicalStressForceDensity_
    );
    volVectorField exactStabilisationForceDensity
    (
        IOobject
        (
            "exactAffineMomentumStabilisationForceDensity",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        momentumStabForceDensity_
    );
    exactPhysicalForceDensity.write();
    exactStabilisationForceDensity.write();

    volScalarField exactPressureConstraintResidual
    (
        IOobject
        (
            "exactAffinePressureConstraintResidual",
            mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        pressureConstraintResidual_
    );
    exactPressureConstraintResidual.write();

    scalar volume = gSum(V);
    scalar physicalSqr = 0;
    scalar stabilisationSqr = 0;
    vector netPhysical = vector::zero;
    vector netStabilisation = vector::zero;
    forAll(V, cellI)
    {
        physicalSqr += V[cellI]*magSqr(physicalStressForceDensity_[cellI]);
        stabilisationSqr += V[cellI]*magSqr(momentumStabForceDensity_[cellI]);
        netPhysical += V[cellI]*physicalStressForceDensity_[cellI];
        netStabilisation += V[cellI]*momentumStabForceDensity_[cellI];
    }
    reduce(physicalSqr, sumOp<scalar>());
    reduce(stabilisationSqr, sumOp<scalar>());
    reduce(netPhysical, sumOp<vector>());
    reduce(netStabilisation, sumOp<vector>());

    Info<< "S4F_EXACT_AFFINE_RESIDUAL,"
        << "gamma=" << gamma << ','
        << "cells=" << mesh.nCells() << ','
        << "rawMomentumL2=" << sqrt(rawMomentumSqr) << ','
        << "rawPressureL2=" << sqrt(rawPressureSqr) << ','
        << "momentumRms=" << sqrt(momentumSqr/(volume + VSMALL)) << ','
        << "momentumLinf=" << momentumLinf << ','
        << "pressureRms=" << sqrt(pressureSqr/(volume + VSMALL)) << ','
        << "pressureLinf=" << pressureLinf << ','
        << "physicalForceRms=" << sqrt(physicalSqr/(volume + VSMALL)) << ','
        << "stabilisationForceRms="
        << sqrt(stabilisationSqr/(volume + VSMALL)) << ','
        << "netPhysical=" << netPhysical << ','
        << "netStabilisation=" << netStabilisation
        << endl;

    AssertPETSc(VecDestroy(&residual));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagTotalDispSolid::nonLinGeomTotalLagTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    nonLinGeomTotalLagTotalDispSolid(typeName, runTime, region)
{}


nonLinGeomTotalLagTotalDispSolid::nonLinGeomTotalLagTotalDispSolid
(
    const word& modelType,
    Time& runTime,
    const word& region
)
:
    solidModel(modelType, runTime, region),
    foamPetscSnesHelper
    (
        "D",
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::CELLS,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
    useBodyForceField_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "useBodyForceField", false
        )
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        det(F_)
    ),
    Jminus1_
    (
        IOobject
        (
            "Jminus1",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        J_ - dimensionedScalar("one", dimless, 1.0)
    ),
    Jgeom_
    (
        IOobject
        (
            "Jgeom",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    JgradMinusJgeom_
    (
        IOobject
        (
            "JgradMinusJgeom",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    JgeomMinus1_
    (
        IOobject
        (
            "JgeomMinus1",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureCompressibilityTerm_
    (
        IOobject
        (
            "pressureCompressibilityTerm",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureStabilisationTerm_
    (
        IOobject
        (
            "pressureStabilisationTerm",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    volumetricConstraintTerm_
    (
        IOobject
        (
            "volumetricConstraintTerm",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureConstraintResidual_
    (
        IOobject
        (
            "pressureConstraintResidual",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    inertiaForceDensity_
    (
        IOobject
        (
            "inertiaForceDensity",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    ),
    dampingForceDensity_
    (
        IOobject
        (
            "dampingForceDensity",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    ),
    bodyForceDensity_
    (
        IOobject
        (
            "bodyForceDensity",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    ),
    physicalStressForceDensity_
    (
        IOobject
        (
            "physicalStressForceDensity",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    ),
    momentumStabForceDensity_
    (
        IOobject
        (
            "momentumStabForceDensity",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    ),
    highJgradMask_
    (
        IOobject
        (
            "highJgradMask",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    lowJgradMask_
    (
        IOobject
        (
            "lowJgradMask",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    highJgeomMask_
    (
        IOobject
        (
            "highJgeomMask",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    lowJgeomMask_
    (
        IOobject
        (
            "lowJgeomMask",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    largeJgradMinusJgeomMask_
    (
        IOobject
        (
            "largeJgradMinusJgeomMask",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    A_
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::d2dt2(D())
    ),
    impK_
    (
        solvePressure()
      ? 2.0*mechanical().shearModulus()
      : mechanical().impK()
    ),
    impKf_(fvc::interpolate(impK_)),
    rImpK_(1.0/impK_),
    rKappaPtr_(),
    dpdtPtr_(),
    rAUfTimeIndex_(-1),
    rAUfDeltaT_(0),
    scaleMixedPetScFields_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "scaleMixedPetScFields", true
        )
    ),
    pressureUnknownScaleType_
    (
        solidModelDict().lookupOrDefault<word>("pressureUnknownScale", "twoMu")
    ),
    pressureUnknownScale_(1.0),
    pressureScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("pressureScaleFactor", 1.0)
    ),
    pressureScaleByTwoMu_
    (
        solidModelDict().lookupOrDefault<Switch>("pressureScaleByTwoMu", true)
    ),
    twoMuRef_(1.0),
    pressureEqnScale_(pressureScaleFactor_),
    pressureRowScaling_
    (
        solidModelDict().lookupOrDefault<word>
        (
            "pressureRowScaling", "legacy"
        )
    ),
    reportMomentumStabDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "reportMomentumStabDiagnostics", false
        )
    ),
    momentumStabReportInterval_
    (
        max
        (
            solidModelDict().lookupOrDefault<label>
            (
                "momentumStabReportInterval", 100
            ),
            label(1)
        )
    ),
    writeMomentumStabFields_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writeMomentumStabFields", false
        )
    ),
    momentumStabilisationBoundaryTreatment_
    (
        solidModelDict().lookupOrDefault<word>
        (
            "momentumStabilisationBoundaryTreatment",
            "legacy"
        )
    ),
    faceStressTreatment_
    (
        solidModelDict().lookupOrDefault<word>
        (
            "faceStressTreatment", "interpolatedCell"
        )
    ),
    preconditionerMaterialTangent_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerMaterialTangent", false
        )
    ),
    preconditionerBoundaryTangent_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerBoundaryTangent", false
        )
    ),
    preconditionerViscousTangent_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerViscousTangent", false
        )
    ),
    preconditionerPassiveNominalTangent_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerPassiveNominalTangent", false
        )
    ),
    preconditionerLeastSquaresPressureCoupling_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerLeastSquaresPressureCoupling", false
        )
    ),
    preconditionerLeastSquaresPressureStabilisation_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "preconditionerLeastSquaresPressureStabilisation", false
        )
    ),
    momentumStabDiagnosticCall_(0),
    reportSnesTrialDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "reportSnesTrialDiagnostics", false
        )
    ),
    petscDomainSafeTrials_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "petscDomainSafeTrials", false
        )
    ),
    lastSnesTrialDiagnosticIteration_(-1),
    snesLineSearchTrialCounter_(0),
    writePressureConstraintDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writePressureConstraintDiagnostics", false
        )
    ),
    writeGeometricJacobianDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writeGeometricJacobianDiagnostics", false
        )
    ),
    writeQuasiStaticDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writeQuasiStaticDiagnostics", false
        )
    ),
    reportExtendedIncompressibilityDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "reportExtendedIncompressibilityDiagnostics", false
        )
    ),
    incompressibilityBadCellThreshold_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "incompressibilityBadCellThreshold", 0.05
        )
    ),
    landProblem3DiagnosticsEnabled_
    (
        landProblem3DiagnosticsSwitch(solidModelDict(), "enabled", false)
    ),
    landWriteMomentumDecomposition_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writeMomentumDecomposition",
            true
        )
      : Switch(false)
    ),
    landWritePressureDecomposition_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writePressureDecomposition",
            true
        )
      : Switch(false)
    ),
    landWriteJacobianComparison_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writeJacobianComparison",
            true
        )
      : Switch(false)
    ),
    landWriteBasalAudit_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writeBasalAudit",
            true
        )
      : Switch(false)
    ),
    landWritePressureTractionCheck_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writePressureTractionCheck",
            true
        )
      : Switch(false)
    ),
    landWriteEveryTimeStep_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writeEveryTimeStep",
            false
        )
      : Switch(false)
    ),
    landWriteAtFinalTime_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsSwitch
        (
            solidModelDict(),
            "writeAtFinalTime",
            true
        )
      : Switch(false)
    ),
    landBasalZReference_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsScalar
        (
            solidModelDict(),
            "basalZReference",
            0.005
        )
      : 0.005
    ),
    landBasalZTolerance_
    (
        landProblem3DiagnosticsEnabled_
      ? landProblem3DiagnosticsScalar
        (
            solidModelDict(),
            "basalZTolerance",
            1e-9
        )
      : 1e-9
    ),
    momentumStabilisationConsistencyTestEnabled_
    (
        optionalSubDictSwitch
        (
            solidModelDict(),
            "momentumStabilisationConsistencyTest",
            "enabled",
            false
        )
    ),
    momentumStabilisationConsistencyWriteFields_
    (
        momentumStabilisationConsistencyTestEnabled_
      ? optionalSubDictSwitch
        (
            solidModelDict(),
            "momentumStabilisationConsistencyTest",
            "writeFields",
            true
        )
      : Switch(false)
    ),
    momentumStabilisationConsistencyRunAtStartup_
    (
        momentumStabilisationConsistencyTestEnabled_
      ? optionalSubDictSwitch
        (
            solidModelDict(),
            "momentumStabilisationConsistencyTest",
            "runAtStartup",
            false
        )
      : Switch(false)
    ),
    affineKinematicsConsistencyTestEnabled_
    (
        optionalSubDictSwitch
        (
            solidModelDict(),
            "affineKinematicsConsistencyTest",
            "enabled",
            false
        )
    ),
    affineKinematicsConsistencyWriteFields_
    (
        affineKinematicsConsistencyTestEnabled_
      ? optionalSubDictSwitch
        (
            solidModelDict(),
            "affineKinematicsConsistencyTest",
            "writeFields",
            true
        )
      : Switch(false)
    ),
    affineKinematicsConsistencyRunAtStartup_
    (
        affineKinematicsConsistencyTestEnabled_
      ? optionalSubDictSwitch
        (
            solidModelDict(),
            "affineKinematicsConsistencyTest",
            "runAtStartup",
            false
        )
      : Switch(false)
    ),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    )
{
    DisRequired();

    if
    (
        momentumStabilisationBoundaryTreatment_ != "legacy"
     && momentumStabilisationBoundaryTreatment_ != "internalFacesOnly"
    )
    {
        FatalErrorInFunction
            << "Unknown momentumStabilisationBoundaryTreatment "
            << momentumStabilisationBoundaryTreatment_ << nl
            << "Valid options are legacy and internalFacesOnly"
            << exit(FatalError);
    }

    Info<< "Momentum stabilisation boundary treatment = "
        << momentumStabilisationBoundaryTreatment_ << endl;

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        mechanical().correct(sigma());
    }

    Info<< "solvePressure = " << solvePressure() << endl;
    Info<< "Incompressibility diagnostic controls:" << nl
        << "    writePressureConstraintDiagnostics = "
        << writePressureConstraintDiagnostics_ << nl
        << "    writeGeometricJacobianDiagnostics = "
        << writeGeometricJacobianDiagnostics_ << nl
        << "    writeQuasiStaticDiagnostics = "
        << writeQuasiStaticDiagnostics_ << nl
        << "    reportExtendedIncompressibilityDiagnostics = "
        << reportExtendedIncompressibilityDiagnostics_ << nl
        << "    incompressibilityBadCellThreshold = "
        << incompressibilityBadCellThreshold_ << endl;

    if (landProblem3DiagnosticsEnabled_)
    {
        Info<< "Land Problem 3 diagnostic controls:" << nl
            << "    enabled = " << landProblem3DiagnosticsEnabled_ << nl
            << "    writeMomentumDecomposition = "
            << landWriteMomentumDecomposition_ << nl
            << "    writePressureDecomposition = "
            << landWritePressureDecomposition_ << nl
            << "    writeJacobianComparison = "
            << landWriteJacobianComparison_ << nl
            << "    writeBasalAudit = "
            << landWriteBasalAudit_ << nl
            << "    writePressureTractionCheck = "
            << landWritePressureTractionCheck_ << nl
            << "    writeEveryTimeStep = "
            << landWriteEveryTimeStep_ << nl
            << "    writeAtFinalTime = "
            << landWriteAtFinalTime_ << nl
            << "    basalZReference = "
            << landBasalZReference_ << nl
            << "    basalZTolerance = "
            << landBasalZTolerance_ << endl;
    }

    if (momentumStabilisationConsistencyTestEnabled_)
    {
        Info<< "Momentum stabilisation consistency-test controls:" << nl
            << "    enabled = "
            << momentumStabilisationConsistencyTestEnabled_ << nl
            << "    writeFields = "
            << momentumStabilisationConsistencyWriteFields_ << nl
            << "    runAtStartup = "
            << momentumStabilisationConsistencyRunAtStartup_ << endl;
    }

    if (affineKinematicsConsistencyTestEnabled_)
    {
        Info<< "Affine kinematics consistency-test controls:" << nl
            << "    enabled = "
            << affineKinematicsConsistencyTestEnabled_ << nl
            << "    writeFields = "
            << affineKinematicsConsistencyWriteFields_ << nl
            << "    runAtStartup = "
            << affineKinematicsConsistencyRunAtStartup_ << endl;
    }

    if (solvePressure())
    {
        if (solutionAlg() != solutionAlgorithm::PETSC_SNES)
        {
            FatalErrorInFunction
                << "The solution algorithm must be "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " when solvePressure is enabled" << abort(FatalError);
        }

        // Ensure p is created and the oldTime is stored
        p().oldTime();

        // Initialise dpdt field
        dpdtPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(p())
            )
        );

        // Use a volume-weighted average of 2*mu as the physical scale
        // of the pressure equation. The pressure-row residual and
        // Jacobian are then multiplied by
        // pressureEqnScale_ = pressureScaleFactor_ * twoMuRef_ so that
        // their natural magnitude is comparable to the momentum block.
        const volScalarField twoMu(2.0*mechanical().shearModulus());
        scalar twoMuV = 0;
        scalar Vtot = 0;
        forAll(twoMu, cellI)
        {
            const scalar Vc = mesh().V()[cellI];
            twoMuV += twoMu[cellI]*Vc;
            Vtot += Vc;
        }
        reduce(twoMuV, sumOp<scalar>());
        reduce(Vtot, sumOp<scalar>());
        twoMuRef_ = twoMuV/Vtot;
        pressureEqnScale_ =
            pressureScaleFactor_*(pressureScaleByTwoMu_ ? twoMuRef_ : 1.0);

        pressureUnknownScale_ = 1.0;
        if (scaleMixedPetScFields_)
        {
            if
            (
                pressureUnknownScaleType_ == "twoMu"
             || pressureUnknownScaleType_ == "2mu"
            )
            {
                pressureUnknownScale_ = twoMuRef_;
            }
            else if
            (
                pressureUnknownScaleType_ == "user"
             || pressureUnknownScaleType_ == "scalar"
            )
            {
                pressureUnknownScale_ =
                    readScalar
                    (
                        solidModelDict().lookup("pressureUnknownScaleValue")
                    );
            }
            else if (pressureUnknownScaleType_ == "none")
            {
                pressureUnknownScale_ = 1.0;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown pressureUnknownScale "
                    << pressureUnknownScaleType_ << nl
                    << "Valid options are twoMu, user, scalar, none"
                    << abort(FatalError);
            }

            if (pressureUnknownScale_ <= VSMALL)
            {
                FatalErrorInFunction
                    << "pressureUnknownScale must be positive, found "
                    << pressureUnknownScale_ << abort(FatalError);
            }
        }

        Info<< "pressureEqnScale = " << pressureEqnScale_
            << ", where pressureScaleFactor = " << pressureScaleFactor_
            << " and 2*mu = " << twoMuRef_ << endl;

        if
        (
            pressureRowScaling_ != "legacy"
         && pressureRowScaling_ != "volumeRmsForce"
        )
        {
            FatalErrorInFunction
                << "Unknown pressureRowScaling " << pressureRowScaling_
                << nl << "Valid options are legacy and volumeRmsForce"
                << abort(FatalError);
        }

        Info<< "PETSc pressure row scaling = " << pressureRowScaling_;
        if (pressureRowScaling_ == "volumeRmsForce")
        {
            const scalar L0 = cbrt(Vtot);
            Info<< ": ||R_p||_2 = pressureEqnScale*L0^2"
                << "*volumeRMS(r_p), L0 = " << L0;
        }
        Info<< endl;

        Info<< "PETSc pressure unknown scale = " << pressureUnknownScale_
            << " (scaleMixedPetScFields = " << scaleMixedPetScFields_
            << ", pressureUnknownScale = " << pressureUnknownScaleType_
            << ")" << endl;
    }

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);
        updateJacobianDiagnostics();

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }

    // Check the gradScheme
    const word gradDScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + D().name() +')')
#else
        mesh().schemesDict().gradScheme("grad(" + D().name() +')')
#endif
    );

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        if
        (
            gradDScheme != "leastSquaresS4f"
         && gradDScheme != "leastSquaresS4fDirichlet"
        )
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithm" << abort(FatalError);
        }

        // Set extrapolateValue to true for solidTraction boundaries
        forAll(D().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
                        boundaryFieldRef(D())[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }

    if (momentumStabilisationConsistencyRunAtStartup_)
    {
        runMomentumStabilisationConsistencyTests();
    }

    if (affineKinematicsConsistencyRunAtStartup_)
    {
        runAffineKinematicsConsistencyTest();
    }
}


void nonLinGeomTotalLagTotalDispSolid::makeRKappa() const
{
    if (rKappaPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rKappaPtr_.set(new volScalarField(1.0/mechanical().bulkModulus()));
}


const volScalarField& nonLinGeomTotalLagTotalDispSolid::rKappa() const
{
    if (rKappaPtr_.empty())
    {
        makeRKappa();
    }

    return rKappaPtr_();
}


void nonLinGeomTotalLagTotalDispSolid::updateRAUfIfStale()
{
    if (!solvePressure())
    {
        // rAUf is unused unless we are solving the mixed system
        return;
    }

    const label tIdx = runTime().timeIndex();
    const scalar dt = runTime().deltaT().value();

    // rAUf depends on (impKf_, rho, mesh, deltaT) only. It is fresh
    // when both the timeIndex and deltaT match what we cached
    if
    (
        rAUfTimeIndex_ >= 0
     && rAUfTimeIndex_ == tIdx
     && mag(dt - rAUfDeltaT_) <= SMALL*max(mag(dt), SMALL)
    )
    {
        return;
    }

    // Build the approximate momentum diagonal. fvm::laplacian and
    // fvm::d2dt2 read only the BC structure of D, not its values, so
    // the resulting diagonal is independent of the current Newton
    // iterate and any MFFD perturbation
    fvVectorMatrix approxMomJ
    (
        fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D())
    );
    approxMomJ.relax();
    rAUf() = -1.0/(fvc::interpolate(approxMomJ.A()));

    rAUfTimeIndex_ = tIdx;
    rAUfDeltaT_ = dt;
}


void nonLinGeomTotalLagTotalDispSolid::updateJacobianDiagnostics()
{
    Jminus1_ = J_ - dimensionedScalar("one", dimless, 1.0);
    Jminus1_.correctBoundaryConditions();
}


void nonLinGeomTotalLagTotalDispSolid::updateFinalStateDiagnostics()
{
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    updateJacobianDiagnostics();

    if
    (
        writeGeometricJacobianDiagnostics_
     || reportExtendedIncompressibilityDiagnostics_
     || (runLandProblem3DiagnosticsNow() && landWriteJacobianComparison_)
    )
    {
        updateGeometricJacobianDiagnostics();
    }

    if
    (
        writePressureConstraintDiagnostics_
     || reportExtendedIncompressibilityDiagnostics_
     || (runLandProblem3DiagnosticsNow() && landWritePressureDecomposition_)
    )
    {
        updatePressureConstraintDiagnostics();
    }

    if
    (
        writeQuasiStaticDiagnostics_
     || reportExtendedIncompressibilityDiagnostics_
     || (runLandProblem3DiagnosticsNow() && landWriteMomentumDecomposition_)
    )
    {
        updateQuasiStaticDiagnostics();
    }

    reportIncompressibilityDiagnostics();

    if (reportExtendedIncompressibilityDiagnostics_)
    {
        reportPressureConstraintDiagnostics();
        reportQuasiStaticDiagnostics();
    }
    else if (writeQuasiStaticDiagnostics_)
    {
        reportQuasiStaticDiagnostics();
    }

    if (writePressureConstraintDiagnostics_)
    {
        writePressureConstraintDiagnostics();
    }

    if (writeGeometricJacobianDiagnostics_)
    {
        writeGeometricJacobianDiagnostics();
    }

    if (writeQuasiStaticDiagnostics_)
    {
        writeQuasiStaticDiagnostics();
    }

    updateLandProblem3Diagnostics();

    if (momentumStabilisationConsistencyTestEnabled_)
    {
        if (!momentumStabilisationConsistencyRunAtStartup_)
        {
            runMomentumStabilisationConsistencyTests();
        }
    }

    if (affineKinematicsConsistencyTestEnabled_)
    {
        if (!affineKinematicsConsistencyRunAtStartup_)
        {
            runAffineKinematicsConsistencyTest();
        }
    }
}


void
nonLinGeomTotalLagTotalDispSolid::updatePressureConstraintDiagnostics()
{
    if (!solvePressure())
    {
        return;
    }

    const volVectorField gradp(fvc::grad(p()));

    pressureStabilisation().updateScalar(p(), &gradp);
    updateRAUfIfStale();

    pressureCompressibilityTerm_ = -p()*rKappa();
    pressureStabilisationTerm_ =
        pressureStabilisation().cellScalar(&rAUf(), true);
    volumetricConstraintTerm_ =
        -mixedVolumetricConstraint(J_);
    pressureConstraintResidual_ =
        pressureCompressibilityTerm_
      + pressureStabilisationTerm_
      + volumetricConstraintTerm_;

    pressureCompressibilityTerm_.correctBoundaryConditions();
    pressureStabilisationTerm_.correctBoundaryConditions();
    volumetricConstraintTerm_.correctBoundaryConditions();
    pressureConstraintResidual_.correctBoundaryConditions();
}


void
nonLinGeomTotalLagTotalDispSolid::calcDeformedCellVolumes
(
    const pointVectorField& pointDField,
    scalarField& deformedCellVolumes
) const
{
    pointField deformedPoints(mesh().points());
    const vectorField& pointDI = pointDField.primitiveField();

    forAll(deformedPoints, pointI)
    {
        deformedPoints[pointI] += pointDI[pointI];
    }

    vectorField deformedFaceCentres(mesh().nFaces(), vector::zero);
    vectorField deformedFaceAreas(mesh().nFaces(), vector::zero);
    vectorField deformedCellCentres(mesh().nCells(), vector::zero);

    primitiveMeshTools::makeFaceCentresAndAreas
    (
        mesh(),
        deformedPoints,
        deformedFaceCentres,
        deformedFaceAreas
    );

    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh(),
        deformedFaceCentres,
        deformedFaceAreas,
        deformedCellCentres,
        deformedCellVolumes
    );
}


void
nonLinGeomTotalLagTotalDispSolid::updateGeometricJacobianDiagnostics()
{
    // Jgrad measures the local volume change implied by the reconstructed
    // cell gradient. Jgeom measures the volume change implied by reconstructed
    // point geometry. Disagreement identifies inconsistency between those two
    // kinematic reconstructions; agreement away from one indicates local volume
    // change in both representations. Jgeom is diagnostic only, not inherently
    // more accurate than Jgrad.
    scalarField deformedCellVolumes(mesh().nCells(), 0.0);

    calcDeformedCellVolumes(pointD(), deformedCellVolumes);

    scalarField& JgeomI = primitiveFieldRef(Jgeom_);
    scalarField& JgradMinusJgeomI = primitiveFieldRef(JgradMinusJgeom_);
    scalarField& JgeomMinus1I = primitiveFieldRef(JgeomMinus1_);

    scalarField& highJgrad = primitiveFieldRef(highJgradMask_);
    scalarField& lowJgrad = primitiveFieldRef(lowJgradMask_);
    scalarField& highJgeom = primitiveFieldRef(highJgeomMask_);
    scalarField& lowJgeom = primitiveFieldRef(lowJgeomMask_);
    scalarField& largeDiff = primitiveFieldRef(largeJgradMinusJgeomMask_);

    const scalarField& V = mesh().V();
    const scalar threshold = incompressibilityBadCellThreshold_;
    const globalIndex globalCells(mesh().nCells());

    label nNonPositiveLocal = 0;
    const label maxWarningsPerProcessor = 20;

    forAll(JgeomI, cellI)
    {
        JgeomI[cellI] = deformedCellVolumes[cellI]/V[cellI];
        JgradMinusJgeomI[cellI] = J_[cellI] - JgeomI[cellI];
        JgeomMinus1I[cellI] = JgeomI[cellI] - 1.0;

        highJgrad[cellI] = (J_[cellI] - 1.0 > threshold) ? 1.0 : 0.0;
        lowJgrad[cellI] = (J_[cellI] - 1.0 < -threshold) ? 1.0 : 0.0;
        highJgeom[cellI] = (JgeomI[cellI] - 1.0 > threshold) ? 1.0 : 0.0;
        lowJgeom[cellI] = (JgeomI[cellI] - 1.0 < -threshold) ? 1.0 : 0.0;
        largeDiff[cellI] =
            (mag(JgradMinusJgeomI[cellI]) > threshold) ? 1.0 : 0.0;

        if (deformedCellVolumes[cellI] <= 0.0)
        {
            if (nNonPositiveLocal < maxWarningsPerProcessor)
            {
                Pout<< "Warning: non-positive deformed geometric cell volume"
                    << nl
                    << "    processor = " << Pstream::myProcNo() << nl
                    << "    local cell ID = " << cellI << nl
                    << "    global cell ID = "
                    << globalCells.toGlobal(cellI) << nl
                    << "    reference volume = " << V[cellI] << nl
                    << "    deformed geometric volume = "
                    << deformedCellVolumes[cellI] << nl
                    << "    Jgrad = " << J_[cellI] << nl
                    << "    Jgeom = " << JgeomI[cellI] << endl;
            }

            nNonPositiveLocal++;
        }
    }

    label nNonPositive = nNonPositiveLocal;
    reduce(nNonPositive, sumOp<label>());

    if (nNonPositive > 0)
    {
        WarningInFunction
            << "Detected " << nNonPositive
            << " cells with non-positive deformed geometric volume. "
            << "Printed at most " << maxWarningsPerProcessor
            << " entries per processor." << endl;
    }

    Jgeom_.correctBoundaryConditions();
    JgradMinusJgeom_.correctBoundaryConditions();
    JgeomMinus1_.correctBoundaryConditions();
    highJgradMask_.correctBoundaryConditions();
    lowJgradMask_.correctBoundaryConditions();
    highJgeomMask_.correctBoundaryConditions();
    lowJgeomMask_.correctBoundaryConditions();
    largeJgradMinusJgeomMask_.correctBoundaryConditions();
}


void
nonLinGeomTotalLagTotalDispSolid::updateQuasiStaticDiagnostics()
{
    inertiaForceDensity_ = -rho()*fvc::d2dt2(D());
    dampingForceDensity_ = -rho()*dampingCoeff()*fvc::ddt(D());
    if (!useBodyForceField_)
    {
        bodyForceDensity_ = rho()*g();
    }

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh().Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    const surfaceVectorField tractionPhysical
    (
        nCurrent & fvc::interpolate(sigma())
    );

    surfaceVectorField forcePhysical(magSfCurrent*tractionPhysical);
    enforceTractionBoundaries(forcePhysical, D(), nCurrent, magSfCurrent);
    physicalStressForceDensity_ = fvc::div(forcePhysical);

    momentumStabilisation().updateVector(D(), &gradD());

    const surfaceVectorField tractionWithStab
    (
        tractionPhysical + impKf_*momentumStabilisation().faceVector()
    );

    surfaceVectorField forceWithStab(magSfCurrent*tractionWithStab);
    enforceTractionBoundaries(forceWithStab, D(), nCurrent, magSfCurrent);

    surfaceVectorField forceStab(forceWithStab - forcePhysical);
    applyMomentumStabilisationBoundaryTreatment(forceStab);

    vectorField stabCellForce(fvc::div(forceStab));
    stabCellForce *= mesh().V();

    momentumStabilisation().projectExtensiveVectorForce(stabCellForce);

    primitiveFieldRef(momentumStabForceDensity_) = stabCellForce/mesh().V();

    inertiaForceDensity_.correctBoundaryConditions();
    dampingForceDensity_.correctBoundaryConditions();
    bodyForceDensity_.correctBoundaryConditions();
    physicalStressForceDensity_.correctBoundaryConditions();
    momentumStabForceDensity_.correctBoundaryConditions();
}


void nonLinGeomTotalLagTotalDispSolid::writePressureConstraintDiagnostics() const
{
    if (!solvePressure() || !mesh().time().writeTime())
    {
        return;
    }

    pressureCompressibilityTerm_.write();
    pressureStabilisationTerm_.write();
    volumetricConstraintTerm_.write();
    pressureConstraintResidual_.write();
}


void nonLinGeomTotalLagTotalDispSolid::writeGeometricJacobianDiagnostics() const
{
    if (!mesh().time().writeTime())
    {
        return;
    }

    Jgeom_.write();
    JgradMinusJgeom_.write();
    JgeomMinus1_.write();
    highJgradMask_.write();
    lowJgradMask_.write();
    highJgeomMask_.write();
    lowJgeomMask_.write();
    largeJgradMinusJgeomMask_.write();
}


void nonLinGeomTotalLagTotalDispSolid::writeQuasiStaticDiagnostics() const
{
    if (!mesh().time().writeTime())
    {
        return;
    }

    inertiaForceDensity_.write();
    dampingForceDensity_.write();
    bodyForceDensity_.write();
    physicalStressForceDensity_.write();
    momentumStabForceDensity_.write();
}


void nonLinGeomTotalLagTotalDispSolid::reportPressureConstraintDiagnostics()
const
{
    if (!solvePressure())
    {
        return;
    }

    const scalarField& V = mesh().V();

    const ScalarDiagnosticStats compressibilityStats =
        scalarStats(pressureCompressibilityTerm_, V);
    const ScalarDiagnosticStats stabilisationStats =
        scalarStats(pressureStabilisationTerm_, V);
    const ScalarDiagnosticStats volumetricStats =
        scalarStats(volumetricConstraintTerm_, V);
    const ScalarDiagnosticStats residualStats =
        scalarStats(pressureConstraintResidual_, V);

    scalar volume = 0;
    scalar cancellationSqr = 0;
    scalar stabSqr = 0;
    scalar volSqr = 0;
    scalar corrNumerator = 0;
    scalar componentErrorMax = 0;
    scalar scaledComponentErrorMax = 0;

    forAll(pressureConstraintResidual_, cellI)
    {
        const scalar Vc = V[cellI];
        const scalar stab = pressureStabilisationTerm_[cellI];
        const scalar vol = volumetricConstraintTerm_[cellI];
        const scalar sumComponents =
            pressureCompressibilityTerm_[cellI] + stab + vol;
        const scalar componentError =
            pressureConstraintResidual_[cellI] - sumComponents;

        volume += Vc;
        cancellationSqr += Vc*sqr(stab + vol);
        stabSqr += Vc*sqr(stab);
        volSqr += Vc*sqr(vol);
        corrNumerator += Vc*stab*(-vol);
        componentErrorMax = max(componentErrorMax, mag(componentError));
        scaledComponentErrorMax =
            max
            (
                scaledComponentErrorMax,
                mag(Vc*pressureEqnScale_*componentError)
            );
    }

    reduce(volume, sumOp<scalar>());
    reduce(cancellationSqr, sumOp<scalar>());
    reduce(stabSqr, sumOp<scalar>());
    reduce(volSqr, sumOp<scalar>());
    reduce(corrNumerator, sumOp<scalar>());
    reduce(componentErrorMax, maxOp<scalar>());
    reduce(scaledComponentErrorMax, maxOp<scalar>());

    const scalar cancellationRms =
        sqrt(cancellationSqr/(volume + VSMALL));
    const scalar stabilisationToVolumetricRms =
        stabilisationStats.rms/(volumetricStats.rms + VSMALL);
    const scalar correlation =
        corrNumerator/(sqrt(stabSqr*volSqr) + VSMALL);

    Info<< "Pressure constraint diagnostics:" << nl;
    printScalarStats("pressureCompressibilityTerm", compressibilityStats);
    printScalarStats("pressureStabilisationTerm", stabilisationStats);
    printScalarStats("volumetricConstraintTerm", volumetricStats);
    printScalarStats("pressureConstraintResidual", residualStats);
    Info<< "    RMS(pressureStabilisationTerm"
        << " + volumetricConstraintTerm) = "
        << cancellationRms << nl
        << "    RMS(pressureStabilisationTerm)"
        << "/RMS(volumetricConstraintTerm) = "
        << stabilisationToVolumetricRms << nl
        << "    corr(pressureStabilisationTerm,"
        << " -volumetricConstraintTerm) = "
        << correlation << nl
        << "    max component-sum error = "
        << componentErrorMax << nl
        << "    max scaled extensive component-sum error = "
        << scaledComponentErrorMax << endl;
}


void nonLinGeomTotalLagTotalDispSolid::reportQuasiStaticDiagnostics() const
{
    const scalarField& V = mesh().V();

    const scalar inertiaRms = vectorRms(inertiaForceDensity_, V);
    const scalar dampingRms = vectorRms(dampingForceDensity_, V);
    const scalar physicalStressRms = vectorRms(physicalStressForceDensity_, V);
    const scalar inertiaPlusDampingRms =
        vectorSumRms(inertiaForceDensity_, dampingForceDensity_, V);

    Info<< "Quasi-static diagnostics:" << nl
        << "    RMS(inertiaForceDensity) = " << inertiaRms << nl
        << "    RMS(dampingForceDensity) = " << dampingRms << nl
        << "    RMS(physicalStressForceDensity) = "
        << physicalStressRms << nl
        << "    RMS(inertiaForceDensity)"
        << "/RMS(physicalStressForceDensity) = "
        << inertiaRms/(physicalStressRms + VSMALL) << nl
        << "    RMS(dampingForceDensity)"
        << "/RMS(physicalStressForceDensity) = "
        << dampingRms/(physicalStressRms + VSMALL) << nl
        << "    RMS(inertiaForceDensity + dampingForceDensity)"
        << "/RMS(physicalStressForceDensity) = "
        << inertiaPlusDampingRms/(physicalStressRms + VSMALL) << endl;
}


Foam::tmp<Foam::volSymmTensorField>
nonLinGeomTotalLagTotalDispSolid::activeCauchyStressDiagnostic() const
{
    const PtrList<mechanicalLaw>& laws = mechanical();

    tmp<volSymmTensorField> tSigmaActive
    (
        new volSymmTensorField
        (
            IOobject
            (
                "landSigmaActiveDiagnostic",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        )
    );
    volSymmTensorField& sigmaActive = tmpRef(tSigmaActive);

    if (laws.size() == 1)
    {
        if (isA<electroMechanicalLaw>(laws[0]))
        {
            const electroMechanicalLaw& activeLaw =
                refCast<const electroMechanicalLaw>(laws[0]);

            if (activeLaw.hasActiveStress())
            {
                sigmaActive = activeLaw.activeCauchyStress();
            }
        }

        return tSigmaActive;
    }

    PtrList<volSymmTensorField> subMeshSigmaActive(laws.size());

    forAll(laws, lawI)
    {
        const fvMesh& lawMesh =
            mechanical().solSubMeshes().subMeshes()[lawI].subMesh();

        if (isA<electroMechanicalLaw>(laws[lawI]))
        {
            const electroMechanicalLaw& activeLaw =
                refCast<const electroMechanicalLaw>(laws[lawI]);

            if (activeLaw.hasActiveStress())
            {
                const tmp<volSymmTensorField> tSubSigmaActive =
                    activeLaw.activeCauchyStress();

                subMeshSigmaActive.set
                (
                    lawI,
                    new volSymmTensorField(tSubSigmaActive())
                );
                continue;
            }
        }

        subMeshSigmaActive.set
        (
            lawI,
            new volSymmTensorField
            (
                IOobject
                (
                    "landSubMeshSigmaActiveDiagnostic",
                    lawMesh.time().timeName(),
                    lawMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lawMesh,
                dimensionedSymmTensor
                (
                    "zero",
                    dimPressure,
                    symmTensor::zero
                )
            )
        );
    }

    mechanical().solSubMeshes().mapSubMeshVolFields<symmTensor>
    (
        subMeshSigmaActive,
        sigmaActive
    );

    return tSigmaActive;
}


void nonLinGeomTotalLagTotalDispSolid::enforceDecomposedTractionBoundaries
(
    surfaceVectorField& passiveForce,
    surfaceVectorField& activeForce,
    surfaceVectorField& pressureForce,
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    forAll(D.boundaryField(), patchI)
    {
        vectorField& passiveP = boundaryFieldRef(passiveForce)[patchI];
        vectorField& activeP = boundaryFieldRef(activeForce)[patchI];
        vectorField& pressureP = boundaryFieldRef(pressureForce)[patchI];

        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            passiveP = vector::zero;
            activeP = vector::zero;

            if (tracPatch.useUndeformedArea())
            {
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                pressureP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                pressureP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfCurrentPatch;
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<slipFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            passiveP = sqr(nPatch) & passiveForce.boundaryField()[patchI];
            activeP = sqr(nPatch) & activeForce.boundaryField()[patchI];
            pressureP = sqr(nPatch) & pressureForce.boundaryField()[patchI];
        }
    }
}


void nonLinGeomTotalLagTotalDispSolid::
applyMomentumStabilisationBoundaryTreatment
(
    surfaceVectorField& stabilisationFaceForce
) const
{
    if (momentumStabilisationBoundaryTreatment_ == "legacy")
    {
        return;
    }

    if (momentumStabilisationBoundaryTreatment_ == "internalFacesOnly")
    {
        // The difference-stencil/Rhie-Chow momentum term is a neighbour-cell
        // stabilisation.  Non-coupled physical boundaries have no opposite
        // cell mode to stabilise; fixed-boundary reactions should therefore
        // come from physical stress, not from this one-sided numerical term.
        forAll(stabilisationFaceForce.boundaryField(), patchI)
        {
            if (!stabilisationFaceForce.mesh().boundary()[patchI].coupled())
            {
                stabilisationFaceForce.boundaryFieldRef()[patchI] =
                    vector::zero;
            }
        }

        return;
    }

    FatalErrorInFunction
        << "Unknown momentumStabilisationBoundaryTreatment "
        << momentumStabilisationBoundaryTreatment_ << exit(FatalError);
}


bool nonLinGeomTotalLagTotalDispSolid::runLandProblem3DiagnosticsNow() const
{
    if (!landProblem3DiagnosticsEnabled_)
    {
        return false;
    }

    if (landWriteEveryTimeStep_)
    {
        return true;
    }

    return
        landWriteAtFinalTime_
     && mag(runTime().value() - runTime().endTime().value()) < SMALL;
}


void nonLinGeomTotalLagTotalDispSolid::updateLandProblem3Diagnostics()
{
    if (!runLandProblem3DiagnosticsNow())
    {
        return;
    }

    Info<< "Land Problem 3 operator diagnostics at t = "
        << runTime().timeName() << nl
        << "    Production momentum residual: "
        << "V*(div(boundary-corrected physical face force)"
        << " + projected momentum stabilisation/V"
        << " + rho*(g - d2dt2(D) - dampingCoeff*ddt(D)))" << nl
        << "    Production pressure residual: "
        << "pressureEqnScale*V*(-p*rKappa"
        << " + pressureStabilisation.cellScalar(rAUf)"
        << " - " << mixedVolumetricConstraintDescription() << ")" << nl
        << "    JCellGrad = det(I + grad(D)^T)" << nl
        << "    mesh motion = static total-Lagrangian mesh; "
        << "deformed geometry uses x = X + pointD" << endl;

    if (landWriteMomentumDecomposition_)
    {
        updateLandMomentumDecompositionDiagnostics();
    }

    if (landWritePressureDecomposition_)
    {
        updateLandPressureDecompositionDiagnostics();
    }

    if (landWriteJacobianComparison_)
    {
        updateLandJacobianComparisonDiagnostics();
    }

    if (landWriteBasalAudit_)
    {
        updateLandBasalAuditDiagnostics();
    }

    if (landWritePressureTractionCheck_)
    {
        updateLandPressureTractionDiagnostics();
    }
}


void
nonLinGeomTotalLagTotalDispSolid::updateLandMomentumDecompositionDiagnostics()
const
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    const tmp<volSymmTensorField> tSigmaActive =
        activeCauchyStressDiagnostic();
    const volSymmTensorField& sigmaActive = tSigmaActive();

    volSymmTensorField sigmaPressure
    (
        IOobject
        (
            "landSigmaPressureDiagnostic",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    if (solvePressure())
    {
        sigmaPressure =
            -p()*dimensionedSymmTensor("I", dimless, I);
    }

    const volSymmTensorField sigmaPassive
    (
        IOobject
        (
            "landSigmaPassiveDiagnostic",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sigma() - sigmaActive - sigmaPressure
    );

    surfaceVectorField passiveStressFaceContribution
    (
        IOobject
        (
            "passiveStressFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSfCurrent*(nCurrent & fvc::interpolate(sigmaPassive))
    );

    surfaceVectorField activeStressFaceContribution
    (
        IOobject
        (
            "activeStressFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSfCurrent*(nCurrent & fvc::interpolate(sigmaActive))
    );

    surfaceVectorField pressureFaceContribution
    (
        IOobject
        (
            "pressureFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSfCurrent*(nCurrent & fvc::interpolate(sigmaPressure))
    );

    enforceDecomposedTractionBoundaries
    (
        passiveStressFaceContribution,
        activeStressFaceContribution,
        pressureFaceContribution,
        D(),
        nCurrent,
        magSfCurrent
    );

    surfaceVectorField forcePhysicalProduction
    (
        magSfCurrent*(nCurrent & fvc::interpolate(sigma()))
    );
    enforceTractionBoundaries
    (
        forcePhysicalProduction,
        D(),
        nCurrent,
        magSfCurrent
    );

    momentumStabilisation().updateVector(D(), &gradD());

    const surfaceVectorField tractionPhysical
    (
        nCurrent & fvc::interpolate(sigma())
    );
    const surfaceVectorField tractionWithStab
    (
        tractionPhysical + impKf_*momentumStabilisation().faceVector()
    );

    surfaceVectorField forceWithStab(magSfCurrent*tractionWithStab);
    enforceTractionBoundaries(forceWithStab, D(), nCurrent, magSfCurrent);

    surfaceVectorField momentumStabilisationFaceContribution
    (
        IOobject
        (
            "momentumStabilisationFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        forceWithStab - forcePhysicalProduction
    );
    applyMomentumStabilisationBoundaryTreatment
    (
        momentumStabilisationFaceContribution
    );

    surfaceVectorField momentumStabilisationInternalFaceContribution
    (
        IOobject
        (
            "momentumStabilisationInternalFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );

    surfaceVectorField momentumStabilisationFixedBoundaryFaceContribution
    (
        IOobject
        (
            "momentumStabilisationFixedBoundaryFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );

    surfaceVectorField momentumStabilisationTractionBoundaryFaceContribution
    (
        IOobject
        (
            "momentumStabilisationTractionBoundaryFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );

    surfaceVectorField momentumStabilisationOtherBoundaryFaceContribution
    (
        IOobject
        (
            "momentumStabilisationOtherBoundaryFaceContribution",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce, vector::zero)
    );

    primitiveFieldRef(momentumStabilisationInternalFaceContribution) =
        primitiveField(momentumStabilisationFaceContribution);

    forAll(mesh.boundary(), patchI)
    {
        const vectorField& sourcePatch =
            momentumStabilisationFaceContribution.boundaryField()[patchI];

        if (mesh.boundary()[patchI].coupled())
        {
            momentumStabilisationInternalFaceContribution
                .boundaryFieldRef()[patchI] = sourcePatch;
        }
        else if (isTractionBoundary(D().boundaryField()[patchI]))
        {
            momentumStabilisationTractionBoundaryFaceContribution
                .boundaryFieldRef()[patchI] = sourcePatch;
        }
        else if (isFixedDisplacementBoundary(D().boundaryField()[patchI]))
        {
            momentumStabilisationFixedBoundaryFaceContribution
                .boundaryFieldRef()[patchI] = sourcePatch;
        }
        else
        {
            momentumStabilisationOtherBoundaryFaceContribution
                .boundaryFieldRef()[patchI] = sourcePatch;
        }
    }

    vectorField stabCellForce(fvc::div(momentumStabilisationFaceContribution));
    stabCellForce *= mesh.V();
    momentumStabilisation().projectExtensiveVectorForce(stabCellForce);

    volVectorField passiveMomentumResidual
    (
        IOobject
        (
            "passiveMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(passiveStressFaceContribution)
    );

    volVectorField activeMomentumResidual
    (
        IOobject
        (
            "activeMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(activeStressFaceContribution)
    );

    volVectorField pressureMomentumResidual
    (
        IOobject
        (
            "pressureMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(pressureFaceContribution)
    );

    volVectorField momentumStabilisationResidual
    (
        IOobject
        (
            "momentumStabilisationResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );
    primitiveFieldRef(momentumStabilisationResidual) =
        stabCellForce/mesh.V();

    volVectorField momentumStabilisationInternalResidual
    (
        IOobject
        (
            "momentumStabilisationInternalResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(momentumStabilisationInternalFaceContribution)
    );

    volVectorField momentumStabilisationFixedBoundaryResidual
    (
        IOobject
        (
            "momentumStabilisationFixedBoundaryResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(momentumStabilisationFixedBoundaryFaceContribution)
    );

    volVectorField momentumStabilisationTractionBoundaryResidual
    (
        IOobject
        (
            "momentumStabilisationTractionBoundaryResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(momentumStabilisationTractionBoundaryFaceContribution)
    );

    volVectorField momentumStabilisationOtherBoundaryResidual
    (
        IOobject
        (
            "momentumStabilisationOtherBoundaryResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(momentumStabilisationOtherBoundaryFaceContribution)
    );

    volVectorField momentumStabilisationTotalResidual
    (
        IOobject
        (
            "momentumStabilisationTotalResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        momentumStabilisationResidual
    );

    const volVectorField momentumStabilisationCategoryResidual
    (
        IOobject
        (
            "momentumStabilisationCategoryResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        momentumStabilisationInternalResidual
      + momentumStabilisationFixedBoundaryResidual
      + momentumStabilisationTractionBoundaryResidual
      + momentumStabilisationOtherBoundaryResidual
    );

    volVectorField totalPhysicalMomentumResidual
    (
        IOobject
        (
            "totalPhysicalMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        passiveMomentumResidual
      + activeMomentumResidual
      + pressureMomentumResidual
    );

    volVectorField totalAssembledMomentumResidual
    (
        IOobject
        (
            "totalAssembledMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        totalPhysicalMomentumResidual
      + momentumStabilisationResidual
      + rho()*(g() - fvc::d2dt2(D()) - dampingCoeff()*fvc::ddt(D()))
    );

    const volVectorField productionSurfaceResidual
    (
        IOobject
        (
            "landProductionSurfaceResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(forcePhysicalProduction)
    );

    const volVectorField productionMomentumResidual
    (
        IOobject
        (
            "landProductionMomentumResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        productionSurfaceResidual
      + momentumStabilisationResidual
      + rho()*(g() - fvc::d2dt2(D()) - dampingCoeff()*fvc::ddt(D()))
    );

    vectorField productionExtensive
    (
        primitiveField(productionMomentumResidual)*mesh.V()
    );
    vectorField decomposedExtensive
    (
        primitiveField(totalAssembledMomentumResidual)*mesh.V()
    );

    scalar maxExtensiveError = 0.0;
    scalar rmsExtensiveErrorSqr = 0.0;
    scalar volume = 0.0;

    forAll(productionExtensive, cellI)
    {
        const vector diff =
            decomposedExtensive[cellI] - productionExtensive[cellI];
        maxExtensiveError = max(maxExtensiveError, mag(diff));
        rmsExtensiveErrorSqr += magSqr(diff);
        volume += mesh.V()[cellI];
    }

    reduce(maxExtensiveError, maxOp<scalar>());
    reduce(rmsExtensiveErrorSqr, sumOp<scalar>());
    reduce(volume, sumOp<scalar>());

    const scalar passiveRms = vectorRms(passiveMomentumResidual, mesh.V());
    const scalar activeRms = vectorRms(activeMomentumResidual, mesh.V());
    const scalar pressureRms = vectorRms(pressureMomentumResidual, mesh.V());
    const scalar stabRms =
        vectorRms(momentumStabilisationResidual, mesh.V());
    const scalar physicalScale =
        passiveRms + activeRms + pressureRms;
    scalar localPhysicalMax = 0.0;

    forAll(passiveMomentumResidual, cellI)
    {
        localPhysicalMax =
            max
            (
                localPhysicalMax,
                mag(passiveMomentumResidual[cellI])
              + mag(activeMomentumResidual[cellI])
              + mag(pressureMomentumResidual[cellI])
            );
    }

    reduce(localPhysicalMax, maxOp<scalar>());

    const scalar ratioEpsilon =
        1e-12*max(localPhysicalMax, VSMALL);

    volScalarField momentumStabToPhysicalRatio
    (
        IOobject
        (
            "momentumStabToPhysicalRatio",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalar nCells = 0.0;
    scalar nAbove001 = 0.0;
    scalar nAbove005 = 0.0;
    scalar nAbove010 = 0.0;
    scalar nAbove050 = 0.0;
    scalar nAbove100 = 0.0;
    scalar maxRatio = 0.0;

    forAll(momentumStabToPhysicalRatio, cellI)
    {
        const scalar denom =
            mag(passiveMomentumResidual[cellI])
          + mag(activeMomentumResidual[cellI])
          + mag(pressureMomentumResidual[cellI])
          + ratioEpsilon;
        const scalar ratio =
            mag(momentumStabilisationResidual[cellI])/denom;

        momentumStabToPhysicalRatio[cellI] = ratio;
        maxRatio = max(maxRatio, ratio);
        nCells += 1.0;
        nAbove001 += (ratio > 0.01) ? 1.0 : 0.0;
        nAbove005 += (ratio > 0.05) ? 1.0 : 0.0;
        nAbove010 += (ratio > 0.10) ? 1.0 : 0.0;
        nAbove050 += (ratio > 0.50) ? 1.0 : 0.0;
        nAbove100 += (ratio > 1.00) ? 1.0 : 0.0;
    }

    reduce(maxRatio, maxOp<scalar>());
    reduce(nCells, sumOp<scalar>());
    reduce(nAbove001, sumOp<scalar>());
    reduce(nAbove005, sumOp<scalar>());
    reduce(nAbove010, sumOp<scalar>());
    reduce(nAbove050, sumOp<scalar>());
    reduce(nAbove100, sumOp<scalar>());

    scalar stabCategoryMaxError = 0.0;
    scalar stabCategoryRmsErrorSqr = 0.0;
    scalar stabCategoryCells = 0.0;

    forAll(momentumStabilisationTotalResidual, cellI)
    {
        const vector diff =
            (
                momentumStabilisationCategoryResidual[cellI]
              - momentumStabilisationTotalResidual[cellI]
            )*mesh.V()[cellI];

        stabCategoryMaxError = max(stabCategoryMaxError, mag(diff));
        stabCategoryRmsErrorSqr += magSqr(diff);
        stabCategoryCells += 1.0;
    }

    reduce(stabCategoryMaxError, maxOp<scalar>());
    reduce(stabCategoryRmsErrorSqr, sumOp<scalar>());
    reduce(stabCategoryCells, sumOp<scalar>());

    const scalar stabInternalRms =
        vectorRms(momentumStabilisationInternalResidual, mesh.V());
    const scalar stabFixedRms =
        vectorRms(momentumStabilisationFixedBoundaryResidual, mesh.V());
    const scalar stabTractionRms =
        vectorRms(momentumStabilisationTractionBoundaryResidual, mesh.V());
    const scalar stabOtherRms =
        vectorRms(momentumStabilisationOtherBoundaryResidual, mesh.V());

    Info<< "Land momentum decomposition diagnostics:" << nl
        << "    max extensive residual reconstruction error = "
        << maxExtensiveError << nl
        << "    RMS extensive residual reconstruction error = "
        << sqrt(rmsExtensiveErrorSqr/(nCells + VSMALL)) << nl
        << "    RMS passive/active/pressure/stabilisation force density = "
        << passiveRms << " / " << activeRms << " / "
        << pressureRms << " / " << stabRms << nl
        << "    RMS stabilisation internal/fixed/traction/other = "
        << stabInternalRms << " / " << stabFixedRms << " / "
        << stabTractionRms << " / " << stabOtherRms << nl
        << "    max/RMS extensive stabilisation category closure error = "
        << stabCategoryMaxError << " / "
        << sqrt(stabCategoryRmsErrorSqr/(stabCategoryCells + VSMALL))
        << nl
        << "    global RMS stabilisation/physical ratio = "
        << stabRms/(physicalScale + ratioEpsilon) << nl
        << "    local ratio epsilon = " << ratioEpsilon << nl
        << "    max local momentumStabToPhysicalRatio = "
        << maxRatio << nl
        << "    cell fractions ratio > 0.01/0.05/0.10/0.50/1.00 = "
        << nAbove001/(nCells + VSMALL) << " / "
        << nAbove005/(nCells + VSMALL) << " / "
        << nAbove010/(nCells + VSMALL) << " / "
        << nAbove050/(nCells + VSMALL) << " / "
        << nAbove100/(nCells + VSMALL) << endl;

    passiveStressFaceContribution.write();
    activeStressFaceContribution.write();
    pressureFaceContribution.write();
    momentumStabilisationFaceContribution.write();
    momentumStabilisationInternalFaceContribution.write();
    momentumStabilisationFixedBoundaryFaceContribution.write();
    momentumStabilisationTractionBoundaryFaceContribution.write();
    momentumStabilisationOtherBoundaryFaceContribution.write();
    passiveMomentumResidual.write();
    activeMomentumResidual.write();
    pressureMomentumResidual.write();
    momentumStabilisationResidual.write();
    momentumStabilisationInternalResidual.write();
    momentumStabilisationFixedBoundaryResidual.write();
    momentumStabilisationTractionBoundaryResidual.write();
    momentumStabilisationOtherBoundaryResidual.write();
    momentumStabilisationTotalResidual.write();
    totalPhysicalMomentumResidual.write();
    totalAssembledMomentumResidual.write();
    momentumStabToPhysicalRatio.write();
}


void
nonLinGeomTotalLagTotalDispSolid::updateLandPressureDecompositionDiagnostics()
const
{
    if (!solvePressure())
    {
        Info<< "Land pressure decomposition diagnostics skipped: "
            << "solvePressure is false" << endl;
        return;
    }

    const fvMesh& mesh = this->mesh();

    volScalarField pressureCompressibilityResidual
    (
        IOobject
        (
            "pressureCompressibilityResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureCompressibilityTerm_
    );

    volScalarField pressureGeometricConstraintResidual
    (
        IOobject
        (
            "pressureGeometricConstraintResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        volumetricConstraintTerm_
    );

    volScalarField pressureStabilisationResidual
    (
        IOobject
        (
            "pressureStabilisationResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureStabilisationTerm_
    );

    volScalarField totalPressureResidual
    (
        IOobject
        (
            "totalPressureResidual",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureConstraintResidual_
    );

    forAll(totalPressureResidual, cellI)
    {
        const scalar scale = mesh.V()[cellI]*pressureEqnScale_;
        pressureCompressibilityResidual[cellI] *= scale;
        pressureGeometricConstraintResidual[cellI] *= scale;
        pressureStabilisationResidual[cellI] *= scale;
        totalPressureResidual[cellI] *= scale;
    }

    volScalarField pressureStabToConstraintRatio
    (
        IOobject
        (
            "pressureStabToConstraintRatio",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalar volume = 0.0;
    scalar maxError = 0.0;
    scalar physicalSqr = 0.0;
    scalar stabSqr = 0.0;
    scalar totalSqr = 0.0;
    scalar constraintMax = 0.0;
    scalar corrNumerator = 0.0;
    scalar nCells = 0.0;
    scalar nStabExceedsPhysical = 0.0;
    scalar nOppositeSigns = 0.0;
    scalar cancellationSqr = 0.0;

    forAll(totalPressureResidual, cellI)
    {
        const scalar Vc = mesh.V()[cellI];
        const scalar compressibility =
            pressureCompressibilityResidual[cellI];
        const scalar geometric =
            pressureGeometricConstraintResidual[cellI];
        const scalar physical = compressibility + geometric;
        const scalar stab = pressureStabilisationResidual[cellI];
        const scalar total = totalPressureResidual[cellI];
        const scalar error = total - (compressibility + geometric + stab);

        constraintMax =
            max(constraintMax, mag(compressibility) + mag(geometric));

        pressureStabToConstraintRatio[cellI] = mag(stab);

        volume += Vc;
        maxError = max(maxError, mag(error));
        physicalSqr += Vc*sqr(physical);
        stabSqr += Vc*sqr(stab);
        totalSqr += Vc*sqr(total);
        corrNumerator += Vc*stab*(-physical);
        cancellationSqr += Vc*sqr(physical + stab);
        nCells += 1.0;
        nStabExceedsPhysical += (mag(stab) > mag(physical)) ? 1.0 : 0.0;
        nOppositeSigns += (stab*physical < 0.0) ? 1.0 : 0.0;
    }

    reduce(volume, sumOp<scalar>());
    reduce(maxError, maxOp<scalar>());
    reduce(physicalSqr, sumOp<scalar>());
    reduce(stabSqr, sumOp<scalar>());
    reduce(totalSqr, sumOp<scalar>());
    reduce(constraintMax, maxOp<scalar>());
    reduce(corrNumerator, sumOp<scalar>());
    reduce(cancellationSqr, sumOp<scalar>());
    reduce(nCells, sumOp<scalar>());
    reduce(nStabExceedsPhysical, sumOp<scalar>());
    reduce(nOppositeSigns, sumOp<scalar>());

    const scalar physicalRms = sqrt(physicalSqr/(volume + VSMALL));
    const scalar stabRms = sqrt(stabSqr/(volume + VSMALL));
    const scalar totalRms = sqrt(totalSqr/(volume + VSMALL));
    const scalar cancellationRms =
        sqrt(cancellationSqr/(volume + VSMALL));
    const scalar correlation =
        corrNumerator/(sqrt(stabSqr*physicalSqr) + VSMALL);
    const scalar ratioEpsilon =
        1e-12*max(constraintMax, VSMALL);

    forAll(pressureStabToConstraintRatio, cellI)
    {
        const scalar physical =
            pressureCompressibilityResidual[cellI]
          + pressureGeometricConstraintResidual[cellI];

        pressureStabToConstraintRatio[cellI] /=
            mag(physical) + ratioEpsilon;
    }

    Info<< "Land pressure decomposition diagnostics:" << nl
        << "    compressibility term = "
        << "pressureEqnScale*V*(-p*rKappa)" << nl
        << "    geometric constraint term = "
        << "pressureEqnScale*V*(-"
        << mixedVolumetricConstraintDescription() << ")" << nl
        << "    stabilisation term = "
        << "pressureEqnScale*V*pressureStabilisation.cellScalar(rAUf)"
        << nl
        << "    max component-sum error = " << maxError << nl
        << "    RMS physical/stabilisation/total = "
        << physicalRms << " / " << stabRms << " / " << totalRms << nl
        << "    RMS(physical + stabilisation) = "
        << cancellationRms << nl
        << "    corr(stabilisation, -physical) = "
        << correlation << nl
        << "    local ratio epsilon = " << ratioEpsilon << nl
        << "    fraction |stabilisation| > |physical| = "
        << nStabExceedsPhysical/(nCells + VSMALL) << nl
        << "    fraction opposite signs = "
        << nOppositeSigns/(nCells + VSMALL) << endl;

    pressureCompressibilityResidual.write();
    pressureGeometricConstraintResidual.write();
    pressureStabilisationResidual.write();
    totalPressureResidual.write();
    pressureStabToConstraintRatio.write();
}


void
nonLinGeomTotalLagTotalDispSolid::updateLandJacobianComparisonDiagnostics()
const
{
    const fvMesh& mesh = this->mesh();

    surfaceTensorField gradDf
    (
        IOobject
        (
            "landGradDf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    const_cast<mechanicalModel&>(mechanical()).grad(D(), pointD(), gradDf);

    surfaceTensorField diagnosticFfFromReconstructedGradD
    (
        IOobject
        (
            "diagnosticFfFromReconstructedGradD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDf.T()
    );

    surfaceScalarField JFaceFromReconstructedGradD
    (
        IOobject
        (
            "JFaceFromReconstructedGradD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(diagnosticFfFromReconstructedGradD)
    );

    surfaceScalarField JInterpolatedCellToFace
    (
        IOobject
        (
            "JInterpolatedCellToFace",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(J_)
    );

    volScalarField JFaceInternalOnlyCellAverage
    (
        IOobject
        (
            "JFaceInternalOnlyCellAverage",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    scalarField faceWeightSum(mesh.nCells(), 0.0);
    scalarField JFaceSum(mesh.nCells(), 0.0);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& magSf = mesh.magSf();

    forAll(neighbour, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];
        const scalar w = magSf[faceI];
        JFaceSum[own] += w*JFaceFromReconstructedGradD[faceI];
        faceWeightSum[own] += w;
        JFaceSum[nei] += w*JFaceFromReconstructedGradD[faceI];
        faceWeightSum[nei] += w;
    }

    forAll(mesh.boundary(), patchI)
    {
        const word& patchName = mesh.boundary()[patchI].name();

        if
        (
            patchName == "fixed"
         || patchName == "inside"
         || patchName == "outside"
        )
        {
            scalar minJ = VGREAT;
            scalar maxJ = -VGREAT;
            scalar sumJ = 0.0;
            scalar sumSqrJ = 0.0;
            scalar nFaces = 0.0;

            const scalarField& patchJ =
                JFaceFromReconstructedGradD.boundaryField()[patchI];

            forAll(patchJ, faceI)
            {
                const scalar J = patchJ[faceI];

                minJ = min(minJ, J);
                maxJ = max(maxJ, J);
                sumJ += J;
                sumSqrJ += sqr(J);
                nFaces += 1.0;
            }

            reduce(minJ, minOp<scalar>());
            reduce(maxJ, maxOp<scalar>());
            reduce(sumJ, sumOp<scalar>());
            reduce(sumSqrJ, sumOp<scalar>());
            reduce(nFaces, sumOp<scalar>());

            if (nFaces > 0.5)
            {
                Info<< "    boundary patch " << patchName
                    << " JFaceFromReconstructedGradD:"
                    << " min = " << minJ
                    << ", max = " << maxJ
                    << ", mean = " << sumJ/nFaces
                    << ", RMS = " << sqrt(sumSqrJ/nFaces) << nl;
            }
        }
    }

    forAll(JFaceInternalOnlyCellAverage, cellI)
    {
        JFaceInternalOnlyCellAverage[cellI] =
            JFaceSum[cellI]/(faceWeightSum[cellI] + VSMALL);
    }

    volScalarField JCellGradMinusJDeformedVolume
    (
        IOobject
        (
            "JCellGradMinusJDeformedVolume",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JgradMinusJgeom_
    );

    volScalarField JFaceInternalAverageMinusJCellGrad
    (
        IOobject
        (
            "JFaceInternalAverageMinusJCellGrad",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JFaceInternalOnlyCellAverage - J_
    );

    volScalarField JFaceInternalAverageMinusJDeformedVolume
    (
        IOobject
        (
            "JFaceInternalAverageMinusJDeformedVolume",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JFaceInternalOnlyCellAverage - Jgeom_
    );

    const ScalarDiagnosticStats cellMinusGeomStats =
        scalarStats(JCellGradMinusJDeformedVolume, mesh.V());
    const ScalarDiagnosticStats faceMinusCellStats =
        scalarStats(JFaceInternalAverageMinusJCellGrad, mesh.V());
    const ScalarDiagnosticStats faceMinusGeomStats =
        scalarStats(JFaceInternalAverageMinusJDeformedVolume, mesh.V());

    Info<< "Land Jacobian comparison diagnostics:" << nl
        << "    JCellGrad = det(I + grad(D)^T)" << nl
        << "    JFaceFromReconstructedGradD = det(I + grad(D)f^T), "
        << "using mechanical().grad(D,pointD,gradDf)" << nl
        << "    JFaceInternalOnlyCellAverage excludes all boundary faces"
        << nl
        << "    JInterpolatedCellToFace = fvc::interpolate(JCellGrad), "
        << "not an independent deformation measure" << nl
        << "    JDeformedVolume = Vdeformed(pointD)/Vreference" << nl;
    printScalarStats("JCellGradMinusJDeformedVolume", cellMinusGeomStats);
    printScalarStats
    (
        "JFaceInternalAverageMinusJCellGrad",
        faceMinusCellStats
    );
    printScalarStats
    (
        "JFaceInternalAverageMinusJDeformedVolume",
        faceMinusGeomStats
    );

    diagnosticFfFromReconstructedGradD.write();
    JFaceFromReconstructedGradD.write();
    JInterpolatedCellToFace.write();
    JFaceInternalOnlyCellAverage.write();
    JCellGradMinusJDeformedVolume.write();
    JFaceInternalAverageMinusJCellGrad.write();
    JFaceInternalAverageMinusJDeformedVolume.write();
}


void
nonLinGeomTotalLagTotalDispSolid::updateLandBasalAuditDiagnostics() const
{
    const fvMesh& mesh = this->mesh();

    surfaceScalarField fixedBaseFaceMask
    (
        IOobject
        (
            "fixedBaseFaceMask",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    volScalarField fixedBaseOwnerCellMask
    (
        IOobject
        (
            "fixedBaseOwnerCellMask",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    const pointMesh& pMesh = pointMesh::New(mesh);
    pointScalarField fixedBasePointMask
    (
        IOobject
        (
            "fixedBasePointMask",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    boolList usedPoints(mesh.nPoints(), false);
    boolList usedOwnerCells(mesh.nCells(), false);

    label totalFixedFaces = 0;
    label nFixedPatches = 0;
    vector minCoord(VGREAT, VGREAT, VGREAT);
    vector maxCoord(-VGREAT, -VGREAT, -VGREAT);
    scalar minZ = VGREAT;
    scalar maxZ = -VGREAT;
    scalar maxPointD = 0.0;
    scalar maxBoundaryD = 0.0;
    label nFacesOffBasalPlane = 0;

    forAll(D().boundaryField(), patchI)
    {
        const fvPatchVectorField& DPatch = D().boundaryField()[patchI];
        const bool fixedPatch =
            isA<fixedDisplacementFvPatchVectorField>(DPatch)
         || isA<fixedDisplacementZeroShearFvPatchVectorField>(DPatch)
         || DPatch.fixesValue();

        if (!fixedPatch)
        {
            continue;
        }

        nFixedPatches++;

        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const labelUList& faceCells = mesh.boundary()[patchI].faceCells();
        scalarField& faceMaskP =
            fixedBaseFaceMask.boundaryFieldRef()[patchI];

        label patchPointCount = 0;
        boolList patchPoints(mesh.nPoints(), false);

        forAll(faceMaskP, faceI)
        {
            faceMaskP[faceI] = 1.0;
            totalFixedFaces++;

            const label cellI = faceCells[faceI];
            fixedBaseOwnerCellMask[cellI] = 1.0;
            usedOwnerCells[cellI] = true;

            const face& f = mesh.faces()[pp.start() + faceI];
            bool faceOnBasalPlane = true;

            forAll(f, fpI)
            {
                const label pointI = f[fpI];
                const point& X = mesh.points()[pointI];

                minCoord.x() = min(minCoord.x(), X.x());
                minCoord.y() = min(minCoord.y(), X.y());
                minCoord.z() = min(minCoord.z(), X.z());
                maxCoord.x() = max(maxCoord.x(), X.x());
                maxCoord.y() = max(maxCoord.y(), X.y());
                maxCoord.z() = max(maxCoord.z(), X.z());
                minZ = min(minZ, X.z());
                maxZ = max(maxZ, X.z());

                if (mag(X.z() - landBasalZReference_) > landBasalZTolerance_)
                {
                    faceOnBasalPlane = false;
                }

                if (!usedPoints[pointI])
                {
                    usedPoints[pointI] = true;
                    fixedBasePointMask[pointI] = 1.0;
                }

                if (!patchPoints[pointI])
                {
                    patchPoints[pointI] = true;
                    patchPointCount++;
                }

                maxPointD = max(maxPointD, mag(pointD()[pointI]));
            }

            if (!faceOnBasalPlane)
            {
                nFacesOffBasalPlane++;
            }
        }

        forAll(DPatch, faceI)
        {
            maxBoundaryD = max(maxBoundaryD, mag(DPatch[faceI]));
        }

        label patchOwnerCount = 0;
        forAll(faceCells, faceI)
        {
            if (usedOwnerCells[faceCells[faceI]])
            {
                patchOwnerCount++;
            }
        }

        Info<< "Land fixed-base patch audit:" << nl
            << "    patch name = " << pp.name() << nl
            << "    D boundary type = " << DPatch.type() << nl
            << "    D fixesValue() = " << DPatch.fixesValue() << nl
            << "    face count = " << pp.size() << nl
            << "    unique patch point count = " << patchPointCount << nl
            << "    owner-cell face count = " << patchOwnerCount << nl;
    }

    label uniquePointCount = 0;
    label uniqueOwnerCellCount = 0;
    forAll(usedPoints, pointI)
    {
        uniquePointCount += usedPoints[pointI] ? 1 : 0;
    }
    forAll(usedOwnerCells, cellI)
    {
        uniqueOwnerCellCount += usedOwnerCells[cellI] ? 1 : 0;
    }

    reduce(totalFixedFaces, sumOp<label>());
    reduce(nFixedPatches, sumOp<label>());
    reduce(uniquePointCount, sumOp<label>());
    reduce(uniqueOwnerCellCount, sumOp<label>());
    reduce(minZ, minOp<scalar>());
    reduce(maxZ, maxOp<scalar>());
    reduce(maxPointD, maxOp<scalar>());
    reduce(maxBoundaryD, maxOp<scalar>());
    reduce(nFacesOffBasalPlane, sumOp<label>());
    reduce(minCoord.x(), minOp<scalar>());
    reduce(minCoord.y(), minOp<scalar>());
    reduce(minCoord.z(), minOp<scalar>());
    reduce(maxCoord.x(), maxOp<scalar>());
    reduce(maxCoord.y(), maxOp<scalar>());
    reduce(maxCoord.z(), maxOp<scalar>());

    Info<< "Land fixed-base audit summary:" << nl
        << "    fixed patch count = " << nFixedPatches << nl
        << "    fixed face count = " << totalFixedFaces << nl
        << "    unique fixed point count = " << uniquePointCount << nl
        << "    unique owner-cell count = " << uniqueOwnerCellCount << nl
        << "    reference-coordinate min = " << minCoord << nl
        << "    reference-coordinate max = " << maxCoord << nl
        << "    min/max z = " << minZ << " / " << maxZ << nl
        << "    basal z reference/tolerance = "
        << landBasalZReference_ << " / " << landBasalZTolerance_ << nl
        << "    faces outside basal z tolerance = "
        << nFacesOffBasalPlane << nl
        << "    max |D boundary| on fixed patches = "
        << maxBoundaryD << nl
        << "    max |pointD| on fixed-patch points = "
        << maxPointD << nl
        << "    mask fields written: fixedBaseFaceMask, "
        << "fixedBaseOwnerCellMask, fixedBasePointMask" << endl;

    fixedBaseFaceMask.write();
    fixedBaseOwnerCellMask.write();
    fixedBasePointMask.write();
}


void
nonLinGeomTotalLagTotalDispSolid::updateLandPressureTractionDiagnostics()
const
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    pointField deformedPoints(mesh.points());
    const vectorField& pointDI = pointD().primitiveField();

    forAll(deformedPoints, pointI)
    {
        deformedPoints[pointI] += pointDI[pointI];
    }

    vectorField deformedFaceCentres(mesh.nFaces(), vector::zero);
    vectorField deformedFaceAreas(mesh.nFaces(), vector::zero);
    primitiveMeshTools::makeFaceCentresAndAreas
    (
        mesh,
        deformedPoints,
        deformedFaceCentres,
        deformedFaceAreas
    );

    vector basalCentre(vector::zero);
    scalar basalWeight = 0.0;

    forAll(D().boundaryField(), patchI)
    {
        const fvPatchVectorField& DPatch = D().boundaryField()[patchI];
        const bool fixedPatch =
            isA<fixedDisplacementFvPatchVectorField>(DPatch)
         || isA<fixedDisplacementZeroShearFvPatchVectorField>(DPatch)
         || DPatch.fixesValue();

        if (!fixedPatch)
        {
            continue;
        }

        const scalarField& area = mesh.boundary()[patchI].magSf();

        forAll(area, faceI)
        {
            basalCentre += area[faceI]*mesh.Cf().boundaryField()[patchI][faceI];
            basalWeight += area[faceI];
        }
    }

    reduce(basalCentre, sumOp<vector>());
    reduce(basalWeight, sumOp<scalar>());
    basalCentre /= (basalWeight + VSMALL);

    vector solverForce(vector::zero);
    vector independentForce(vector::zero);
    vector solverMomentOrigin(vector::zero);
    vector independentMomentOrigin(vector::zero);
    vector solverMomentBasal(vector::zero);
    vector independentMomentBasal(vector::zero);
    scalar referenceArea = 0.0;
    scalar currentArea = 0.0;
    label faceCount = 0;

    forAll(D().boundaryField(), patchI)
    {
        if
        (
            !isA<solidTractionFvPatchVectorField>
            (
                D().boundaryField()[patchI]
            )
        )
        {
            continue;
        }

        const solidTractionFvPatchVectorField& tracPatch =
            refCast<const solidTractionFvPatchVectorField>
            (
                D().boundaryField()[patchI]
            );

        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const vectorField& nPatch = nCurrent.boundaryField()[patchI];
        const scalarField& magSfCurrentPatch =
            magSfCurrent.boundaryField()[patchI];
        const scalarField& magSfRefPatch = mesh.boundary()[patchI].magSf();
        const scalarField& pressurePatch = tracPatch.pressure();
        const vectorField& tractionPatch = tracPatch.traction();

        vector patchSolverForce(vector::zero);
        vector patchIndependentForce(vector::zero);
        scalar patchReferenceArea = 0.0;
        scalar patchCurrentArea = 0.0;

        forAll(pressurePatch, faceI)
        {
            const label faceID = pp.start() + faceI;
            const scalar solverArea =
                tracPatch.useUndeformedArea()
              ? magSfRefPatch[faceI]
              : magSfCurrentPatch[faceI];

            const vector solverFaceForce =
                (
                    tractionPatch[faceI]
                  - nPatch[faceI]*pressurePatch[faceI]
                )*solverArea;

            const vector independentFaceForce =
                -pressurePatch[faceI]*nPatch[faceI]*magSfCurrentPatch[faceI];

            const vector C = deformedFaceCentres[faceID];

            solverForce += solverFaceForce;
            independentForce += independentFaceForce;
            solverMomentOrigin += C ^ solverFaceForce;
            independentMomentOrigin += C ^ independentFaceForce;
            solverMomentBasal += (C - basalCentre) ^ solverFaceForce;
            independentMomentBasal += (C - basalCentre) ^ independentFaceForce;
            patchSolverForce += solverFaceForce;
            patchIndependentForce += independentFaceForce;
            referenceArea += magSfRefPatch[faceI];
            currentArea += magSfCurrentPatch[faceI];
            patchReferenceArea += magSfRefPatch[faceI];
            patchCurrentArea += magSfCurrentPatch[faceI];
            faceCount++;
        }

        Info<< "Land pressure traction patch check:" << nl
            << "    patch = " << pp.name() << nl
            << "    useUndeformedArea = " << tracPatch.useUndeformedArea()
            << nl
            << "    face count = " << pressurePatch.size() << nl
            << "    reference/current area = "
            << patchReferenceArea << " / " << patchCurrentArea << nl
            << "    solver force = " << patchSolverForce << nl
            << "    independent -p*nCurrent*Acurrent force = "
            << patchIndependentForce << nl;
    }

    reduce(solverForce, sumOp<vector>());
    reduce(independentForce, sumOp<vector>());
    reduce(solverMomentOrigin, sumOp<vector>());
    reduce(independentMomentOrigin, sumOp<vector>());
    reduce(solverMomentBasal, sumOp<vector>());
    reduce(independentMomentBasal, sumOp<vector>());
    reduce(referenceArea, sumOp<scalar>());
    reduce(currentArea, sumOp<scalar>());
    reduce(faceCount, sumOp<label>());

    const vector forceDiff = solverForce - independentForce;
    const vector momentOriginDiff =
        solverMomentOrigin - independentMomentOrigin;
    const vector momentBasalDiff =
        solverMomentBasal - independentMomentBasal;

    Info<< "Land pressure traction force/moment diagnostics:" << nl
        << "    basal moment centre = " << basalCentre << nl
        << "    face count = " << faceCount << nl
        << "    reference endocardial/traction area = "
        << referenceArea << nl
        << "    current endocardial/traction area = "
        << currentArea << nl
        << "    solver force = " << solverForce << nl
        << "    independent force = " << independentForce << nl
        << "    absolute force difference = " << mag(forceDiff) << nl
        << "    relative force difference = "
        << mag(forceDiff)/(mag(independentForce) + VSMALL) << nl
        << "    solver moment about origin = "
        << solverMomentOrigin << nl
        << "    independent moment about origin = "
        << independentMomentOrigin << nl
        << "    absolute origin-moment difference = "
        << mag(momentOriginDiff) << nl
        << "    relative origin-moment difference = "
        << mag(momentOriginDiff)/(mag(independentMomentOrigin) + VSMALL)
        << nl
        << "    solver moment about basal centre = "
        << solverMomentBasal << nl
        << "    independent moment about basal centre = "
        << independentMomentBasal << nl
        << "    absolute basal-moment difference = "
        << mag(momentBasalDiff) << nl
        << "    relative basal-moment difference = "
        << mag(momentBasalDiff)/(mag(independentMomentBasal) + VSMALL)
        << nl
        << "    net force direction = "
        << solverForce/(mag(solverForce) + VSMALL) << nl
        << "    net moment direction about basal centre = "
        << solverMomentBasal/(mag(solverMomentBasal) + VSMALL) << nl
        << "    all solidTraction patch faces were visited once" << endl;
}


void
nonLinGeomTotalLagTotalDispSolid::runMomentumStabilisationConsistencyTests()
const
{
    const fvMesh& mesh = this->mesh();

    Info<< "Momentum stabilisation consistency tests:" << nl
        << "    raw face operator = scaleFactor*(snGrad(D)"
        << " - n & interpolate(gradD))" << nl
        << "    scaled diagnostic face contribution = "
        << "magSfCurrent*impKf*rawFaceOperator" << nl
        << "    production-corrected contribution is isolated after "
        << "traction-boundary enforcement" << nl
        << "    boundary treatment = "
        << momentumStabilisationBoundaryTreatment_ << nl
        << "    cell residual = div(production-corrected contribution)" << nl
        << "    these tests do not alter the production solution" << endl;

    volVectorField prescribedD
    (
        IOobject
        (
            "prescribedD_work",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero),
        diagnosticVectorPatchTypes(mesh, D())
    );

    volTensorField prescribedGradD
    (
        IOobject
        (
            "grad(prescribedD_work)",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero),
        diagnosticTensorPatchTypes(mesh, gradD())
    );

    const vector bConstant(1.1e-4, -2.0e-4, 0.7e-4);

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] = bConstant;
        prescribedGradD[cellI] = tensor::zero;
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];
        DPatch = bConstant;
        gradPatch = tensor::zero;
    }
    runMomentumStabilisationConsistencyTest
    (
        "constant",
        prescribedD,
        prescribedGradD,
        true
    );

    const vector bTranslation(-0.4e-4, 1.6e-4, -1.3e-4);

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] = bTranslation;
        prescribedGradD[cellI] = tensor::zero;
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];
        DPatch = bTranslation;
        gradPatch = tensor::zero;
    }
    runMomentumStabilisationConsistencyTest
    (
        "rigidTranslation",
        prescribedD,
        prescribedGradD,
        true
    );

    const vector omega(0.17, -0.11, 0.23);
    const tensor infinitesimalRotationMap(skewDisplacementMap(omega));
    const tensor infinitesimalRotationGrad
    (
        gradDFromDisplacementMap(infinitesimalRotationMap)
    );

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] =
            infinitesimalRotationMap & mesh.C()[cellI];
        prescribedGradD[cellI] = infinitesimalRotationGrad;
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];

        forAll(DPatch, faceI)
        {
            DPatch[faceI] = infinitesimalRotationMap & CfPatch[faceI];
            gradPatch[faceI] = infinitesimalRotationGrad;
        }
    }
    runMomentumStabilisationConsistencyTest
    (
        "infinitesimalRigidRotation",
        prescribedD,
        prescribedGradD,
        true
    );

    const tensor finiteRotationMap(zRotation(0.35) - I);
    const tensor finiteRotationGrad
    (
        gradDFromDisplacementMap(finiteRotationMap)
    );

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] = finiteRotationMap & mesh.C()[cellI];
        prescribedGradD[cellI] = finiteRotationGrad;
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];

        forAll(DPatch, faceI)
        {
            DPatch[faceI] = finiteRotationMap & CfPatch[faceI];
            gradPatch[faceI] = finiteRotationGrad;
        }
    }
    runMomentumStabilisationConsistencyTest
    (
        "finiteRigidRotation",
        prescribedD,
        prescribedGradD,
        true
    );

    const tensor affineMap
    (
        0.08,  0.11, -0.03,
       -0.04, -0.06,  0.09,
        0.02, -0.05,  0.07
    );
    const vector affineB(0.3e-4, -0.2e-4, 0.5e-4);
    const tensor affineGrad(gradDFromDisplacementMap(affineMap));

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] =
            affineDisplacement(affineMap, affineB, mesh.C()[cellI]);
        prescribedGradD[cellI] = affineGrad;
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];

        forAll(DPatch, faceI)
        {
            DPatch[faceI] =
                affineDisplacement(affineMap, affineB, CfPatch[faceI]);
            gradPatch[faceI] = affineGrad;
        }
    }
    runMomentumStabilisationConsistencyTest
    (
        "generalAffine",
        prescribedD,
        prescribedGradD,
        true
    );

    forAll(prescribedD, cellI)
    {
        prescribedD[cellI] = quadraticDisplacement(mesh.C()[cellI]);
        prescribedGradD[cellI] = quadraticGradD(mesh.C()[cellI]);
    }
    prescribedD.correctBoundaryConditions();
    prescribedGradD.correctBoundaryConditions();
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        vectorField& DPatch = prescribedD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = prescribedGradD.boundaryFieldRef()[patchI];

        forAll(DPatch, faceI)
        {
            DPatch[faceI] = quadraticDisplacement(CfPatch[faceI]);
            gradPatch[faceI] = quadraticGradD(CfPatch[faceI]);
        }
    }
    runMomentumStabilisationConsistencyTest
    (
        "smoothQuadratic",
        prescribedD,
        prescribedGradD,
        false
    );
}


void
nonLinGeomTotalLagTotalDispSolid::runMomentumStabilisationConsistencyTest
(
    const word& testName,
    const volVectorField& prescribedD,
    const volTensorField& prescribedGradD,
    const bool shouldVanish
) const
{
    const fvMesh& mesh = this->mesh();

    momentumStabilisation().updateVector(prescribedD, &prescribedGradD);

    surfaceVectorField rawFaceStabilisation
    (
        IOobject
        (
            "momentumStabilisationRawFace_" + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        momentumStabilisation().faceVector()
    );

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    surfaceVectorField rawFaceContribution
    (
        IOobject
        (
            "rawMomentumStabilisationFaceContribution_" + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSfCurrent*impKf_*rawFaceStabilisation
    );

    const surfaceVectorField tractionPhysical
    (
        nCurrent & fvc::interpolate(sigma())
    );

    surfaceVectorField forcePhysical(magSfCurrent*tractionPhysical);
    enforceTractionBoundaries(forcePhysical, D(), nCurrent, magSfCurrent);

    const surfaceVectorField tractionWithStab
    (
        tractionPhysical + impKf_*rawFaceStabilisation
    );

    surfaceVectorField forceWithStab(magSfCurrent*tractionWithStab);
    enforceTractionBoundaries(forceWithStab, D(), nCurrent, magSfCurrent);

    surfaceVectorField productionCorrectedFaceContribution
    (
        IOobject
        (
            "productionBoundaryCorrectedMomentumStabilisationFaceContribution_"
          + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        forceWithStab - forcePhysical
    );
    applyMomentumStabilisationBoundaryTreatment
    (
        productionCorrectedFaceContribution
    );

    vectorField productionCellForce
    (
        fvc::div(productionCorrectedFaceContribution)
    );
    productionCellForce *= mesh.V();

    momentumStabilisation().projectExtensiveVectorForce(productionCellForce);

    volVectorField residual
    (
        IOobject
        (
            "momentumStabilisationResidual_" + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );
    primitiveFieldRef(residual) = productionCellForce/mesh.V();
    residual.correctBoundaryConditions();

    surfaceVectorField internalFaceContribution
    (
        IOobject
        (
            "momentumStabilisationInternalFaceContribution_" + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        productionCorrectedFaceContribution
    );

    forAll(mesh.boundary(), patchI)
    {
        if (!mesh.boundary()[patchI].coupled())
        {
            internalFaceContribution.boundaryFieldRef()[patchI] =
                vector::zero;
        }
    }

    const volVectorField internalResidual
    (
        IOobject
        (
            "momentumStabilisationInternalResidual_" + testName,
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(internalFaceContribution)
    );

    scalar rawStabInternalSqr = 0.0;
    scalar rawStabInternalMax = 0.0;
    scalar rawForceInternalSqr = 0.0;
    scalar rawForceInternalMax = 0.0;
    scalar productionInternalSqr = 0.0;
    scalar productionInternalMax = 0.0;
    scalar nInternalFaces = 0.0;

    forAll(rawFaceStabilisation, faceI)
    {
        rawStabInternalSqr += magSqr(rawFaceStabilisation[faceI]);
        rawStabInternalMax =
            max(rawStabInternalMax, mag(rawFaceStabilisation[faceI]));
        rawForceInternalSqr += magSqr(rawFaceContribution[faceI]);
        rawForceInternalMax =
            max(rawForceInternalMax, mag(rawFaceContribution[faceI]));
        productionInternalSqr +=
            magSqr(productionCorrectedFaceContribution[faceI]);
        productionInternalMax =
            max
            (
                productionInternalMax,
                mag(productionCorrectedFaceContribution[faceI])
            );
        nInternalFaces += 1.0;
    }

    scalar rawStabBoundarySqr = 0.0;
    scalar rawStabBoundaryMax = 0.0;
    scalar rawForceBoundarySqr = 0.0;
    scalar rawForceBoundaryMax = 0.0;
    scalar productionBoundarySqr = 0.0;
    scalar productionBoundaryMax = 0.0;
    scalar productionFixedBoundarySqr = 0.0;
    scalar productionFixedBoundaryMax = 0.0;
    scalar nFixedBoundaryFaces = 0.0;
    scalar productionTractionBoundarySqr = 0.0;
    scalar productionTractionBoundaryMax = 0.0;
    scalar nTractionBoundaryFaces = 0.0;
    scalar productionOtherBoundarySqr = 0.0;
    scalar productionOtherBoundaryMax = 0.0;
    scalar nOtherBoundaryFaces = 0.0;
    scalar nBoundaryFaces = 0.0;
    vector rawBoundaryResultant(vector::zero);
    vector productionBoundaryResultant(vector::zero);

    forAll(mesh.boundary(), patchI)
    {
        const vectorField& rawPatch =
            rawFaceStabilisation.boundaryField()[patchI];
        const vectorField& rawForcePatch =
            rawFaceContribution.boundaryField()[patchI];
        const vectorField& productionPatch =
            productionCorrectedFaceContribution.boundaryField()[patchI];
        const bool fixedPatch =
            !mesh.boundary()[patchI].coupled()
         && isFixedDisplacementBoundary(D().boundaryField()[patchI]);
        const bool tractionPatch =
            !mesh.boundary()[patchI].coupled()
         && isTractionBoundary(D().boundaryField()[patchI]);
        const bool otherPatch =
            !mesh.boundary()[patchI].coupled()
         && !fixedPatch
         && !tractionPatch;

        vector patchRawResultant(vector::zero);
        vector patchProductionResultant(vector::zero);
        scalar patchRawSqr = 0.0;
        scalar patchRawMax = 0.0;
        scalar patchProductionSqr = 0.0;
        scalar patchProductionMax = 0.0;
        scalar nPatchFaces = 0.0;

        forAll(productionPatch, faceI)
        {
            rawStabBoundarySqr += magSqr(rawPatch[faceI]);
            rawStabBoundaryMax =
                max(rawStabBoundaryMax, mag(rawPatch[faceI]));
            rawForceBoundarySqr += magSqr(rawForcePatch[faceI]);
            rawForceBoundaryMax =
                max(rawForceBoundaryMax, mag(rawForcePatch[faceI]));
            productionBoundarySqr += magSqr(productionPatch[faceI]);
            productionBoundaryMax =
                max(productionBoundaryMax, mag(productionPatch[faceI]));

            if (fixedPatch)
            {
                productionFixedBoundarySqr += magSqr(productionPatch[faceI]);
                productionFixedBoundaryMax =
                    max
                    (
                        productionFixedBoundaryMax,
                        mag(productionPatch[faceI])
                    );
                nFixedBoundaryFaces += 1.0;
            }
            else if (tractionPatch)
            {
                productionTractionBoundarySqr +=
                    magSqr(productionPatch[faceI]);
                productionTractionBoundaryMax =
                    max
                    (
                        productionTractionBoundaryMax,
                        mag(productionPatch[faceI])
                    );
                nTractionBoundaryFaces += 1.0;
            }
            else if (otherPatch)
            {
                productionOtherBoundarySqr += magSqr(productionPatch[faceI]);
                productionOtherBoundaryMax =
                    max
                    (
                        productionOtherBoundaryMax,
                        mag(productionPatch[faceI])
                    );
                nOtherBoundaryFaces += 1.0;
            }

            patchRawSqr += magSqr(rawForcePatch[faceI]);
            patchRawMax = max(patchRawMax, mag(rawForcePatch[faceI]));
            patchProductionSqr += magSqr(productionPatch[faceI]);
            patchProductionMax =
                max(patchProductionMax, mag(productionPatch[faceI]));

            if (!mesh.boundary()[patchI].coupled())
            {
                rawBoundaryResultant += rawForcePatch[faceI];
                productionBoundaryResultant += productionPatch[faceI];
            }
            patchRawResultant += rawForcePatch[faceI];
            patchProductionResultant += productionPatch[faceI];
            nBoundaryFaces += 1.0;
            nPatchFaces += 1.0;
        }

        reduce(patchRawResultant, sumOp<vector>());
        reduce(patchProductionResultant, sumOp<vector>());
        reduce(patchRawSqr, sumOp<scalar>());
        reduce(patchRawMax, maxOp<scalar>());
        reduce(patchProductionSqr, sumOp<scalar>());
        reduce(patchProductionMax, maxOp<scalar>());
        reduce(nPatchFaces, sumOp<scalar>());

        if (nPatchFaces > 0.5)
        {
            Info<< "    " << testName << " patch "
                << mesh.boundary()[patchI].name()
                << " raw force RMS/Linf/net = "
                << sqrt(patchRawSqr/(nPatchFaces + VSMALL)) << " / "
                << patchRawMax << " / " << patchRawResultant << nl
                << "    " << testName << " patch "
                << mesh.boundary()[patchI].name()
                << " production-corrected force RMS/Linf/net = "
                << sqrt(patchProductionSqr/(nPatchFaces + VSMALL)) << " / "
                << patchProductionMax << " / "
                << patchProductionResultant << nl;
        }
    }

    vector residualIntegral(vector::zero);
    vector residualMoment(vector::zero);
    vector internalResidualIntegral(vector::zero);
    scalar residualSqr = 0.0;
    scalar residualMax = 0.0;
    scalar residualCellForceSqr = 0.0;
    scalar residualCellForceMax = 0.0;
    scalar nResidualCells = 0.0;
    scalar volume = 0.0;

    forAll(residual, cellI)
    {
        const scalar Vc = mesh.V()[cellI];
        const vector cellForce = Vc*residual[cellI];
        residualIntegral += cellForce;
        residualMoment += mesh.C()[cellI] ^ cellForce;
        internalResidualIntegral += Vc*internalResidual[cellI];
        residualSqr += Vc*magSqr(residual[cellI]);
        residualMax = max(residualMax, mag(residual[cellI]));
        residualCellForceSqr += magSqr(cellForce);
        residualCellForceMax = max(residualCellForceMax, mag(cellForce));
        nResidualCells += 1.0;
        volume += Vc;
    }

    reduce(rawStabInternalSqr, sumOp<scalar>());
    reduce(rawStabInternalMax, maxOp<scalar>());
    reduce(rawForceInternalSqr, sumOp<scalar>());
    reduce(rawForceInternalMax, maxOp<scalar>());
    reduce(productionInternalSqr, sumOp<scalar>());
    reduce(productionInternalMax, maxOp<scalar>());
    reduce(nInternalFaces, sumOp<scalar>());
    reduce(rawStabBoundarySqr, sumOp<scalar>());
    reduce(rawStabBoundaryMax, maxOp<scalar>());
    reduce(rawForceBoundarySqr, sumOp<scalar>());
    reduce(rawForceBoundaryMax, maxOp<scalar>());
    reduce(productionBoundarySqr, sumOp<scalar>());
    reduce(productionBoundaryMax, maxOp<scalar>());
    reduce(productionFixedBoundarySqr, sumOp<scalar>());
    reduce(productionFixedBoundaryMax, maxOp<scalar>());
    reduce(nFixedBoundaryFaces, sumOp<scalar>());
    reduce(productionTractionBoundarySqr, sumOp<scalar>());
    reduce(productionTractionBoundaryMax, maxOp<scalar>());
    reduce(nTractionBoundaryFaces, sumOp<scalar>());
    reduce(productionOtherBoundarySqr, sumOp<scalar>());
    reduce(productionOtherBoundaryMax, maxOp<scalar>());
    reduce(nOtherBoundaryFaces, sumOp<scalar>());
    reduce(nBoundaryFaces, sumOp<scalar>());
    reduce(rawBoundaryResultant, sumOp<vector>());
    reduce(productionBoundaryResultant, sumOp<vector>());
    reduce(residualIntegral, sumOp<vector>());
    reduce(residualMoment, sumOp<vector>());
    reduce(internalResidualIntegral, sumOp<vector>());
    reduce(residualSqr, sumOp<scalar>());
    reduce(residualMax, maxOp<scalar>());
    reduce(residualCellForceSqr, sumOp<scalar>());
    reduce(residualCellForceMax, maxOp<scalar>());
    reduce(nResidualCells, sumOp<scalar>());
    reduce(volume, sumOp<scalar>());

    const scalar rawStabInternalRms =
        sqrt(rawStabInternalSqr/(nInternalFaces + VSMALL));
    const scalar rawStabBoundaryRms =
        sqrt(rawStabBoundarySqr/(nBoundaryFaces + VSMALL));
    const scalar rawForceInternalRms =
        sqrt(rawForceInternalSqr/(nInternalFaces + VSMALL));
    const scalar rawForceBoundaryRms =
        sqrt(rawForceBoundarySqr/(nBoundaryFaces + VSMALL));
    const scalar productionInternalRms =
        sqrt(productionInternalSqr/(nInternalFaces + VSMALL));
    const scalar productionBoundaryRms =
        sqrt(productionBoundarySqr/(nBoundaryFaces + VSMALL));
    const scalar productionFixedBoundaryRms =
        sqrt(productionFixedBoundarySqr/(nFixedBoundaryFaces + VSMALL));
    const scalar productionTractionBoundaryRms =
        sqrt
        (
            productionTractionBoundarySqr
           /(nTractionBoundaryFaces + VSMALL)
        );
    const scalar productionOtherBoundaryRms =
        sqrt(productionOtherBoundarySqr/(nOtherBoundaryFaces + VSMALL));
    const scalar residualRms =
        sqrt(residualSqr/(volume + VSMALL));
    const scalar residualCellForceRms =
        sqrt(residualCellForceSqr/(nResidualCells + VSMALL));

    const scalar vanishScale =
        max
        (
            max(productionInternalMax, productionBoundaryMax),
            residualCellForceMax
        );
    const bool vanishPassed =
        !shouldVanish || vanishScale < 1e-10;

    Info<< "Momentum stabilisation consistency test '" << testName
        << "':" << nl
        << "    should vanish for exact input = " << shouldVanish << nl
        << "    raw face stabilisation internal RMS/Linf = "
        << rawStabInternalRms << " / " << rawStabInternalMax << nl
        << "    raw face stabilisation boundary RMS/Linf = "
        << rawStabBoundaryRms << " / " << rawStabBoundaryMax << nl
        << "    raw scaled force internal RMS/Linf = "
        << rawForceInternalRms << " / " << rawForceInternalMax << nl
        << "    raw scaled force boundary RMS/Linf = "
        << rawForceBoundaryRms << " / " << rawForceBoundaryMax << nl
        << "    production-corrected force internal RMS/Linf = "
        << productionInternalRms << " / " << productionInternalMax << nl
        << "    production-corrected force boundary RMS/Linf = "
        << productionBoundaryRms << " / " << productionBoundaryMax << nl
        << "    production-corrected fixed-boundary RMS/Linf = "
        << productionFixedBoundaryRms << " / "
        << productionFixedBoundaryMax << nl
        << "    production-corrected traction-boundary RMS/Linf = "
        << productionTractionBoundaryRms << " / "
        << productionTractionBoundaryMax << nl
        << "    production-corrected other-boundary RMS/Linf = "
        << productionOtherBoundaryRms << " / "
        << productionOtherBoundaryMax << nl
        << "    final production cell residual RMS/Linf = "
        << residualRms << " / " << residualMax << nl
        << "    final production extensive cell force RMS/Linf = "
        << residualCellForceRms << " / " << residualCellForceMax << nl
        << "    residual volume integral = "
        << residualIntegral << nl
        << "    residual moment about origin = "
        << residualMoment << nl
        << "    internal/coupled-face conservation error = "
        << mag(internalResidualIntegral) << nl
        << "    raw non-coupled boundary resultant = "
        << rawBoundaryResultant << nl
        << "    production-corrected non-coupled boundary resultant = "
        << productionBoundaryResultant << nl
        << "    consistency status = "
        << (vanishPassed ? "pass/report" : "FAILED expected-zero test")
        << endl;

    if (mesh.foundObject<volScalarField>("apexLayerIndex"))
    {
        const volScalarField& apexLayer =
            mesh.lookupObject<volScalarField>("apexLayerIndex");

        label minLayer = labelMax;
        label maxLayer = -labelMax;

        forAll(apexLayer, cellI)
        {
            const label layerI = label(floor(apexLayer[cellI] + 0.5));
            minLayer = min(minLayer, layerI);
            maxLayer = max(maxLayer, layerI);
        }

        reduce(minLayer, minOp<label>());
        reduce(maxLayer, maxOp<label>());

        if (minLayer <= maxLayer && (maxLayer - minLayer) < 100)
        {
            for (label layerI = minLayer; layerI <= maxLayer; ++layerI)
            {
                scalar layerVolume = 0.0;
                scalar layerSqr = 0.0;
                scalar layerMax = 0.0;
                vector layerIntegral(vector::zero);

                forAll(apexLayer, cellI)
                {
                    const label curLayer =
                        label(floor(apexLayer[cellI] + 0.5));

                    if (curLayer != layerI)
                    {
                        continue;
                    }

                    const scalar Vc = mesh.V()[cellI];
                    layerVolume += Vc;
                    layerSqr += Vc*magSqr(residual[cellI]);
                    layerMax = max(layerMax, mag(residual[cellI]));
                    layerIntegral += Vc*residual[cellI];
                }

                reduce(layerVolume, sumOp<scalar>());
                reduce(layerSqr, sumOp<scalar>());
                reduce(layerMax, maxOp<scalar>());
                reduce(layerIntegral, sumOp<vector>());

                if (layerVolume > VSMALL)
                {
                    Info<< "    apexLayerIndex " << layerI
                        << " residual RMS/Linf/integral = "
                        << sqrt(layerSqr/(layerVolume + VSMALL))
                        << " / " << layerMax
                        << " / " << layerIntegral << nl;
                }
            }
        }
    }

    if (momentumStabilisationConsistencyWriteFields_)
    {
        volVectorField prescribedDWrite
        (
            IOobject
            (
                "prescribedD_" + testName,
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            prescribedD
        );

        prescribedDWrite.write();
        rawFaceStabilisation.write();
        rawFaceContribution.write();
        productionCorrectedFaceContribution.write();
        residual.write();
    }
}


void
nonLinGeomTotalLagTotalDispSolid::runAffineKinematicsConsistencyTest() const
{
    const fvMesh& mesh = this->mesh();

    const tensor Fprescribed
    (
        1.08,  0.12, -0.04,
        0.03,  0.93,  0.05,
       -0.02,  0.07,  1.11
    );
    const tensor displacementMap(Fprescribed - I);
    const tensor exactGradD(gradDFromDisplacementMap(displacementMap));
    const scalar Jexact = det(Fprescribed);

    volVectorField affineD
    (
        IOobject
        (
            "affineKinematicsPrescribedD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero),
        diagnosticVectorPatchTypes(mesh, D())
    );

    volTensorField exactGradDField
    (
        IOobject
        (
            "grad(affineKinematicsPrescribedD)",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero),
        diagnosticTensorPatchTypes(mesh, gradD())
    );

    forAll(affineD, cellI)
    {
        affineD[cellI] = displacementMap & mesh.C()[cellI];
        exactGradDField[cellI] = exactGradD;
    }

    affineD.correctBoundaryConditions();
    exactGradDField.correctBoundaryConditions();

    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        vectorField& DPatch = affineD.boundaryFieldRef()[patchI];
        tensorField& gradPatch = exactGradDField.boundaryFieldRef()[patchI];

        forAll(DPatch, faceI)
        {
            DPatch[faceI] = displacementMap & CfPatch[faceI];
            gradPatch[faceI] = exactGradD;
        }
    }

    pointVectorField exactPointD
    (
        IOobject
        (
            "affineKinematicsExactPointD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointD()
    );

    forAll(exactPointD, pointI)
    {
        exactPointD[pointI] = displacementMap & mesh.points()[pointI];
    }

    pointVectorField prescribedPointD
    (
        IOobject
        (
            "affineKinematicsPrescribedPointD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        exactPointD
    );

    volTensorField reconstructedGradD
    (
        IOobject
        (
            "affineKinematicsReconstructedGradD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    const_cast<mechanicalModel&>(mechanical()).grad
    (
        affineD,
        prescribedPointD,
        reconstructedGradD
    );

    volTensorField gradDError
    (
        IOobject
        (
            "affineKinematicsGradDError",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        reconstructedGradD - exactGradDField
    );

    volScalarField JCellGradAffine
    (
        IOobject
        (
            "affineKinematicsJCellGrad",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(I + reconstructedGradD.T())
    );

    volScalarField JCellGradError
    (
        IOobject
        (
            "affineKinematicsJCellGradError",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JCellGradAffine
      - dimensionedScalar("Jexact", dimless, Jexact)
    );

    scalarField deformedVolumesPrescribed(mesh.nCells(), 0.0);
    scalarField deformedVolumesExact(mesh.nCells(), 0.0);
    calcDeformedCellVolumes(prescribedPointD, deformedVolumesPrescribed);
    calcDeformedCellVolumes(exactPointD, deformedVolumesExact);

    volScalarField JDeformedVolumeAffine
    (
        IOobject
        (
            "affineKinematicsJDeformedVolume",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    volScalarField JDeformedVolumeExactPointD
    (
        IOobject
        (
            "affineKinematicsJDeformedVolumeExactPointD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    forAll(JDeformedVolumeAffine, cellI)
    {
        JDeformedVolumeAffine[cellI] =
            deformedVolumesPrescribed[cellI]/mesh.V()[cellI];
        JDeformedVolumeExactPointD[cellI] =
            deformedVolumesExact[cellI]/mesh.V()[cellI];
    }

    volScalarField JDeformedVolumeError
    (
        IOobject
        (
            "affineKinematicsJDeformedVolumeError",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JDeformedVolumeAffine
      - dimensionedScalar("Jexact", dimless, Jexact)
    );

    volScalarField JDeformedVolumeExactPointDError
    (
        IOobject
        (
            "affineKinematicsJDeformedVolumeExactPointDError",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        JDeformedVolumeExactPointD
      - dimensionedScalar("Jexact", dimless, Jexact)
    );

    scalar gradSqr = 0.0;
    scalar gradMax = 0.0;
    scalar volume = 0.0;
    forAll(gradDError, cellI)
    {
        const scalar Vc = mesh.V()[cellI];
        gradSqr += Vc*magSqr(gradDError[cellI]);
        gradMax = max(gradMax, mag(gradDError[cellI]));
        volume += Vc;
    }

    reduce(gradSqr, sumOp<scalar>());
    reduce(gradMax, maxOp<scalar>());
    reduce(volume, sumOp<scalar>());

    scalar pointSqr = 0.0;
    scalar pointMax = 0.0;
    scalar nPoints = 0.0;

    forAll(prescribedPointD, pointI)
    {
        const vector diff =
            prescribedPointD[pointI] - exactPointD[pointI];

        pointSqr += magSqr(diff);
        pointMax = max(pointMax, mag(diff));
        nPoints += 1.0;
    }

    reduce(pointSqr, sumOp<scalar>());
    reduce(pointMax, maxOp<scalar>());
    reduce(nPoints, sumOp<scalar>());

    const ScalarDiagnosticStats JCellStats =
        scalarStats(JCellGradError, mesh.V());
    const ScalarDiagnosticStats JGeomStats =
        scalarStats(JDeformedVolumeError, mesh.V());
    const ScalarDiagnosticStats JGeomExactPointStats =
        scalarStats(JDeformedVolumeExactPointDError, mesh.V());

    Info<< "Affine kinematics consistency test:" << nl
        << "    Fprescribed = " << Fprescribed << nl
        << "    det(Fprescribed) = " << Jexact << nl
        << "    expected grad(D) = Fprescribed.T - I = "
        << exactGradD << nl
        << "    grad(D) reconstruction RMS/Linf error = "
        << sqrt(gradSqr/(volume + VSMALL)) << " / " << gradMax << nl
        << "    prescribed pointD RMS/Linf error = "
        << sqrt(pointSqr/(nPoints + VSMALL)) << " / " << pointMax << nl;
    printScalarStats("JCellGrad - det(Fprescribed)", JCellStats);
    printScalarStats("JDeformedVolume - det(Fprescribed)", JGeomStats);
    printScalarStats
    (
        "JDeformedVolumeExactPointD - det(Fprescribed)",
        JGeomExactPointStats
    );

    const cellShapeList& cellShapes = mesh.cellShapes();
    wordList cellTypeNames(3);
    cellTypeNames[0] = "hex";
    cellTypeNames[1] = "wedge";
    cellTypeNames[2] = "other";

    forAll(cellTypeNames, typeI)
    {
        scalar typeVolume = 0.0;
        scalar typeJCellSqr = 0.0;
        scalar typeJGeomSqr = 0.0;
        scalar typeJCellMax = 0.0;
        scalar typeJGeomMax = 0.0;

        forAll(cellShapes, cellI)
        {
            const word cellModelName = cellShapes[cellI].model().name();
            const bool inType =
                (
                    cellTypeNames[typeI] == "other"
                  ? (cellModelName != "hex" && cellModelName != "wedge")
                  : (cellModelName == cellTypeNames[typeI])
                );

            if (!inType)
            {
                continue;
            }

            const scalar Vc = mesh.V()[cellI];
            typeVolume += Vc;
            typeJCellSqr += Vc*sqr(JCellGradError[cellI]);
            typeJGeomSqr += Vc*sqr(JDeformedVolumeError[cellI]);
            typeJCellMax = max(typeJCellMax, mag(JCellGradError[cellI]));
            typeJGeomMax = max(typeJGeomMax, mag(JDeformedVolumeError[cellI]));
        }

        reduce(typeVolume, sumOp<scalar>());
        reduce(typeJCellSqr, sumOp<scalar>());
        reduce(typeJGeomSqr, sumOp<scalar>());
        reduce(typeJCellMax, maxOp<scalar>());
        reduce(typeJGeomMax, maxOp<scalar>());

        if (typeVolume > VSMALL)
        {
            Info<< "    cell type " << cellTypeNames[typeI]
                << " JCell/JDeformed RMS errors = "
                << sqrt(typeJCellSqr/(typeVolume + VSMALL))
                << " / "
                << sqrt(typeJGeomSqr/(typeVolume + VSMALL))
                << ", Linf = "
                << typeJCellMax << " / " << typeJGeomMax << nl;
        }
    }

    boolList boundaryAdjacent(mesh.nCells(), false);
    forAll(mesh.boundary(), patchI)
    {
        const labelUList& faceCells = mesh.boundary()[patchI].faceCells();

        forAll(faceCells, faceI)
        {
            boundaryAdjacent[faceCells[faceI]] = true;
        }
    }

    for (label adjacentI = 0; adjacentI < 2; ++adjacentI)
    {
        scalar regionVolume = 0.0;
        scalar regionJCellSqr = 0.0;
        scalar regionJGeomSqr = 0.0;
        scalar regionJCellMax = 0.0;
        scalar regionJGeomMax = 0.0;

        const bool wantAdjacent = adjacentI == 1;

        forAll(boundaryAdjacent, cellI)
        {
            if (boundaryAdjacent[cellI] != wantAdjacent)
            {
                continue;
            }

            const scalar Vc = mesh.V()[cellI];
            regionVolume += Vc;
            regionJCellSqr += Vc*sqr(JCellGradError[cellI]);
            regionJGeomSqr += Vc*sqr(JDeformedVolumeError[cellI]);
            regionJCellMax = max(regionJCellMax, mag(JCellGradError[cellI]));
            regionJGeomMax = max(regionJGeomMax, mag(JDeformedVolumeError[cellI]));
        }

        reduce(regionVolume, sumOp<scalar>());
        reduce(regionJCellSqr, sumOp<scalar>());
        reduce(regionJGeomSqr, sumOp<scalar>());
        reduce(regionJCellMax, maxOp<scalar>());
        reduce(regionJGeomMax, maxOp<scalar>());

        if (regionVolume > VSMALL)
        {
            Info<< "    "
                << (wantAdjacent ? "boundary-adjacent" : "interior")
                << " cells JCell/JDeformed RMS errors = "
                << sqrt(regionJCellSqr/(regionVolume + VSMALL))
                << " / "
                << sqrt(regionJGeomSqr/(regionVolume + VSMALL))
                << ", Linf = "
                << regionJCellMax << " / " << regionJGeomMax << nl;
        }
    }

    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            continue;
        }

        const vectorField& CfPatch = mesh.Cf().boundaryField()[patchI];
        const vectorField& DPatch = affineD.boundaryField()[patchI];
        const tensorField& gradPatch = exactGradDField.boundaryField()[patchI];
        const labelUList& faceCells = mesh.boundary()[patchI].faceCells();

        scalar patchDSqr = 0.0;
        scalar patchDMax = 0.0;
        scalar patchGradSqr = 0.0;
        scalar patchGradMax = 0.0;
        scalar patchJCellSqr = 0.0;
        scalar patchJCellMax = 0.0;
        scalar patchJGeomSqr = 0.0;
        scalar patchJGeomMax = 0.0;
        scalar nPatchFaces = 0.0;

        forAll(faceCells, faceI)
        {
            const vector DError =
                DPatch[faceI] - (displacementMap & CfPatch[faceI]);
            const tensor gradError = gradPatch[faceI] - exactGradD;
            const label cellI = faceCells[faceI];

            patchDSqr += magSqr(DError);
            patchDMax = max(patchDMax, mag(DError));
            patchGradSqr += magSqr(gradError);
            patchGradMax = max(patchGradMax, mag(gradError));
            patchJCellSqr += sqr(JCellGradError[cellI]);
            patchJCellMax =
                max(patchJCellMax, mag(JCellGradError[cellI]));
            patchJGeomSqr += sqr(JDeformedVolumeError[cellI]);
            patchJGeomMax =
                max(patchJGeomMax, mag(JDeformedVolumeError[cellI]));
            nPatchFaces += 1.0;
        }

        reduce(patchDSqr, sumOp<scalar>());
        reduce(patchDMax, maxOp<scalar>());
        reduce(patchGradSqr, sumOp<scalar>());
        reduce(patchGradMax, maxOp<scalar>());
        reduce(patchJCellSqr, sumOp<scalar>());
        reduce(patchJCellMax, maxOp<scalar>());
        reduce(patchJGeomSqr, sumOp<scalar>());
        reduce(patchJGeomMax, maxOp<scalar>());
        reduce(nPatchFaces, sumOp<scalar>());

        if (nPatchFaces > 0.5)
        {
            Info<< "    boundary patch " << mesh.boundary()[patchI].name()
                << " D/grad(D) RMS errors = "
                << sqrt(patchDSqr/(nPatchFaces + VSMALL)) << " / "
                << sqrt(patchGradSqr/(nPatchFaces + VSMALL))
                << ", Linf = " << patchDMax << " / " << patchGradMax << nl
                << "    boundary patch " << mesh.boundary()[patchI].name()
                << " owner-cell JCell/JDeformed RMS errors = "
                << sqrt(patchJCellSqr/(nPatchFaces + VSMALL))
                << " / "
                << sqrt(patchJGeomSqr/(nPatchFaces + VSMALL))
                << ", Linf = "
                << patchJCellMax << " / " << patchJGeomMax << nl;
        }
    }

    if (mesh.foundObject<volScalarField>("apexLayerIndex"))
    {
        const volScalarField& apexLayer =
            mesh.lookupObject<volScalarField>("apexLayerIndex");

        label minLayer = labelMax;
        label maxLayer = -labelMax;

        forAll(apexLayer, cellI)
        {
            const label layerI = label(floor(apexLayer[cellI] + 0.5));
            minLayer = min(minLayer, layerI);
            maxLayer = max(maxLayer, layerI);
        }

        reduce(minLayer, minOp<label>());
        reduce(maxLayer, maxOp<label>());

        if (minLayer <= maxLayer && (maxLayer - minLayer) < 100)
        {
            for (label layerI = minLayer; layerI <= maxLayer; ++layerI)
            {
                scalar layerVolume = 0.0;
                scalar layerJCellSqr = 0.0;
                scalar layerJGeomSqr = 0.0;
                scalar layerJCellMax = 0.0;
                scalar layerJGeomMax = 0.0;

                forAll(apexLayer, cellI)
                {
                    const label curLayer =
                        label(floor(apexLayer[cellI] + 0.5));

                    if (curLayer != layerI)
                    {
                        continue;
                    }

                    const scalar Vc = mesh.V()[cellI];
                    layerVolume += Vc;
                    layerJCellSqr += Vc*sqr(JCellGradError[cellI]);
                    layerJGeomSqr += Vc*sqr(JDeformedVolumeError[cellI]);
                    layerJCellMax =
                        max(layerJCellMax, mag(JCellGradError[cellI]));
                    layerJGeomMax =
                        max(layerJGeomMax, mag(JDeformedVolumeError[cellI]));
                }

                reduce(layerVolume, sumOp<scalar>());
                reduce(layerJCellSqr, sumOp<scalar>());
                reduce(layerJGeomSqr, sumOp<scalar>());
                reduce(layerJCellMax, maxOp<scalar>());
                reduce(layerJGeomMax, maxOp<scalar>());

                if (layerVolume > VSMALL)
                {
                    Info<< "    apexLayerIndex " << layerI
                        << " JCell/JDeformed RMS errors = "
                        << sqrt(layerJCellSqr/(layerVolume + VSMALL))
                        << " / "
                        << sqrt(layerJGeomSqr/(layerVolume + VSMALL))
                        << ", Linf = "
                        << layerJCellMax << " / " << layerJGeomMax << nl;
                }
            }
        }
    }

    if (affineKinematicsConsistencyWriteFields_)
    {
        affineD.write();
        exactGradDField.write();
        prescribedPointD.write();
        exactPointD.write();
        reconstructedGradD.write();
        gradDError.write();
        JCellGradAffine.write();
        JCellGradError.write();
        JDeformedVolumeAffine.write();
        JDeformedVolumeExactPointD.write();
        JDeformedVolumeError.write();
        JDeformedVolumeExactPointDError.write();
    }
}


void nonLinGeomTotalLagTotalDispSolid::reportScalarExtremumLocation
(
    const word& name,
    const volScalarField& field,
    const scalar extremumValue
) const
{
    const globalIndex globalCells(mesh().nCells());
    const scalar tolerance =
        10*SMALL*max(scalar(1), mag(extremumValue));

    forAll(field, cellI)
    {
        if (mag(field[cellI] - extremumValue) <= tolerance)
        {
            Pout<< "    " << name << " location:" << nl
                << "        processor = " << Pstream::myProcNo() << nl
                << "        local cell ID = " << cellI << nl
                << "        global cell ID = "
                << globalCells.toGlobal(cellI) << nl
                << "        value = " << field[cellI] << nl
                << "        reference cell centre = " << mesh().C()[cellI]
                << nl
                << "        reference cell volume = " << mesh().V()[cellI]
                << endl;
        }
    }
}


void nonLinGeomTotalLagTotalDispSolid::reportIncompressibilityDiagnostics() const
{
    scalar minJ = VGREAT;
    scalar maxJ = -VGREAT;
    scalar volume = 0.0;
    scalar volumeJ = 0.0;
    scalar volumeJminus1Sqr = 0.0;
    scalar maxMagJminus1 = 0.0;

    const scalarField& V = mesh().V();

    forAll(J_, cellI)
    {
        const scalar Jc = J_[cellI];
        const scalar Jminus1 = Jc - 1.0;
        const scalar Vc = V[cellI];

        minJ = min(minJ, Jc);
        maxJ = max(maxJ, Jc);
        volume += Vc;
        volumeJ += Vc*Jc;
        volumeJminus1Sqr += Vc*sqr(Jminus1);
        maxMagJminus1 = max(maxMagJminus1, mag(Jminus1));
    }

    reduce(minJ, minOp<scalar>());
    reduce(maxJ, maxOp<scalar>());
    reduce(volume, sumOp<scalar>());
    reduce(volumeJ, sumOp<scalar>());
    reduce(volumeJminus1Sqr, sumOp<scalar>());
    reduce(maxMagJminus1, maxOp<scalar>());

    Info<< "Incompressibility diagnostics:" << nl
        << "    min(J)                  = " << minJ << nl
        << "    max(J)                  = " << maxJ << nl
        << "    volume-weighted mean(J) = "
        << volumeJ/(volume + VSMALL) << nl
        << "    L2(J-1)                 = "
        << sqrt(volumeJminus1Sqr/(volume + VSMALL)) << nl
        << "    Linf(J-1)               = "
        << maxMagJminus1 << endl;

    if
    (
        !writeGeometricJacobianDiagnostics_
     && !reportExtendedIncompressibilityDiagnostics_
    )
    {
        return;
    }

    scalar volumeJgeom = 0.0;
    scalar volumeJgeomMinus1Sqr = 0.0;
    scalar maxMagJgeomMinus1 = 0.0;
    scalar jgradMinusJgeomMean = 0.0;
    scalar jgradMinusJgeomSqr = 0.0;
    scalar jgradMinusJgeomLinf = 0.0;
    scalar corrNumerator = 0.0;
    scalar jgradMinus1Sqr = 0.0;
    scalar jgeomMinus1Sqr = 0.0;

    scalar minJgeom = VGREAT;
    scalar maxJgeom = -VGREAT;

    forAll(Jgeom_, cellI)
    {
        const scalar Vc = V[cellI];
        const scalar jgradMinus1 = J_[cellI] - 1.0;
        const scalar jgeom = Jgeom_[cellI];
        const scalar jgeomMinus1 = jgeom - 1.0;
        const scalar diff = JgradMinusJgeom_[cellI];

        minJgeom = min(minJgeom, jgeom);
        maxJgeom = max(maxJgeom, jgeom);
        volumeJgeom += Vc*jgeom;
        volumeJgeomMinus1Sqr += Vc*sqr(jgeomMinus1);
        maxMagJgeomMinus1 = max(maxMagJgeomMinus1, mag(jgeomMinus1));

        jgradMinusJgeomMean += Vc*diff;
        jgradMinusJgeomSqr += Vc*sqr(diff);
        jgradMinusJgeomLinf = max(jgradMinusJgeomLinf, mag(diff));

        corrNumerator += Vc*jgradMinus1*jgeomMinus1;
        jgradMinus1Sqr += Vc*sqr(jgradMinus1);
        jgeomMinus1Sqr += Vc*sqr(jgeomMinus1);
    }

    reduce(minJgeom, minOp<scalar>());
    reduce(maxJgeom, maxOp<scalar>());
    reduce(volumeJgeom, sumOp<scalar>());
    reduce(volumeJgeomMinus1Sqr, sumOp<scalar>());
    reduce(maxMagJgeomMinus1, maxOp<scalar>());
    reduce(jgradMinusJgeomMean, sumOp<scalar>());
    reduce(jgradMinusJgeomSqr, sumOp<scalar>());
    reduce(jgradMinusJgeomLinf, maxOp<scalar>());
    reduce(corrNumerator, sumOp<scalar>());
    reduce(jgradMinus1Sqr, sumOp<scalar>());
    reduce(jgeomMinus1Sqr, sumOp<scalar>());

    const scalar corr =
        corrNumerator/(sqrt(jgradMinus1Sqr*jgeomMinus1Sqr) + VSMALL);

    Info<< "Geometric incompressibility diagnostics:" << nl
        << "    Jgrad = det(I + grad(D)^T)" << nl
        << "        min                 = " << minJ << nl
        << "        max                 = " << maxJ << nl
        << "        volume-weighted mean= "
        << volumeJ/(volume + VSMALL) << nl
        << "        L2(Jgrad-1)         = "
        << sqrt(volumeJminus1Sqr/(volume + VSMALL)) << nl
        << "        Linf(Jgrad-1)       = "
        << maxMagJminus1 << nl
        << "    Jgeom = Vdeformed(pointD)/Vreference" << nl
        << "        min                 = " << minJgeom << nl
        << "        max                 = " << maxJgeom << nl
        << "        volume-weighted mean= "
        << volumeJgeom/(volume + VSMALL) << nl
        << "        L2(Jgeom-1)         = "
        << sqrt(volumeJgeomMinus1Sqr/(volume + VSMALL)) << nl
        << "        Linf(Jgeom-1)       = "
        << maxMagJgeomMinus1 << nl
        << "    JgradMinusJgeom:" << nl
        << "        volume-weighted mean= "
        << jgradMinusJgeomMean/(volume + VSMALL) << nl
        << "        volume-weighted RMS = "
        << sqrt(jgradMinusJgeomSqr/(volume + VSMALL)) << nl
        << "        Linf                = "
        << jgradMinusJgeomLinf << nl
        << "        corr(Jgrad-1,Jgeom-1)= " << corr << endl;

    reportScalarExtremumLocation("min(Jgrad)", J_, minJ);
    reportScalarExtremumLocation("max(Jgrad)", J_, maxJ);
    reportScalarExtremumLocation("min(Jgeom)", Jgeom_, minJgeom);
    reportScalarExtremumLocation("max(Jgeom)", Jgeom_, maxJgeom);
}


#ifdef USE_PETSC

void nonLinGeomTotalLagTotalDispSolid::reportSnesLineSearchTrial
(
    const Vec x,
    const Vec residual,
    const bool forceReport
)
{
    if (!reportSnesTrialDiagnostics_ && !forceReport)
    {
        return;
    }

    SNESLineSearch lineSearch = nullptr;
    AssertPETSc(SNESGetLineSearch(snes(), &lineSearch));

    Vec acceptedX = nullptr;
    Vec step = nullptr;
    Vec candidateX = nullptr;
    AssertPETSc
    (
        SNESLineSearchGetVecs
        (
            lineSearch,
            &acceptedX,
            nullptr,
            &step,
            &candidateX,
            nullptr
        )
    );

    // Matrix-free KSP perturbations also call formResidual. Only the PETSc
    // line-search work vector is a nonlinear step trial.
    if (x != candidateX && !forceReport)
    {
        return;
    }

    SNES snesObject = nullptr;
    AssertPETSc(SNESLineSearchGetSNES(lineSearch, &snesObject));

    PetscInt snesIteration = -1;
    AssertPETSc(SNESGetIterationNumber(snesObject, &snesIteration));

    const bool firstReportedTrial =
        lastSnesTrialDiagnosticIteration_ < 0;
    if (lastSnesTrialDiagnosticIteration_ != snesIteration)
    {
        lastSnesTrialDiagnosticIteration_ = snesIteration;
        snesLineSearchTrialCounter_ = 0;
    }
    ++snesLineSearchTrialCounter_;

    PetscInt blockSize = 0;
    PetscInt localSize = 0;
    AssertPETSc(VecGetBlockSize(candidateX, &blockSize));
    AssertPETSc(VecGetLocalSize(candidateX, &localSize));

    const PetscScalar* acceptedValues = nullptr;
    const PetscScalar* stepValues = nullptr;
    const PetscScalar* candidateValues = nullptr;
    AssertPETSc(VecGetArrayRead(acceptedX, &acceptedValues));
    AssertPETSc(VecGetArrayRead(step, &stepValues));
    AssertPETSc(VecGetArrayRead(candidateX, &candidateValues));

    scalar stepNormSqr = 0;
    scalar attemptedIncrementNormSqr = 0;
    scalar maxDisplacementIncrement = 0;
    for (PetscInt i = 0; i < localSize; ++i)
    {
        stepNormSqr += sqr(PetscRealPart(stepValues[i]));
        attemptedIncrementNormSqr +=
            sqr
            (
                PetscRealPart(candidateValues[i])
              - PetscRealPart(acceptedValues[i])
            );
    }

    for (PetscInt i = 0; i < localSize; i += blockSize)
    {
        scalar displacementIncrementSqr = 0;
        for (PetscInt cmpt = 0; cmpt < min(PetscInt(3), blockSize); ++cmpt)
        {
            displacementIncrementSqr +=
                sqr
                (
                    PetscRealPart(candidateValues[i + cmpt])
                  - PetscRealPart(acceptedValues[i + cmpt])
                );
        }
        maxDisplacementIncrement =
            max(maxDisplacementIncrement, sqrt(displacementIncrementSqr));
    }

    AssertPETSc(VecRestoreArrayRead(candidateX, &candidateValues));
    AssertPETSc(VecRestoreArrayRead(step, &stepValues));
    AssertPETSc(VecRestoreArrayRead(acceptedX, &acceptedValues));

    reduce(stepNormSqr, sumOp<scalar>());
    reduce(attemptedIncrementNormSqr, sumOp<scalar>());
    reduce(maxDisplacementIncrement, maxOp<scalar>());

    const scalar attemptedStepLength =
        sqrt(attemptedIncrementNormSqr)/(sqrt(stepNormSqr) + VSMALL);

    scalar momentumResidualNormSqr = 0;
    scalar pressureResidualNormSqr = 0;
    if (residual)
    {
        const PetscScalar* residualValues = nullptr;
        AssertPETSc(VecGetArrayRead(residual, &residualValues));
        for (PetscInt i = 0; i < localSize; i += blockSize)
        {
            for
            (
                PetscInt cmpt = 0;
                cmpt < min(PetscInt(3), blockSize);
                ++cmpt
            )
            {
                momentumResidualNormSqr +=
                    sqr(PetscRealPart(residualValues[i + cmpt]));
            }
            if (blockSize > 3)
            {
                pressureResidualNormSqr +=
                    sqr(PetscRealPart(residualValues[i + 3]));
            }
        }
        AssertPETSc(VecRestoreArrayRead(residual, &residualValues));
    }
    else
    {
        momentumResidualNormSqr =
            std::numeric_limits<scalar>::quiet_NaN();
        pressureResidualNormSqr =
            std::numeric_limits<scalar>::quiet_NaN();
    }

    if (residual)
    {
        reduce(momentumResidualNormSqr, sumOp<scalar>());
        reduce(pressureResidualNormSqr, sumOp<scalar>());
    }

    scalar minJ = VGREAT;
    label minJCell = -1;
    forAll(J_, cellI)
    {
        if (!std::isfinite(J_[cellI]) || J_[cellI] < minJ)
        {
            minJ = J_[cellI];
            minJCell = cellI;
            if (!std::isfinite(minJ))
            {
                break;
            }
        }
    }

    if (firstReportedTrial)
    {
        Info<< "S4F_SNES_TRIAL_HEADER,"
            << "snesIteration,lineSearchTrial,attemptedStepLength,"
            << "minJ,minJCell,minJReferenceX,minJReferenceY,"
            << "minJReferenceZ,maxDisplacementIncrement,"
            << "momentumResidualNorm,pressureResidualNorm,"
            << "kspReason,kspIterations,admissible,forcedReport" << endl;
    }

    KSP ksp = nullptr;
    AssertPETSc(SNESGetKSP(snesObject, &ksp));
    KSPConvergedReason kspReason = KSP_CONVERGED_ITERATING;
    PetscInt kspIterations = 0;
    AssertPETSc(KSPGetConvergedReason(ksp, &kspReason));
    AssertPETSc(KSPGetIterationNumber(ksp, &kspIterations));

    const vector minJLocation =
        minJCell >= 0
      ? mesh().C()[minJCell]
      : vector(GREAT, GREAT, GREAT);
    Info<< "S4F_SNES_TRIAL,"
        << snesIteration << ','
        << snesLineSearchTrialCounter_ << ','
        << attemptedStepLength << ','
        << minJ << ','
        << minJCell << ','
        << minJLocation.x() << ','
        << minJLocation.y() << ','
        << minJLocation.z() << ','
        << maxDisplacementIncrement << ','
        << sqrt(momentumResidualNormSqr) << ','
        << sqrt(pressureResidualNormSqr) << ','
        << label(kspReason) << ','
        << kspIterations << ','
        << Switch(residual != nullptr) << ','
        << Switch(forceReport)
        << endl;
}


bool nonLinGeomTotalLagTotalDispSolid::unpackSolution
(
    const Vec x,
    const bool allowDomainError
)
{
    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions on D
    D.correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Enforce the boundary conditions again for any conditions that use gradD
    D.correctBoundaryConditions();

    // Increment of displacement
    DD() = D - D.oldTime();

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Total deformation gradient
    F_ = I + gradD().T();

    // Jacobian of the deformation gradient
    J_ = det(F_);
    updateJacobianDiagnostics();

    bool allJFinite = true;
    forAll(J_, cellI)
    {
        allJFinite = allJFinite && std::isfinite(J_[cellI]);
    }
    reduce(allJFinite, andOp<bool>());
    const scalar minJ = gMin(J_.primitiveField());
    const bool invalidJ = !allJFinite || minJ <= 0;

    bool recoverableLineSearchTrial = false;
    if (invalidJ && petscDomainSafeTrials_ && allowDomainError)
    {
        SNESLineSearch lineSearch = nullptr;
        AssertPETSc(SNESGetLineSearch(snes(), &lineSearch));
        Vec candidateX = nullptr;
        AssertPETSc
        (
            SNESLineSearchGetVecs
            (
                lineSearch,
                nullptr,
                nullptr,
                nullptr,
                &candidateX,
                nullptr
            )
        );
        recoverableLineSearchTrial = x == candidateX;
    }

    if (invalidJ)
    {
        if (reportSnesTrialDiagnostics_ || recoverableLineSearchTrial)
        {
            reportSnesLineSearchTrial(x, nullptr, true);
        }

        if (recoverableLineSearchTrial)
        {
            Info<< "S4F_SNES_FUNCTION_DOMAIN_ERROR,minJ=" << minJ
                << ",trialRejected=true" << endl;
            return false;
        }

        if (petscDomainSafeTrials_)
        {
            FatalErrorInFunction
                << "Invalid accepted PETSc state: minimum deformation "
                << "Jacobian J = " << minJ
                << "; J must be finite and positive"
                << exit(FatalError);
        }
    }

    // Inverse is only evaluated after a domain-safe trial has passed the
    // determinant admissibility check.
    Finv_ = inv(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Copy the scaled pressure unknown pHat from x into the
        // physical p field via p = pressureUnknownScale_ * pHat
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        scalarField pHat(pI.size(), 0.0);
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pHat, blockSize_ - 1
        );
        pI = pressureUnknownScale_*pHat;

        // Enforce the boundary conditions on p
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        applyMixedPressureStressSplit();
    }

    return true;
}


void nonLinGeomTotalLagTotalDispSolid::packSolution(Vec x)
{
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        primitiveField(D()),
        x,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    if (solvePressure())
    {
        // Insert the scaled pressure unknown pHat = p/pressureUnknownScale_
        scalarField pHat(primitiveField(p()));
        pHat /= pressureUnknownScale_;
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pHat,
            x,
            blockSize_ - 1
        );
    }
}

#endif // USE_PETSC


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagTotalDispSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    // {
    //     // Not yet implmented, although coupledUnsLinGeomLinearElasticSolid
    //     // could be combined with PETSc to achieve this.. todo!
    //     return evolveImplicitCoupled();
    // }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    // else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    // {
    //     return evolveExplicit();
    // }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
               ]
            // << solidModel::solutionAlgorithmNames_
            //    [
            //        solidModel::solutionAlgorithm::EXPLICIT
            //    ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


#ifdef USE_PETSC

label nonLinGeomTotalLagTotalDispSolid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


namespace
{

void addDirectNominalTangentToPETSc
(
    Mat jac,
    const List<tensor>& nominalTangent,
    const fvMesh& mesh,
    const label blockSize,
    const bool twoD,
    const globalIndex& globalCells,
    const bool transposeGradientIndex,
    const word& description
)
{
    if (nominalTangent.size() != 9*mesh.nFaces())
    {
        FatalErrorInFunction
            << description << " tangent has " << nominalTangent.size()
            << " entries; expected " << 9*mesh.nFaces()
            << exit(FatalError);
    }

    const vectorField& Sf = mesh.Sf();
    const scalarField& deltaCoeffs = mesh.deltaCoeffs();
    const label nDComponents = twoD ? 2 : 3;
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> tangentValues(nCoeffCmpts, 0.0);
    scalar minFaceNorm = GREAT;
    scalar maxFaceNorm = 0.0;
    scalar sumFaceNorm = 0.0;
    label nFaces = 0;

    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        const scalar magSf = mag(Sf[faceI]);
        if (magSf <= VSMALL)
        {
            continue;
        }

        const vector n = Sf[faceI]/magSf;
        const scalar deltaCoeff = deltaCoeffs[faceI];
        tensor K(tensor::zero);

        for (label rowI = 0; rowI < nDComponents; ++rowI)
        {
            for (label colI = 0; colI < nDComponents; ++colI)
            {
                scalar coefficient = 0.0;
                for (label gradDirectionI = 0;
                     gradDirectionI < 3; ++gradDirectionI)
                {
                    const label tangentComponent =
                        transposeGradientIndex
                      ? 3*gradDirectionI + colI
                      : 3*colI + gradDirectionI;
                    const tensor& dPdF =
                        nominalTangent[9*faceI + tangentComponent];
                    coefficient +=
                        (dPdF & Sf[faceI])[rowI]
                       *n[gradDirectionI]*deltaCoeff;
                }
                K(rowI, colI) = coefficient;
            }
        }

        const scalar faceNorm = mag(K);
        minFaceNorm = min(minFaceNorm, faceNorm);
        maxFaceNorm = max(maxFaceNorm, faceNorm);
        sumFaceNorm += faceNorm;
        ++nFaces;

        const label globalBlocks[2] =
        {
            globalCells.toGlobal(mesh.owner()[faceI]),
            globalCells.toGlobal
            (
                mesh.neighbour()[faceI]
            )
        };

        for (label rowBlockI = 0; rowBlockI < 2; ++rowBlockI)
        {
            for (label colBlockI = 0; colBlockI < 2; ++colBlockI)
            {
                forAll(tangentValues, valueI)
                {
                    tangentValues[valueI] = 0.0;
                }

                const scalar blockSign =
                    rowBlockI == colBlockI ? -1.0 : 1.0;
                for (label rowI = 0; rowI < nDComponents; ++rowI)
                {
                    for (label colI = 0; colI < nDComponents; ++colI)
                    {
                        tangentValues[rowI*blockSize + colI] =
                            blockSign*K(rowI, colI);
                    }
                }

                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac,
                        1,
                        &globalBlocks[rowBlockI],
                        1,
                        &globalBlocks[colBlockI],
                        tangentValues.cdata(),
                        ADD_VALUES
                    )
                );
            }
        }
    }

    reduce(minFaceNorm, minOp<scalar>());
    reduce(maxFaceNorm, maxOp<scalar>());
    reduce(sumFaceNorm, sumOp<scalar>());
    reduce(nFaces, sumOp<label>());
    Info<< description << ": faces = " << nFaces
        << ", norm range = " << minFaceNorm << " to " << maxFaceNorm
        << ", mean = "
        << sumFaceNorm/(nFaces > 0 ? scalar(nFaces) : 1.0)
        << ", internal owner-neighbour stencil only" << endl;
}

} // End anonymous namespace


label nonLinGeomTotalLagTotalDispSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label nonLinGeomTotalLagTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Copy x into the D field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    if (!unpackSolution(x, true))
    {
        return foamPetscSnesHelper::recoverableFunctionDomainError;
    }

    // Take a non-const reference to D for local use below
    volVectorField& D = const_cast<volVectorField&>(this->D());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above, its BCs corrected and the pressure component of sigma
        // updated; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        // Calculate the pressure gradient
        const volVectorField gradp(fvc::grad(p));

        // Re-calculate the pressure stabilisation parameter
        pressureStabilisation().updateScalar(p, &gradp);

        // Refresh rAUf (the positive face-interpolated reciprocal of
        // the approximate momentum equation diagonal -- the solid
        // analogue of rAUf in pressure-velocity coupling, units [m^2/Pa])
        // only when the mesh or deltaT have changed. The diagonal is
        // value-independent of D and p, so this is safe under PETSc
        // matrix-free finite-difference perturbations
        updateRAUfIfStale();

        // Calculate the pressure equation residual using the runtime model's
        // finite-strain constraint g(J). The generic model returns
        // 0.5*(J^2-1)/J; isolated derived models may override it.
        const tmp<volScalarField> tConstraint =
            mixedVolumetricConstraint(J_);
        scalarField pressureResidual
        (
          - p*rKappa()
          + pressureStabilisation().cellScalar(&rAUf(), true)
          - tConstraint()
        );

        // Make residual extensive
        pressureResidual *= mesh.V();

        // Apply the physical row scaling. pressureEqnScale_ already
        // bakes in both the user-facing pressureScaleFactor and the
        // 2*mu physical scale.
        if (pressureEqnScale_ != 1.0)
        {
            pressureResidual *= pressureEqnScale_;
        }

        if (pressureRowScaling_ == "volumeRmsForce")
        {
            scalar Vtot = gSum(mesh.V());
            const scalar L0 = cbrt(Vtot);
            const scalar L0Sqr = sqr(L0);
            forAll(pressureResidual, cellI)
            {
                // The legacy row is alpha*V_c*r_c.  Multiplication by
                // L0^2/sqrt(Vtot*V_c) gives
                // alpha*L0^2*sqrt(V_c/Vtot)*r_c, hence its PETSc 2-norm
                // is exactly alpha*L0^2 times the physical volume RMS.
                pressureResidual[cellI] *=
                    L0Sqr/sqrt(Vtot*mesh.V()[cellI]);
            }
        }

        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
        );
    }

    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    const tmp<volVectorField> tMomentumSurfaceForceDensity
    (
        momentumSurfaceForceDensity(D, nCurrent, magSfCurrent)
    );

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual(tMomentumSurfaceForceDensity());

    if (useBodyForceField_)
    {
        residual += bodyForceDensity_.internalField();
    }
    else
    {
        residual += rho()*g();
    }

    residual -= rho()
       *(
            fvc::d2dt2(D) + dampingCoeff()*fvc::ddt(D)
        );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    //residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    reportSnesLineSearchTrial(x, f);

    return 0;
}


label nonLinGeomTotalLagTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the D field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    unpackSolution(x, false);

    // Take a non-const reference to D for local use below
    volVectorField& D = const_cast<volVectorField&>(this->D());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above and its BCs corrected; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        {
            // Refresh rAUf only when the mesh or deltaT have changed
            // (the diagonal of the approximate momentum Jacobian is
            // independent of D and p values)
            updateRAUfIfStale();

            scalar compactPressureStabScale = 1.0;
            if (preconditionerLeastSquaresPressureStabilisation_)
            {
                const scalar jacobianScale =
                    pressureStabilisation().scaleFactorJacobian();
                if (mag(jacobianScale) <= VSMALL)
                {
                    FatalErrorInFunction
                        << "The exact least-squares pressure stabilisation "
                        << "Jacobian requires a non-zero "
                        << "scaleFactorJacobian so its cached compact "
                        << "laplacian can be rescaled to the production "
                        << "scaleFactor" << abort(FatalError);
                }
                compactPressureStabScale =
                    pressureStabilisation().scaleFactor()/jacobianScale;
            }

            fvScalarMatrix approxPressureJ
            (
              - pressureEqnScale_*pressureUnknownScale_*fvm::Sp(rKappa(), p)
              + pressureEqnScale_*pressureUnknownScale_
               *compactPressureStabScale
               *pressureStabilisation().scalarJacobian(p, &rAUf())
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );

            if (preconditionerLeastSquaresPressureStabilisation_)
            {
                const fvMesh& mesh = this->mesh();

                const word gradPScheme
                (
                    mesh.gradScheme("grad(" + p.name() + ')')
                );

                if
                (
                    gradPScheme != "leastSquaresS4f"
                 && gradPScheme != "leastSquaresS4fDirichlet"
                )
                {
                    FatalErrorInFunction
                        << "preconditionerLeastSquaresPressureStabilisation "
                        << "requires grad(" << p.name()
                        << ") leastSquaresS4f or leastSquaresS4fDirichlet; "
                        << "got " << gradPScheme
                        << exit(FatalError);
                }

                boolList useBoundaryFaceValues(mesh.boundary().size(), false);
                if (gradPScheme == "leastSquaresS4fDirichlet")
                {
                    forAll(useBoundaryFaceValues, patchI)
                    {
                        useBoundaryFaceValues[patchI] =
                            p.boundaryField()[patchI].fixesValue();
                    }
                }

#ifdef OPENFOAM_COM
                word lsCacheName("leastSquaresVectors");
                if (gradPScheme == "leastSquaresS4fDirichlet")
                {
                    lsCacheName += "Dirichlet";
                }
                lsCacheName += p.name();
                const leastSquaresS4fVectors& lsv =
                    leastSquaresS4fVectors::New
                    (
                        lsCacheName, mesh, useBoundaryFaceValues
                    );
#else
                const leastSquaresS4fVectors& lsv =
                    leastSquaresS4fVectors::New
                    (
                        mesh, useBoundaryFaceValues
                    );
#endif

                const surfaceVectorField& ownLs = lsv.pVectors();
                const surfaceVectorField& neiLs = lsv.nVectors();
                const surfaceScalarField& weights = mesh.weights();
                const surfaceScalarField& magSf = mesh.magSf();
                const surfaceScalarField& gamma = rAUf();
                const surfaceVectorField& nonOrthCorrectionVectors =
                    mesh.nonOrthCorrectionVectors();
                const labelUList& owner = mesh.owner();
                const labelUList& neighbour = mesh.neighbour();

                // The compact fvm::laplacian block above represents the
                // direct snGrad part of the production stabilisation.  The
                // remaining production dependency is the corrected-snGrad
                // gradient contribution minus the Rhie-Chow interpolated
                // gradient contribution.  Its coefficient is therefore
                // (nonOrthCorrectionVectors - n), not just -n.
                //
                // This option assembles the exact derivative of the
                // production stabilisation, so both its compact and wide
                // pieces use scaleFactor. scaleFactorJacobian remains the
                // legacy approximate-preconditioner coefficient only.
                const scalar pressureStabScale =
                    pressureEqnScale_
                   *pressureUnknownScale_
                   *pressureStabilisation().scaleFactor();

                // Store the complete local least-squares gradient stencil for
                // every cell.  Local columns and global off-rank columns are
                // kept separately so local PETSc insertion can use the
                // globalIndex map without confusing local labels with global
                // pressure unknown IDs.
                List<DynamicList<label>> gradientLocalColumns(mesh.nCells());
                List<DynamicList<vector>> gradientLocalCoefficients
                (
                    mesh.nCells()
                );
                List<DynamicList<label>> gradientGlobalColumns(mesh.nCells());
                List<DynamicList<vector>> gradientGlobalCoefficients
                (
                    mesh.nCells()
                );

                forAll(owner, sourceFaceI)
                {
                    const label sourceOwner = owner[sourceFaceI];
                    const label sourceNeighbour = neighbour[sourceFaceI];

                    gradientLocalColumns[sourceOwner].append(sourceOwner);
                    gradientLocalCoefficients[sourceOwner].append
                    (
                        -ownLs[sourceFaceI]
                    );
                    gradientLocalColumns[sourceOwner].append(sourceNeighbour);
                    gradientLocalCoefficients[sourceOwner].append
                    (
                        ownLs[sourceFaceI]
                    );

                    gradientLocalColumns[sourceNeighbour].append
                    (
                        sourceOwner
                    );
                    gradientLocalCoefficients[sourceNeighbour].append
                    (
                        neiLs[sourceFaceI]
                    );
                    gradientLocalColumns[sourceNeighbour].append
                    (
                        sourceNeighbour
                    );
                    gradientLocalCoefficients[sourceNeighbour].append
                    (
                        -neiLs[sourceFaceI]
                    );
                }

                const PtrList<labelList>* neiProcGlobalIDsPtr = nullptr;
                if (Pstream::parRun())
                {
                    neiProcGlobalIDsPtr =
                        &this->processorNeighbourGlobalIDs(mesh);
                }

                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& patch = mesh.boundary()[patchI];
                    const labelUList& faceCells = patch.faceCells();
                    const fvsPatchVectorField& patchOwnLs =
                        ownLs.boundaryField()[patchI];

                    if (patch.type() == "processor")
                    {
                        const labelList& patchNeiGlobal =
                            (*neiProcGlobalIDsPtr)[patchI];

                        forAll(faceCells, faceI)
                        {
                            const label cellI = faceCells[faceI];
                            gradientLocalColumns[cellI].append(cellI);
                            gradientLocalCoefficients[cellI].append
                            (
                                -patchOwnLs[faceI]
                            );
                            gradientGlobalColumns[cellI].append
                            (
                                patchNeiGlobal[faceI]
                            );
                            gradientGlobalCoefficients[cellI].append
                            (
                                patchOwnLs[faceI]
                            );
                        }
                    }
                    else if (patch.coupled())
                    {
                        FatalErrorInFunction
                            << "Coupled boundaries (except processors) are "
                            << "not supported by "
                            << "preconditionerLeastSquaresPressureStabilisation"
                            << exit(FatalError);
                    }
                    else if
                    (
                        isA<symmetryPolyPatch>(patch.patch())
#ifdef OPENFOAM_NOT_EXTEND
                     || isA<symmetryPlanePolyPatch>(patch.patch())
#endif
                    )
                    {
                        // The scalar least-squares gradient uses the mirrored
                        // value p_cell - p_cell on symmetry faces, hence no
                        // pressure-column contribution is added here.
                    }
                    else if (useBoundaryFaceValues[patchI])
                    {
                        forAll(faceCells, faceI)
                        {
                            const label cellI = faceCells[faceI];
                            gradientLocalColumns[cellI].append(cellI);
                            gradientLocalCoefficients[cellI].append
                            (
                                -patchOwnLs[faceI]
                            );
                        }
                    }
                }

                List<PetscScalar> pValues(blockSize_*blockSize_, 0.0);

                auto insertPressurePressure =
                [&]
                (
                    const label rowCell,
                    const label colCell,
                    const scalar coefficient
                )
                {
                    forAll(pValues, valueI)
                    {
                        pValues[valueI] = 0.0;
                    }

                    pValues
                    [
                        (blockSize_ - 1)*blockSize_ + blockSize_ - 1
                    ] = coefficient;

                    const label globalRow =
                        foamPetscSnesHelper::globalCells().toGlobal(rowCell);
                    const label globalCol =
                        foamPetscSnesHelper::globalCells().toGlobal(colCell);

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalRow,
                            1,
                            &globalCol,
                            pValues.cdata(),
                            ADD_VALUES
                        )
                    );
                };

                auto insertPressurePressureGlobal =
                [&]
                (
                    const label rowCell,
                    const label globalCol,
                    const scalar coefficient
                )
                {
                    forAll(pValues, valueI)
                    {
                        pValues[valueI] = 0.0;
                    }

                    pValues
                    [
                        (blockSize_ - 1)*blockSize_ + blockSize_ - 1
                    ] = coefficient;

                    const label globalRow =
                        foamPetscSnesHelper::globalCells().toGlobal(rowCell);

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalRow,
                            1,
                            &globalCol,
                            pValues.cdata(),
                            ADD_VALUES
                        )
                    );
                };

                // Exchange the complete gradient stencil of each processor
                // boundary cell.  A processor-face interpolated gradient
                // depends on both local and remote gradients; the remote
                // gradient itself can contain neighbour-of-neighbour and
                // further off-rank columns.  Sending the stencil, rather
                // than only the remote face-cell ID, preserves that graph.
                List<List<labelList>> remoteGradientColumns
                (
                    mesh.boundary().size()
                );
                List<List<vectorField>> remoteGradientCoefficients
                (
                    mesh.boundary().size()
                );

                if (Pstream::parRun())
                {
                    List<DynamicList<label>> procPatchIndices
                    (
                        Pstream::nProcs()
                    );

                    forAll(mesh.boundary(), patchI)
                    {
                        const fvPatch& patch = mesh.boundary()[patchI];
                        if (patch.type() == "processor")
                        {
                            procPatchIndices
                            [refCast<const processorFvPatch>(patch).neighbProcNo()]
                            .append(patchI);
                        }
                    }

                    PstreamBuffers stencilBuffers
                    (
                        Pstream::commsTypes::nonBlocking
                    );

                    forAll(procPatchIndices, procI)
                    {
                        if (procI == Pstream::myProcNo())
                        {
                            continue;
                        }

                        const DynamicList<label>& patchIndices =
                            procPatchIndices[procI];
                        if (patchIndices.empty())
                        {
                            continue;
                        }

                        UOPstream toProc(procI, stencilBuffers);
                        toProc << patchIndices.size();

                        forAll(patchIndices, patchOrdinal)
                        {
                            const label patchI = patchIndices[patchOrdinal];
                            const labelUList& faceCells =
                                mesh.boundary()[patchI].faceCells();

                            toProc << mesh.boundary()[patchI].size();

                            forAll(faceCells, faceI)
                            {
                                const label cellI = faceCells[faceI];
                                labelList globalColumns
                                (
                                    gradientLocalColumns[cellI].size()
                                  + gradientGlobalColumns[cellI].size()
                                );
                                vectorField coefficients
                                (
                                    globalColumns.size(), vector::zero
                                );

                                label next = 0;
                                forAll
                                (
                                    gradientLocalColumns[cellI],
                                    entryI
                                )
                                {
                                    globalColumns[next] =
                                        foamPetscSnesHelper::globalCells()
                                       .toGlobal
                                        (
                                            gradientLocalColumns[cellI][entryI]
                                        );
                                    coefficients[next++] =
                                        gradientLocalCoefficients[cellI]
                                        [entryI];
                                }
                                forAll
                                (
                                    gradientGlobalColumns[cellI],
                                    entryI
                                )
                                {
                                    globalColumns[next] =
                                        gradientGlobalColumns[cellI][entryI];
                                    coefficients[next++] =
                                        gradientGlobalCoefficients[cellI]
                                        [entryI];
                                }

                                toProc << globalColumns << coefficients;
                            }
                        }
                    }

                    stencilBuffers.finishedSends();

                    forAll(procPatchIndices, procI)
                    {
                        if (procI == Pstream::myProcNo())
                        {
                            continue;
                        }

                        const DynamicList<label>& patchIndices =
                            procPatchIndices[procI];
                        if (patchIndices.empty())
                        {
                            continue;
                        }

                        UIPstream fromProc(procI, stencilBuffers);
                        label nRemotePatches = 0;
                        fromProc >> nRemotePatches;

                        if (nRemotePatches != patchIndices.size())
                        {
                            FatalErrorInFunction
                                << "Processor least-squares stencil patch "
                                << "count mismatch with processor " << procI
                                << ": received " << nRemotePatches
                                << ", expected " << patchIndices.size()
                                << exit(FatalError);
                        }

                        for (label patchOrdinal = 0;
                             patchOrdinal < nRemotePatches;
                             ++patchOrdinal)
                        {
                            const label patchI = patchIndices[patchOrdinal];
                            label nFaces = 0;
                            fromProc >> nFaces;

                            const label expectedFaces =
                                mesh.boundary()[patchI].size();
                            if (nFaces != expectedFaces)
                            {
                                FatalErrorInFunction
                                    << "Processor least-squares stencil face "
                                    << "count mismatch on patch "
                                    << mesh.boundary()[patchI].name()
                                    << ": received " << nFaces
                                    << ", expected " << expectedFaces
                                    << exit(FatalError);
                            }

                            remoteGradientColumns[patchI].setSize(nFaces);
                            remoteGradientCoefficients[patchI].setSize(nFaces);
                            for (label faceI = 0; faceI < nFaces; ++faceI)
                            {
                                fromProc
                                    >> remoteGradientColumns[patchI][faceI]
                                    >> remoteGradientCoefficients[patchI]
                                    [faceI];
                            }
                        }
                    }
                }

                auto addLocalGradientStencil =
                [&]
                (
                    const label targetRow,
                    const label gradientCell,
                    const scalar interpolationWeight,
                    const scalar divergenceSign,
                    const scalar faceScale,
                    const vector& gradientFaceCoefficient
                )
                {
                    forAll
                    (
                        gradientLocalColumns[gradientCell],
                        entryI
                    )
                    {
                        const scalar coefficient =
                            divergenceSign
                           *faceScale
                           *interpolationWeight
                           *(
                                gradientFaceCoefficient
                              & gradientLocalCoefficients[gradientCell]
                                [entryI]
                            );
                        insertPressurePressure
                        (
                            targetRow,
                            gradientLocalColumns[gradientCell][entryI],
                            coefficient
                        );
                    }

                    forAll
                    (
                        gradientGlobalColumns[gradientCell],
                        entryI
                    )
                    {
                        const scalar coefficient =
                            divergenceSign
                           *faceScale
                           *interpolationWeight
                           *(
                                gradientFaceCoefficient
                              & gradientGlobalCoefficients[gradientCell]
                                [entryI]
                            );
                        insertPressurePressureGlobal
                        (
                            targetRow,
                            gradientGlobalColumns[gradientCell][entryI],
                            coefficient
                        );
                    }
                };

                auto addGlobalGradientStencil =
                [&]
                (
                    const label targetRow,
                    const labelList& columns,
                    const vectorField& coefficients,
                    const scalar interpolationWeight,
                    const scalar divergenceSign,
                    const scalar faceScale,
                    const vector& gradientFaceCoefficient
                )
                {
                    forAll(columns, entryI)
                    {
                        const scalar coefficient =
                            divergenceSign
                           *faceScale
                           *interpolationWeight
                           *(
                                gradientFaceCoefficient
                              & coefficients[entryI]
                            );
                        insertPressurePressureGlobal
                        (
                            targetRow,
                            columns[entryI],
                            coefficient
                        );
                    }
                };

                // Add the complete corrected-snGrad/interpolated-gradient
                // dependency for internal faces.
                forAll(owner, targetFaceI)
                {
                    const label targetOwner = owner[targetFaceI];
                    const label targetNeighbour = neighbour[targetFaceI];
                    const vector n = mesh.Sf()[targetFaceI]/magSf[targetFaceI];
                    const vector gradientFaceCoefficient =
                        nonOrthCorrectionVectors[targetFaceI] - n;
                    const scalar faceScale =
                        pressureStabScale
                       *gamma[targetFaceI]
                       *magSf[targetFaceI];

                    addLocalGradientStencil
                    (
                        targetOwner,
                        targetOwner,
                        weights[targetFaceI],
                        1.0,
                        faceScale,
                        gradientFaceCoefficient
                    );
                    addLocalGradientStencil
                    (
                        targetOwner,
                        targetNeighbour,
                        1.0 - weights[targetFaceI],
                        1.0,
                        faceScale,
                        gradientFaceCoefficient
                    );
                    addLocalGradientStencil
                    (
                        targetNeighbour,
                        targetOwner,
                        weights[targetFaceI],
                        -1.0,
                        faceScale,
                        gradientFaceCoefficient
                    );
                    addLocalGradientStencil
                    (
                        targetNeighbour,
                        targetNeighbour,
                        1.0 - weights[targetFaceI],
                        -1.0,
                        faceScale,
                        gradientFaceCoefficient
                    );
                }

                // Add the same dependency on physical and processor boundary
                // faces.  fvc::div contributes a boundary flux to the owner
                // cell only.  Processor-face interpolation includes the
                // communicated remote cell-gradient stencil.
                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& patch = mesh.boundary()[patchI];
                    const labelUList& faceCells = patch.faceCells();
                    const fvsPatchScalarField& patchGamma =
                        gamma.boundaryField()[patchI];
                    const fvsPatchScalarField& patchMagSf =
                        magSf.boundaryField()[patchI];
                    const fvsPatchVectorField& patchSf =
                        mesh.Sf().boundaryField()[patchI];
                    const fvsPatchVectorField& patchCorrection =
                        nonOrthCorrectionVectors.boundaryField()[patchI];
                    const fvsPatchScalarField& patchWeights =
                        weights.boundaryField()[patchI];

                    forAll(faceCells, faceI)
                    {
                        const label targetCell = faceCells[faceI];
                        const vector n =
                            patchSf[faceI]/patchMagSf[faceI];
                        const vector gradientFaceCoefficient =
                            patchCorrection[faceI] - n;
                        const scalar faceScale =
                            pressureStabScale
                           *patchGamma[faceI]
                           *patchMagSf[faceI];

                        if (patch.type() == "processor")
                        {
                            addLocalGradientStencil
                            (
                                targetCell,
                                targetCell,
                                patchWeights[faceI],
                                1.0,
                                faceScale,
                                gradientFaceCoefficient
                            );
                            addGlobalGradientStencil
                            (
                                targetCell,
                                remoteGradientColumns[patchI][faceI],
                                remoteGradientCoefficients[patchI][faceI],
                                1.0 - patchWeights[faceI],
                                1.0,
                                faceScale,
                                gradientFaceCoefficient
                            );
                        }
                        else
                        {
                            addLocalGradientStencil
                            (
                                targetCell,
                                targetCell,
                                1.0,
                                1.0,
                                faceScale,
                                gradientFaceCoefficient
                            );
                        }
                    }
                }

                Info<< "Least-squares pressure stabilisation P_pp "
                    << "full wide correction enabled; internal, boundary and "
                    << "processor gradient stencil added to the compact "
                    << "Jacobian" << endl;
            }

            // Insert D-in-p equation coefficients matching the
            // linearisation of -g(J) about J=1. Both the generic and
            // Gultekin constraints have dg/dJ=1 at J=1, so the established
            // approximate-Jacobian stencil remains unchanged. Matrix-free
            // directional derivatives use the complete nonlinear residual.
            const scalar constraintDerivativeAtReference =
                mixedVolumetricConstraintDerivative(1.0);
            //
            if (preconditionerLeastSquaresPressureCoupling_)
            {
                const fvMesh& mesh = this->mesh();
                const word gradDScheme
                (
                    mesh.gradScheme("grad(" + D.name() + ')')
                );

                if
                (
                    gradDScheme != "leastSquaresS4f"
                 && gradDScheme != "leastSquaresS4fDirichlet"
                )
                {
                    FatalErrorInFunction
                        << "preconditionerLeastSquaresPressureCoupling "
                        << "requires grad(" << D.name()
                        << ") leastSquaresS4f or leastSquaresS4fDirichlet; "
                        << "got " << gradDScheme
                        << exit(FatalError);
                }

                boolList useBoundaryFaceValues(mesh.boundary().size(), false);
                if (gradDScheme == "leastSquaresS4fDirichlet")
                {
                    forAll(useBoundaryFaceValues, patchI)
                    {
                        useBoundaryFaceValues[patchI] =
                            D.boundaryField()[patchI].fixesValue()
                         && !isA<solidTractionFvPatchVectorField>
                            (D.boundaryField()[patchI]);
                    }
                }

#ifdef OPENFOAM_COM
                word lsCacheName("leastSquaresVectors");
                if (gradDScheme == "leastSquaresS4fDirichlet")
                {
                    lsCacheName += "Dirichlet";
                }
                lsCacheName += D.name();
                const leastSquaresS4fVectors& lsv =
                    leastSquaresS4fVectors::New
                    (
                        lsCacheName, mesh, useBoundaryFaceValues
                    );
#else
                const leastSquaresS4fVectors& lsv =
                    leastSquaresS4fVectors::New
                    (
                        mesh, useBoundaryFaceValues
                    );
#endif

                const surfaceVectorField& ownLs = lsv.pVectors();
                const surfaceVectorField& neiLs = lsv.nVectors();
                const scalarField& cellJ = J_;
                const tensorField& cellFinv = Finv_;
                const scalarField& cellV = mesh.V();
                const label nDComponents = solidModel::twoD() ? 2 : 3;
                const label nCoeffCmpts = blockSize_*blockSize_;
                List<PetscScalar> values(nCoeffCmpts, 0.0);

                auto insertPressureDisplacementBlock =
                [&]
                (
                    const label rowCell,
                    const label colCell,
                    const vector& coefficient
                )
                {
                    forAll(values, valueI)
                    {
                        values[valueI] = 0.0;
                    }

                    for (label cmptI = 0; cmptI < nDComponents; ++cmptI)
                    {
                        values
                        [
                            (blockSize_ - 1)*blockSize_ + cmptI
                        ] = coefficient[cmptI];
                    }

                    const label globalRow =
                        foamPetscSnesHelper::globalCells().toGlobal(rowCell);
                    const label globalCol =
                        foamPetscSnesHelper::globalCells().toGlobal(colCell);

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalRow,
                            1,
                            &globalCol,
                            values.cdata(),
                            ADD_VALUES
                        )
                    );
                };

                auto insertPressureDisplacementGlobalBlock =
                [&]
                (
                    const label rowCell,
                    const label globalCol,
                    const vector& coefficient
                )
                {
                    forAll(values, valueI)
                    {
                        values[valueI] = 0.0;
                    }

                    for (label cmptI = 0; cmptI < nDComponents; ++cmptI)
                    {
                        values
                        [
                            (blockSize_ - 1)*blockSize_ + cmptI
                        ] = coefficient[cmptI];
                    }

                    const label globalRow =
                        foamPetscSnesHelper::globalCells().toGlobal(rowCell);

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalRow,
                            1,
                            &globalCol,
                            values.cdata(),
                            ADD_VALUES
                        )
                    );
                };

                const labelUList& own = mesh.owner();
                const labelUList& nei = mesh.neighbour();

                forAll(own, faceI)
                {
                    const label ownCell = own[faceI];
                    const label neiCell = nei[faceI];

                    const scalar ownScale =
                        pressureEqnScale_
                       *mixedVolumetricConstraintDerivative(cellJ[ownCell])
                       *cellV[ownCell]
                       *cellJ[ownCell];
                    const scalar neiScale =
                        pressureEqnScale_
                       *mixedVolumetricConstraintDerivative(cellJ[neiCell])
                       *cellV[neiCell]
                       *cellJ[neiCell];

                    const vector ownCoeff =
                        ownScale*(cellFinv[ownCell].T() & ownLs[faceI]);
                    const vector neiCoeff =
                        neiScale*(cellFinv[neiCell].T() & neiLs[faceI]);

                    // R_p contains -g(J).  The owner and neighbour
                    // least-squares gradients use D_neighbour-D_owner,
                    // with the neighbour gradient carrying the opposite
                    // sign.  Thus each cell receives a positive diagonal
                    // coefficient and a negative face-neighbour coefficient.
                    insertPressureDisplacementBlock
                    (
                        ownCell, ownCell, ownCoeff
                    );
                    insertPressureDisplacementBlock
                    (
                        ownCell, neiCell, -ownCoeff
                    );
                    insertPressureDisplacementBlock
                    (
                        neiCell, neiCell, neiCoeff
                    );
                    insertPressureDisplacementBlock
                    (
                        neiCell, ownCell, -neiCoeff
                    );
                }

                const PtrList<labelList>& neiProcGlobalIDs =
                    this->processorNeighbourGlobalIDs(mesh);

                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& patch = mesh.boundary()[patchI];
                    const labelUList& faceCells = patch.faceCells();
                    const fvsPatchVectorField& patchOwnLs =
                        ownLs.boundaryField()[patchI];

                    if (patch.type() == "processor")
                    {
                        const labelList& patchNeiGlobal =
                            neiProcGlobalIDs[patchI];

                        forAll(faceCells, faceI)
                        {
                            const label cellI = faceCells[faceI];
                            const scalar cellScale =
                                pressureEqnScale_
                               *mixedVolumetricConstraintDerivative
                                (
                                    cellJ[cellI]
                                )
                               *cellV[cellI]
                               *cellJ[cellI];
                            const vector coefficient =
                                cellScale
                               *(cellFinv[cellI].T() & patchOwnLs[faceI]);

                            insertPressureDisplacementBlock
                            (
                                cellI, cellI, coefficient
                            );
                            insertPressureDisplacementGlobalBlock
                            (
                                cellI,
                                patchNeiGlobal[faceI],
                                -coefficient
                            );
                        }
                    }
                    else if (patch.coupled())
                    {
                        FatalErrorInFunction
                            << "Coupled boundaries (except processors) are "
                            << "not supported by "
                            << "preconditionerLeastSquaresPressureCoupling"
                            << exit(FatalError);
                    }
                    else if
                    (
                        isA<symmetryPolyPatch>(patch.patch())
#ifdef OPENFOAM_NOT_EXTEND
                     || isA<symmetryPlanePolyPatch>(patch.patch())
#endif
                    )
                    {
                        const vectorField n(patch.nf());
                        forAll(faceCells, faceI)
                        {
                            const label cellI = faceCells[faceI];
                            const scalar cellScale =
                                pressureEqnScale_
                               *mixedVolumetricConstraintDerivative
                                (
                                    cellJ[cellI]
                                )
                               *cellV[cellI]
                               *cellJ[cellI];
                            const vector q =
                                cellFinv[cellI].T() & patchOwnLs[faceI];
                            const vector coefficient =
                                -cellScale
                               *((-2.0*sqr(n[faceI])) & q);
                            insertPressureDisplacementBlock
                            (
                                cellI, cellI, coefficient
                            );
                        }
                    }
                    else if (useBoundaryFaceValues[patchI])
                    {
                        forAll(faceCells, faceI)
                        {
                            const label cellI = faceCells[faceI];
                            const scalar cellScale =
                                pressureEqnScale_
                               *mixedVolumetricConstraintDerivative
                                (
                                    cellJ[cellI]
                                )
                               *cellV[cellI]
                               *cellJ[cellI];
                            const vector coefficient =
                                cellScale
                               *(cellFinv[cellI].T() & patchOwnLs[faceI]);
                            insertPressureDisplacementBlock
                            (
                                cellI, cellI, coefficient
                            );
                        }
                    }
                }

                Info<< "Least-squares pressure P_pD preconditioner block "
                    << "enabled; current J/F^-1 and direct LS stencil used"
                    << endl;
            }
            else
            {
                // InsertFvmDivU's sign convention: scale=+1 assembles
                // `-V*div(U)`. So we pass `+pressureEqnScale_` to get
                // J_pD = -alpha*V*div(D).
                foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
                (
                    p,
                    D,
                    jac,
                    blockSize_ - 1,             // row offset (p row)
                    0,                          // column offset (D columns)
                    solidModel::twoD() ? 2 : 3, // number of D components
                    pressureEqnScale_*constraintDerivativeAtReference
                                                // helper returns -V*div with +1
                );
            }

            // Insert p-in-D term. InsertFvmGrad's updated sign
            // convention: scale=+1 assembles `-V*grad(p)`. Apply the
            // pressure-unknown scale so the column corresponds to the
            // scaled unknown pHat = p/pressureUnknownScale_.
            // Assemble P_Dp with the Gauss finite-volume stencil that is
            // structurally paired with InsertFvmDivUIntoPETScMatrix above.
            // The legacy least-squares insertion omits the physical
            // traction-boundary contribution and can therefore produce an
            // essentially empty pressure-to-momentum action for a constant
            // pressure perturbation, even though the production residual
            // retains the corresponding face-force imbalance.
            foamPetscSnesHelper::InsertFvmGradPGaussIntoPETScMatrix
            (
                p,
                jac,
                0,                          // row offset
                blockSize_ - 1,             // column offset
                solidModel::twoD() ? 2 : 3, // number of D components
                pressureUnknownScale_       // scale (helper returns -V*grad with +1)
            );

            // solidTraction patches replace the physical pressure face force
            // `-p*SfCurrent` with their prescribed/trial traction.  The
            // Gauss pressure stencil above represents the pressure face
            // force before that replacement, so add the opposite derivative
            // on each traction boundary face.  This is a compact P-only
            // correction: the production residual and boundary physics are
            // unchanged.  SfCurrent is J*F^-T & Sf, matching the residual's
            // current pressure force measure for arbitrary F and J.
            const surfaceVectorField SfCurrent
            (
                fvc::interpolate(J_*Finv_.T()) & this->mesh().Sf()
            );
            const label nDComponents = solidModel::twoD() ? 2 : 3;
            const label nCoeffCmpts = blockSize_*blockSize_;
            List<PetscScalar> tractionPressureValues
            (
                nCoeffCmpts,
                0.0
            );

            forAll(D.boundaryField(), patchI)
            {
                if
                (
                    !isA<solidTractionFvPatchVectorField>
                    (
                        D.boundaryField()[patchI]
                    )
                )
                {
                    continue;
                }

                const labelUList& faceCells =
                    D.boundaryField()[patchI].patch().faceCells();
                const vectorField& patchSf =
                    SfCurrent.boundaryField()[patchI];

                forAll(faceCells, faceI)
                {
                    forAll(tractionPressureValues, valueI)
                    {
                        tractionPressureValues[valueI] = 0.0;
                    }

                    const label globalBlockRowI =
                        foamPetscSnesHelper::globalCells().toGlobal
                        (
                            faceCells[faceI]
                        );

                    for (label cmptI = 0; cmptI < nDComponents; ++cmptI)
                    {
                        tractionPressureValues
                        [
                            cmptI*blockSize_ + blockSize_ - 1
                        ] =
                            pressureUnknownScale_*patchSf[faceI][cmptI];
                    }

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalBlockRowI,
                            1,
                            &globalBlockRowI,
                            tractionPressureValues.cdata(),
                            ADD_VALUES
                        )
                    );
                }
            }
        }
    }

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        [&]() -> const fvVectorMatrix&
        {
            if (!preconditionerMaterialTangent_)
            {
                return momentumStabilisation().vectorJacobian(D, &impKf_);
            }

            const PtrList<mechanicalLaw>& mechanicalLaws =
                static_cast<const PtrList<mechanicalLaw>&>(mechanical());

            if (mechanicalLaws.size() != 1)
            {
                FatalErrorInFunction
                    << "preconditionerMaterialTangent currently requires "
                    << "exactly one mechanical law; got "
                    << mechanicalLaws.size()
                    << exit(FatalError);
            }

            List<mat66> passiveTangent(this->mesh().nFaces());
            mechanicalLaws[0].passiveMaterialTangentField(passiveTangent);

            surfaceScalarField materialKf
            (
                IOobject
                (
                    "preconditionerMaterialKf",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimPressure, 0.0)
            );

            scalar minK = VGREAT;
            scalar maxK = 0.0;
            scalar sumK = 0.0;
            label nK = 0;

            forAll(materialKf, faceI)
            {
                // A positive row-sum bound is used as the scalar diffusivity
                // of the compact Laplacian. It bounds the three momentum
                // stress rows while retaining the existing impK baseline.
                scalar tangentBound = 0.0;
                for (label rowI = 0; rowI < 3; ++rowI)
                {
                    scalar rowSum = 0.0;
                    for (label colI = 0; colI < 6; ++colI)
                    {
                        rowSum += mag(passiveTangent[faceI](rowI, colI));
                    }
                    tangentBound = max(tangentBound, rowSum);
                }

                const scalar baseline = impKf_[faceI];
                const scalar value = max(baseline, tangentBound);
                materialKf[faceI] = value;
                minK = min(minK, value);
                maxK = max(maxK, value);
                sumK += value;
                ++nK;
            }

            forAll(materialKf.boundaryField(), patchI)
            {
                scalarField& patchK = materialKf.boundaryFieldRef()[patchI];
                const mat66* tangentPtr =
                    passiveTangent.cdata()
                  + this->mesh().boundaryMesh()[patchI].start();

                forAll(patchK, faceI)
                {
                    scalar tangentBound = 0.0;
                    for (label rowI = 0; rowI < 3; ++rowI)
                    {
                        scalar rowSum = 0.0;
                        for (label colI = 0; colI < 6; ++colI)
                        {
                            rowSum += mag(tangentPtr[faceI](rowI, colI));
                        }
                        tangentBound = max(tangentBound, rowSum);
                    }

                    const scalar baseline = impKf_.boundaryField()[patchI][faceI];
                    const scalar value = max(baseline, tangentBound);
                    patchK[faceI] = value;
                    minK = min(minK, value);
                    maxK = max(maxK, value);
                    sumK += value;
                    ++nK;
                }
            }

            reduce(minK, minOp<scalar>());
            reduce(maxK, maxOp<scalar>());
            reduce(sumK, sumOp<scalar>());
            reduce(nK, sumOp<label>());
            Info<< "Passive material preconditioner Kf range = "
                << minK << " to " << maxK
                << ", mean = " << sumK/(nK > 0 ? scalar(nK) : 1.0)
                << endl;

            return momentumStabilisation().vectorJacobian
            (
                D,
                &materialKf,
                true
            );
        }()
      - rho()*fvm::d2dt2(D)
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D);
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    // Convert fvMatrix matrix to PETSc matrix
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        approxJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
    );

    if (preconditionerPassiveNominalTangent_)
    {
        if (faceStressTreatment_ != "directConstitutive")
        {
            FatalErrorInFunction
                << "preconditionerPassiveNominalTangent is only supported "
                << "for faceStressTreatment directConstitutive; got "
                << faceStressTreatment_ << exit(FatalError);
        }

        if (preconditionerMaterialTangent_)
        {
            FatalErrorInFunction
                << "preconditionerPassiveNominalTangent cannot be combined "
                << "with the diagnostic scalar preconditionerMaterialTangent"
                << exit(FatalError);
        }

        const PtrList<mechanicalLaw>& mechanicalLaws =
            static_cast<const PtrList<mechanicalLaw>&>(mechanical());

        if (mechanicalLaws.size() != 1)
        {
            FatalErrorInFunction
                << "preconditionerPassiveNominalTangent currently requires "
                << "exactly one mechanical law; got "
                << mechanicalLaws.size() << exit(FatalError);
        }

        List<tensor> passiveNominalTangent;
        if
        (
            !mechanicalLaws[0].passiveNominalTangentField
            (
                passiveNominalTangent
            )
        )
        {
            FatalErrorInFunction
                << "preconditionerPassiveNominalTangent was requested, but "
                << mechanicalLaws[0].type()
                << " does not provide a passive nominal-force tangent"
                << exit(FatalError);
        }

        addDirectNominalTangentToPETSc
        (
            jac,
            passiveNominalTangent,
            this->mesh(),
            blockSize_,
            solidModel::twoD(),
            foamPetscSnesHelper::globalCells(),
            false,
            "Direct passive nominal P_DD tangent"
        );
    }

    if (preconditionerViscousTangent_)
    {
        if (faceStressTreatment_ != "directConstitutive")
        {
            FatalErrorInFunction
                << "preconditionerViscousTangent is only supported for "
                << "faceStressTreatment directConstitutive; got "
                << faceStressTreatment_ << exit(FatalError);
        }

        const PtrList<mechanicalLaw>& mechanicalLaws =
            static_cast<const PtrList<mechanicalLaw>&>(mechanical());

        if (mechanicalLaws.size() != 1)
        {
            FatalErrorInFunction
                << "preconditionerViscousTangent currently requires "
                << "exactly one mechanical law; got "
                << mechanicalLaws.size() << exit(FatalError);
        }

        List<tensor> viscousNominalTangent;
        if
        (
            !mechanicalLaws[0].viscousNominalTangentField
            (
                viscousNominalTangent
            )
         || viscousNominalTangent.size() != 9*this->mesh().nFaces()
        )
        {
            FatalErrorInFunction
                << "preconditionerViscousTangent was requested, but the "
                << mechanicalLaws[0].type()
                << " mechanical law does not provide a complete viscous "
                << "nominal-force tangent" << exit(FatalError);
        }

        const vectorField& Sf = this->mesh().Sf();
        const scalarField& deltaCoeffs = this->mesh().deltaCoeffs();
        const label nDComponents = solidModel::twoD() ? 2 : 3;
        const label nCoeffCmpts = blockSize_*blockSize_;
        List<PetscScalar> tangentValues(nCoeffCmpts, 0.0);
        scalar minFaceNorm = GREAT;
        scalar maxFaceNorm = 0.0;
        scalar sumFaceNorm = 0.0;
        label nFaces = 0;

        for (label faceI = 0; faceI < this->mesh().nInternalFaces(); ++faceI)
        {
            const scalar magSf = mag(Sf[faceI]);
            if (magSf <= VSMALL)
            {
                continue;
            }

            const vector n = Sf[faceI]/magSf;
            const scalar deltaCoeff = deltaCoeffs[faceI];
            tensor K(tensor::zero);

            // Use the same owner-neighbour normal-gradient stencil as the
            // compact finite-volume matrix, while retaining all output and
            // displacement-component coupling in the face tangent.
            for (label rowI = 0; rowI < nDComponents; ++rowI)
            {
                for (label colI = 0; colI < nDComponents; ++colI)
                {
                    scalar coefficient = 0.0;

                    for (label gradDirectionI = 0;
                         gradDirectionI < 3;
                         ++gradDirectionI)
                    {
                        const tensor& dPdF =
                            viscousNominalTangent
                            [9*faceI + 3*colI + gradDirectionI];
                        coefficient +=
                            (dPdF & Sf[faceI])[rowI]
                           *n[gradDirectionI]*deltaCoeff;
                    }

                    K(rowI, colI) = coefficient;
                }
            }

            const scalar faceNorm = mag(K);
            minFaceNorm = min(minFaceNorm, faceNorm);
            maxFaceNorm = max(maxFaceNorm, faceNorm);
            sumFaceNorm += faceNorm;
            ++nFaces;

            const label globalOwner =
                foamPetscSnesHelper::globalCells().toGlobal
                (
                    this->mesh().owner()[faceI]
                );
            const label globalNeighbour =
                foamPetscSnesHelper::globalCells().toGlobal
                (
                    this->mesh().neighbour()[faceI]
                );

            const label globalBlocks[2] =
            {
                globalOwner, globalNeighbour
            };

            for (label rowBlockI = 0; rowBlockI < 2; ++rowBlockI)
            {
                for (label colBlockI = 0; colBlockI < 2; ++colBlockI)
                {
                    forAll(tangentValues, valueI)
                    {
                        tangentValues[valueI] = 0.0;
                    }

                    const scalar blockSign =
                        rowBlockI == colBlockI ? -1.0 : 1.0;

                    for (label rowI = 0; rowI < nDComponents; ++rowI)
                    {
                        for (label colI = 0; colI < nDComponents; ++colI)
                        {
                            tangentValues
                            [
                                rowI*blockSize_ + colI
                            ] = blockSign*K(rowI, colI);
                        }
                    }

                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac,
                            1,
                            &globalBlocks[rowBlockI],
                            1,
                            &globalBlocks[colBlockI],
                            tangentValues.cdata(),
                            ADD_VALUES
                        )
                    );
                }
            }
        }

        reduce(minFaceNorm, minOp<scalar>());
        reduce(maxFaceNorm, maxOp<scalar>());
        reduce(sumFaceNorm, sumOp<scalar>());
        reduce(nFaces, sumOp<label>());
        Info<< "Direct viscous nominal P_DD tangent: faces = " << nFaces
            << ", norm range = " << minFaceNorm << " to " << maxFaceNorm
            << ", mean = "
            << sumFaceNorm/(nFaces > 0 ? scalar(nFaces) : 1.0)
            << ", internal owner-neighbour stencil only" << endl;
    }

    if (preconditionerBoundaryTangent_)
    {
        // The production residual contains the reference-area boundary
        // force directly through enforceTractionBoundaries().  Its
        // displacement derivative is therefore a block-local contribution
        // to P_DD, not a second physical force.  Calculate the current
        // coefficient of fvc::ddt(D) from the configured Euler/backward
        // scheme, without constructing another matrix or changing any
        // boundary-field state during Jacobian assembly.
        const dictionary& ddtSchemes =
            this->mesh().schemesDict().subDict("ddtSchemes");
        const word ddtScheme(ddtSchemes.lookup("default"));
        const scalar deltaT = this->mesh().time().deltaTValue();
        scalar c0 = 1.0/deltaT;

        if (ddtScheme == "backward")
        {
            const volVectorField& Dref = D;
            const bool haveSecondHistory =
                Dref.oldTime().timeIndex()
             != Dref.oldTime().oldTime().timeIndex();

            if (haveSecondHistory)
            {
                const scalar deltaT0 =
                    this->mesh().time().deltaT0Value();
                const scalar coefft =
                    1.0 + deltaT/(deltaT + deltaT0);
                c0 = coefft/deltaT;
            }
        }
        else if (ddtScheme != "Euler")
        {
            FatalErrorInFunction
                << "preconditionerBoundaryTangent supports only Euler and "
                << "backward ddt schemes; got " << ddtScheme
                << exit(FatalError);
        }

        const label nDComponents = solidModel::twoD() ? 2 : 3;
        const label nCoeffCmpts = blockSize_*blockSize_;
        List<PetscScalar> boundaryTangentValues
        (
            nCoeffCmpts,
            0.0
        );
        scalar minCoefficient = VGREAT;
        scalar maxCoefficient = 0.0;
        label nFaces = 0;

        forAll(D.boundaryField(), patchI)
        {
            if
            (
                !isA<arosticaSpringDashpotTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                )
            )
            {
                continue;
            }

            const arosticaSpringDashpotTractionFvPatchVectorField& patchBc =
                refCast<const arosticaSpringDashpotTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const bool normalCondition =
                isA<arosticaNormalSpringDashpotTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );
            const bool vectorCondition =
                isA<arosticaVectorSpringDashpotTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            if (!normalCondition && !vectorCondition)
            {
                FatalErrorInFunction
                    << "preconditionerBoundaryTangent requires an Aróstica "
                    << "normal or vector spring-dashpot patch; patch "
                    << D.boundaryField()[patchI].patch().name()
                    << " has type " << D.boundaryField()[patchI].type()
                    << exit(FatalError);
            }

            const labelUList& faceCells =
                D.boundaryField()[patchI].patch().faceCells();
            // fvPatch::nf() returns a tmp field; retain an owning copy for
            // the duration of the insertion loop.
            const vectorField n(D.boundaryField()[patchI].patch().nf());
            const scalarField& area =
                D.mesh().boundary()[patchI].magSf();

            forAll(faceCells, faceI)
            {
                const label cellI = faceCells[faceI];
                const scalar coefficient =
                    patchBc.springCoefficient().value()
                  + patchBc.dashpotCoefficient().value()*c0;

                minCoefficient = min(minCoefficient, coefficient);
                maxCoefficient = max(maxCoefficient, coefficient);
                ++nFaces;

                for (label rowI = 0; rowI < nDComponents; ++rowI)
                {
                    for (label colI = 0; colI < nDComponents; ++colI)
                    {
                        const scalar directionFactor =
                            normalCondition
                          ? n[faceI][rowI]*n[faceI][colI]
                          : (rowI == colI ? 1.0 : 0.0);

                        boundaryTangentValues
                        [
                            rowI*blockSize_ + colI
                        ] = -area[faceI]*coefficient*directionFactor;
                    }
                }

                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(cellI);

                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac,
                        1,
                        &globalBlockRowI,
                        1,
                        &globalBlockRowI,
                        boundaryTangentValues.cdata(),
                        ADD_VALUES
                    )
                );

                forAll(boundaryTangentValues, valueI)
                {
                    boundaryTangentValues[valueI] = 0.0;
                }
            }
        }

        reduce(minCoefficient, minOp<scalar>());
        reduce(maxCoefficient, maxOp<scalar>());
        reduce(nFaces, sumOp<label>());
        Info<< "Aróstica boundary preconditioner tangent coefficient range = "
            << minCoefficient << " to " << maxCoefficient
            << ", faces = " << nFaces << endl;
    }

    if (solvePressure() && pressureRowScaling_ == "volumeRmsForce")
    {
        // All pressure-row blocks have already been assembled with the
        // legacy alpha*V scaling. Apply the same volume-RMS row factor used
        // by formResidual to P_pp and P_pD in one operation, preventing any
        // block from silently receiving a different scaling.
        AssertPETSc(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        AssertPETSc(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

        Vec leftScale;
        AssertPETSc(VecDuplicate(x, &leftScale));
        AssertPETSc(VecSet(leftScale, 1.0));

        PetscScalar* values = nullptr;
        PetscInt localSize = 0;
        AssertPETSc(VecGetLocalSize(leftScale, &localSize));
        AssertPETSc(VecGetArray(leftScale, &values));

        const scalar Vtot = gSum(mesh().V());
        const scalar L0Sqr = sqr(cbrt(Vtot));
        forAll(mesh().V(), cellI)
        {
            const label localPressureI = cellI*blockSize_ + blockSize_ - 1;
            if (localPressureI >= localSize)
            {
                FatalErrorInFunction
                    << "Pressure row index " << localPressureI
                    << " exceeds local PETSc vector size " << localSize
                    << abort(FatalError);
            }
            values[localPressureI] =
                L0Sqr/sqrt(Vtot*mesh().V()[cellI]);
        }

        AssertPETSc(VecRestoreArray(leftScale, &values));
        AssertPETSc(VecAssemblyBegin(leftScale));
        AssertPETSc(VecAssemblyEnd(leftScale));
        AssertPETSc(MatDiagonalScale(jac, leftScale, nullptr));
        AssertPETSc(VecDestroy(&leftScale));
    }

    return 0;
}

#endif // USE_PETSC

tmp<vectorField> nonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
