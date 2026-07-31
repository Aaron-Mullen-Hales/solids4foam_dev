/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Application
    ArosticaHolzapfelOgdenElasticSmoke

Description
    Instantiates the Phase 1 Aróstica law from a runtime dictionary and
    exercises both its cell and direct face stress paths.

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
        solidModel("arosticaElasticSmoke", runTime)
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


void checkUniformStress
(
    const volSymmTensorField& cellStress,
    const surfaceSymmTensorField& faceStress
)
{
    scalar maximumCellFaceDifference = 0.0;
    scalar maximumCellMagnitude = 0.0;
    scalar maximumFaceMagnitude = 0.0;

    forAll(cellStress, cellI)
    {
        maximumCellMagnitude =
            max(maximumCellMagnitude, mag(cellStress[cellI]));
    }

    forAll(cellStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = cellStress.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            maximumCellMagnitude = max(maximumCellMagnitude, mag(patch[faceI]));
        }
    }

    forAll(faceStress, faceI)
    {
        maximumFaceMagnitude =
            max(maximumFaceMagnitude, mag(faceStress[faceI]));
        maximumCellFaceDifference =
            max(maximumCellFaceDifference, mag(faceStress[faceI] - cellStress[0]));
    }

    forAll(faceStress.boundaryField(), patchI)
    {
        const symmTensorField& patch = faceStress.boundaryField()[patchI];

        forAll(patch, faceI)
        {
            maximumFaceMagnitude =
                max(maximumFaceMagnitude, mag(patch[faceI]));
            maximumCellFaceDifference =
                max
                (
                    maximumCellFaceDifference,
                    mag(patch[faceI] - cellStress[0])
                );
        }
    }

    const scalar tolerance =
        1e-10*max(scalar(1.0), maximumCellMagnitude);

    if
    (
        maximumCellMagnitude <= SMALL
     || maximumFaceMagnitude <= SMALL
     || maximumCellFaceDifference > tolerance
    )
    {
        FatalErrorInFunction
            << "Cell/direct-face constitutive path mismatch" << nl
            << "    maximum cell stress = " << maximumCellMagnitude << nl
            << "    maximum face stress = " << maximumFaceMagnitude << nl
            << "    maximum cell/face difference = "
            << maximumCellFaceDifference << nl
            << "    tolerance = " << tolerance
            << abort(FatalError);
    }

    Info<< "PASS: cell/direct-face constitutive paths" << nl
        << "    maximum cell stress = " << maximumCellMagnitude << nl
        << "    maximum face stress = " << maximumFaceMagnitude << nl
        << "    maximum cell/face difference = "
        << maximumCellFaceDifference << endl;
}


void checkPassiveNominalTangent
(
    mechanicalLaw& law,
    surfaceTensorField& gradDf,
    surfaceSymmTensorField& faceStress
)
{
    const label patchI = 0;
    const label patchFaceI = 0;
    const label faceI = gradDf.mesh().nInternalFaces();
    const tensor baseGrad
    (
        0.0, 0.15, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    );
    gradDf.boundaryFieldRef()[patchI][patchFaceI] = baseGrad;

    law.correct(faceStress);
    const tensor baseF =
        law.faceDeformationGradient().boundaryField()[patchI][patchFaceI];
    const tensor baseP =
        det(baseF)
       *(faceStress.boundaryField()[patchI][patchFaceI] & inv(baseF).T());

    List<tensor> tangent;
    if (!law.passiveNominalTangentField(tangent) || tangent.size() <= faceI*9 + 8)
    {
        FatalErrorInFunction
            << "Passive nominal tangent was not provided by " << law.type()
            << abort(FatalError);
    }

    const scalar epsilon = 1.0e-6;
    tensor direction
    (
        1.0, -0.7, 0.3,
        0.2, 0.5, -0.4,
        -0.1, 0.6, 0.8
    );
    const tensor gradPlus = baseGrad + epsilon*direction.T();
    const tensor gradMinus = baseGrad - epsilon*direction.T();

    gradDf.boundaryFieldRef()[patchI][patchFaceI] = gradPlus;
    law.correct(faceStress);
    const tensor Fplus =
        law.faceDeformationGradient().boundaryField()[patchI][patchFaceI];
    const tensor Pplus =
        det(Fplus)
       *(faceStress.boundaryField()[patchI][patchFaceI] & inv(Fplus).T());

    gradDf.boundaryFieldRef()[patchI][patchFaceI] = gradMinus;
    law.correct(faceStress);
    const tensor Fminus =
        law.faceDeformationGradient().boundaryField()[patchI][patchFaceI];
    const tensor Pminus =
        det(Fminus)
       *(faceStress.boundaryField()[patchI][patchFaceI] & inv(Fminus).T());

    gradDf.boundaryFieldRef()[patchI][patchFaceI] = baseGrad;
    law.correct(faceStress);

    tensor tangentAction(tensor::zero);
    for (label componentI = 0; componentI < 9; ++componentI)
    {
        tangentAction += tangent[9*faceI + componentI]*direction[componentI];
    }
    const tensor fdAction = (Pplus - Pminus)/(2.0*epsilon);
    const scalar error = mag(tangentAction - fdAction)
       /(mag(fdAction) + VSMALL);

    if (error > 1e-5 || mag(baseP) <= SMALL)
    {
        FatalErrorInFunction
            << "Passive nominal-force tangent mismatch"
            << " relative error = " << error
            << ", prestressed nominal force = " << mag(baseP)
            << abort(FatalError);
    }

    Info<< "PASS: prestressed passive nominal-force tangent" << nl
        << "    nominal-force magnitude = " << mag(baseP) << nl
        << "    directional FD relative error = " << error << endl;
}

} // End anonymous namespace


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
            "lawProperties.arostica",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<mechanicalLaw> law = mechanicalLaw::NewNonLinGeomMechLaw
    (
        "arosticaElasticSmoke",
        mesh,
        lawDictionary,
        nonLinearGeometry::TOTAL_LAGRANGIAN
    );

    if (law->type() != "ArosticaHolzapfelOgdenElastic")
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
            "grad(D)f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "gradDf",
            dimless,
            tensor
            (
                0.0, 0.15, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        )
    );

    volSymmTensorField cellStress
    (
        IOobject
        (
            "cellStress",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    surfaceSymmTensorField faceStress
    (
        IOobject
        (
            "faceStress",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    law->correct(cellStress);
    law->correct(faceStress);
    checkUniformStress(cellStress, faceStress);
    checkPassiveNominalTangent(*law, gradDf, faceStress);

    const tmp<volScalarField> tBulk = law->bulkModulus();
    const tmp<volScalarField> tImpK = law->impK();

    if
    (
        tBulk().dimensions() != dimPressure
     || tImpK().dimensions() != dimPressure
     || tBulk()[0] <= 0.0
     || tImpK()[0] <= 0.0
    )
    {
        FatalErrorInFunction
            << "Invalid bulkModulus() or impK() returned by law"
            << abort(FatalError);
    }

    Info<< "PASS: runtime selection and modulus APIs for "
        << law->type() << nl
        << "    bulkModulus = " << tBulk()[0] << nl
        << "    impK = " << tImpK()[0] << endl;

    return 0;
}


// ************************************************************************* //
