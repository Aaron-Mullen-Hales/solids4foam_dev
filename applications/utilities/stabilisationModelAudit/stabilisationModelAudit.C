/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | solids4foam: solid mechanics and fluid-solid
   \\    /   O peration     | interactions toolbox for OpenFOAM
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    stabilisationModelAudit

Description
    Directly exercises a runtime-selected stabilisationModel on synthetic
    displacement fields without solving the nonlinear momentum equations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "stabilisationModel.H"
#include "RhieChowHighPassStab.H"
#include "compatibilityFunctions.H"

#include <cmath>

using namespace Foam;

namespace
{
    template<class Type>
    bool optionReadIfPresent(const argList& args, const word& opt, Type& value)
    {
    #ifdef OPENFOAM_COM
        return args.readIfPresent(opt, value);
    #else
        return args.optionReadIfPresent(opt, value);
    #endif
    }


    struct Norms
    {
        scalar maxValue;
        scalar rmsValue;
        label count;
    };


    Norms vectorNorms(const vectorField& field)
    {
        scalar maxValue = 0.0;
        scalar sumSqr = 0.0;
        label count = field.size();

        if (field.size())
        {
            maxValue = max(mag(field));
            sumSqr = sum(magSqr(field));
        }

        reduce(maxValue, maxOp<scalar>());
        reduce(sumSqr, sumOp<scalar>());
        reduce(count, sumOp<label>());

        return Norms
        {
            maxValue,
            count ? Foam::sqrt(sumSqr/scalar(count)) : 0.0,
            count
        };
    }


    Norms boundaryVectorNorms(const surfaceVectorField& field)
    {
        scalar maxValue = 0.0;
        scalar sumSqr = 0.0;
        label count = 0;

        const surfaceVectorField::Boundary& bf = field.boundaryField();

        forAll(bf, patchI)
        {
            const vectorField& patchField = bf[patchI];
            count += patchField.size();

            if (patchField.size())
            {
                maxValue = max(maxValue, max(mag(patchField)));
                sumSqr += sum(magSqr(patchField));
            }
        }

        reduce(maxValue, maxOp<scalar>());
        reduce(sumSqr, sumOp<scalar>());
        reduce(count, sumOp<label>());

        return Norms
        {
            maxValue,
            count ? Foam::sqrt(sumSqr/scalar(count)) : 0.0,
            count
        };
    }


    vector globalMin(const vectorField& field)
    {
        vector minValue(GREAT, GREAT, GREAT);

        forAll(field, i)
        {
            minValue.x() = min(minValue.x(), field[i].x());
            minValue.y() = min(minValue.y(), field[i].y());
            minValue.z() = min(minValue.z(), field[i].z());
        }

        reduce(minValue.x(), minOp<scalar>());
        reduce(minValue.y(), minOp<scalar>());
        reduce(minValue.z(), minOp<scalar>());

        return minValue;
    }


    vector globalMax(const vectorField& field)
    {
        vector maxValue(-GREAT, -GREAT, -GREAT);

        forAll(field, i)
        {
            maxValue.x() = max(maxValue.x(), field[i].x());
            maxValue.y() = max(maxValue.y(), field[i].y());
            maxValue.z() = max(maxValue.z(), field[i].z());
        }

        reduce(maxValue.x(), maxOp<scalar>());
        reduce(maxValue.y(), maxOp<scalar>());
        reduce(maxValue.z(), maxOp<scalar>());

        return maxValue;
    }


    vector safeSpan(const vector& minC, const vector& maxC)
    {
        vector span(maxC - minC);

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (mag(span[cmpt]) < SMALL)
            {
                span[cmpt] = 1.0;
            }
        }

        return span;
    }


    wordList auditPatchTypes(const fvMesh& mesh, const word& physicalPatchType)
    {
        wordList patchTypes(mesh.boundary().size(), physicalPatchType);

        forAll(patchTypes, patchI)
        {
            if (mesh.boundary()[patchI].coupled())
            {
                patchTypes[patchI] = mesh.boundary()[patchI].type();
            }
        }

        return patchTypes;
    }


    tensor rigidRotationGradient(const vector& omega)
    {
        return tensor
        (
            0.0,        omega.z(), -omega.y(),
           -omega.z(),  0.0,        omega.x(),
            omega.y(), -omega.x(),  0.0
        );
    }


    vector checkerboardValue
    (
        const vector& x,
        const vector& minC,
        const vector& span,
        const vector& cells,
        const vector& amplitude
    )
    {
        label parity = 0;

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            const scalar nCells = max(cells[cmpt], scalar(1));
            const scalar xi = (x[cmpt] - minC[cmpt])/span[cmpt];
            const label index =
                min(max(label(std::floor(nCells*xi)), label(0)), label(nCells - 1));

            parity += index;
        }

        const scalar sign = (parity % 2) ? -1.0 : 1.0;

        return sign*amplitude;
    }


    vector smoothValue
    (
        const vector& x,
        const vector& minC,
        const vector& span,
        const vector& amplitude,
        const vector& frequency
    )
    {
        const scalar twoPi = 2.0*constant::mathematical::pi;
        const scalar xi = (x.x() - minC.x())/span.x();
        const scalar eta = (x.y() - minC.y())/span.y();
        const scalar zeta = (x.z() - minC.z())/span.z();

        const scalar kx = twoPi*frequency.x();
        const scalar ky = twoPi*frequency.y();
        const scalar kz = twoPi*frequency.z();

        return vector
        (
            amplitude.x()*Foam::sin(kx*xi)*Foam::sin(ky*eta),
            amplitude.y()*Foam::cos(ky*eta)*Foam::sin(kz*zeta),
            amplitude.z()*Foam::sin(kz*zeta)*Foam::cos(kx*xi)
        );
    }


    vector quadraticValue
    (
        const vector& x,
        const vector& minC,
        const vector& span,
        const vector& amplitude
    )
    {
        const scalar xi = (x.x() - minC.x())/span.x();
        const scalar eta = (x.y() - minC.y())/span.y();
        const scalar zeta = (x.z() - minC.z())/span.z();

        return vector
        (
            amplitude.x()
           *(
                1.0 + 0.7*xi - 0.4*eta + 0.2*zeta
              + 0.3*sqr(xi) - 0.2*sqr(eta) + 0.15*sqr(zeta)
              + 0.25*xi*eta - 0.17*xi*zeta + 0.11*eta*zeta
            ),
            amplitude.y()
           *(
               -0.3 + 0.2*xi + 0.8*eta - 0.5*zeta
              -0.12*sqr(xi) + 0.27*sqr(eta) - 0.19*sqr(zeta)
              + 0.09*xi*eta + 0.21*xi*zeta - 0.13*eta*zeta
            ),
            amplitude.z()
           *(
                0.4 - 0.6*xi + 0.1*eta + 0.5*zeta
              + 0.22*sqr(xi) + 0.16*sqr(eta) - 0.31*sqr(zeta)
              - 0.14*xi*eta + 0.18*xi*zeta + 0.07*eta*zeta
            )
        );
    }


    tensor quadraticGradient
    (
        const vector& x,
        const vector& minC,
        const vector& span,
        const vector& amplitude
    )
    {
        const scalar xi = (x.x() - minC.x())/span.x();
        const scalar eta = (x.y() - minC.y())/span.y();
        const scalar zeta = (x.z() - minC.z())/span.z();

        const scalar dDxdx =
            amplitude.x()*(0.7 + 0.6*xi + 0.25*eta - 0.17*zeta)/span.x();
        const scalar dDxdy =
            amplitude.x()*(-0.4 - 0.4*eta + 0.25*xi + 0.11*zeta)/span.y();
        const scalar dDxdz =
            amplitude.x()*(0.2 + 0.3*zeta - 0.17*xi + 0.11*eta)/span.z();

        const scalar dDydx =
            amplitude.y()*(0.2 - 0.24*xi + 0.09*eta + 0.21*zeta)/span.x();
        const scalar dDydy =
            amplitude.y()*(0.8 + 0.54*eta + 0.09*xi - 0.13*zeta)/span.y();
        const scalar dDydz =
            amplitude.y()*(-0.5 - 0.38*zeta + 0.21*xi - 0.13*eta)/span.z();

        const scalar dDzdx =
            amplitude.z()*(-0.6 + 0.44*xi - 0.14*eta + 0.18*zeta)/span.x();
        const scalar dDzdy =
            amplitude.z()*(0.1 + 0.32*eta - 0.14*xi + 0.07*zeta)/span.y();
        const scalar dDzdz =
            amplitude.z()*(0.5 - 0.62*zeta + 0.18*xi + 0.07*eta)/span.z();

        return tensor
        (
            dDxdx, dDydx, dDzdx,
            dDxdy, dDydy, dDzdy,
            dDxdz, dDydz, dDzdz
        );
    }


    tensor smoothGradient
    (
        const vector& x,
        const vector& minC,
        const vector& span,
        const vector& amplitude,
        const vector& frequency
    )
    {
        const scalar twoPi = 2.0*constant::mathematical::pi;
        const scalar xi = (x.x() - minC.x())/span.x();
        const scalar eta = (x.y() - minC.y())/span.y();
        const scalar zeta = (x.z() - minC.z())/span.z();

        const scalar kx = twoPi*frequency.x();
        const scalar ky = twoPi*frequency.y();
        const scalar kz = twoPi*frequency.z();

        const scalar dDxdx =
            amplitude.x()*Foam::cos(kx*xi)*(kx/span.x())*Foam::sin(ky*eta);
        const scalar dDxdy =
            amplitude.x()*Foam::sin(kx*xi)*Foam::cos(ky*eta)*(ky/span.y());
        const scalar dDxdz = 0.0;

        const scalar dDydx = 0.0;
        const scalar dDydy =
           -amplitude.y()*Foam::sin(ky*eta)*(ky/span.y())*Foam::sin(kz*zeta);
        const scalar dDydz =
            amplitude.y()*Foam::cos(ky*eta)*Foam::cos(kz*zeta)*(kz/span.z());

        const scalar dDzdx =
           -amplitude.z()*Foam::sin(kz*zeta)*Foam::sin(kx*xi)*(kx/span.x());
        const scalar dDzdy = 0.0;
        const scalar dDzdz =
            amplitude.z()*Foam::cos(kz*zeta)*(kz/span.z())*Foam::cos(kx*xi);

        return tensor
        (
            dDxdx, dDydx, dDzdx,
            dDxdy, dDydy, dDzdy,
            dDxdz, dDydz, dDzdz
        );
    }


    void setSyntheticFields
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        volVectorField& D,
        volTensorField& gradD
    )
    {
        const vector b(dict.lookupOrDefault<vector>("translation", vector::zero));
        const vector omega(dict.lookupOrDefault<vector>("omega", vector(0.04, -0.03, 0.02)));
        const tensor affineTensor
        (
            dict.lookupOrDefault<tensor>
            (
                "affineTensor",
                tensor(0.07, -0.05, 0.04, -0.02, -0.03, 0.06, 0.05, 0.01, 0.02)
            )
        );
        const vector smoothAmplitude
        (
            dict.lookupOrDefault<vector>("smoothAmplitude", vector(0.02, -0.015, 0.01))
        );
        const vector smoothFrequency
        (
            dict.lookupOrDefault<vector>("smoothFrequency", vector::one)
        );
        const vector checkerboardAmplitude
        (
            dict.lookupOrDefault<vector>("checkerboardAmplitude", vector(0.01, -0.006, 0.004))
        );
        const vector checkerboardCells
        (
            dict.lookupOrDefault<vector>("checkerboardCells", vector(8, 8, 8))
        );

        const vector minC = globalMin(mesh.C());
        const vector maxC = globalMax(mesh.C());
        const vector span = safeSpan(minC, maxC);

        tensor affineGrad(tensor::zero);

        if (fieldName == "translation")
        {
            affineGrad = tensor::zero;
        }
        else if (fieldName == "rigidRotation")
        {
            affineGrad = rigidRotationGradient(omega);
        }
        else if (fieldName == "affine")
        {
            affineGrad = affineTensor;
        }

        vectorField& DI = primitiveFieldRef(D);
        tensorField& gradDI = primitiveFieldRef(gradD);

        forAll(DI, cellI)
        {
            const vector& C = mesh.C()[cellI];

            if (fieldName == "translation")
            {
                DI[cellI] = b;
                gradDI[cellI] = tensor::zero;
            }
            else if (fieldName == "rigidRotation" || fieldName == "affine")
            {
                DI[cellI] = (C & affineGrad) + b;
                gradDI[cellI] = affineGrad;
            }
            else if (fieldName == "smooth")
            {
                DI[cellI] = smoothValue(C, minC, span, smoothAmplitude, smoothFrequency);
                gradDI[cellI] =
                    smoothGradient(C, minC, span, smoothAmplitude, smoothFrequency);
            }
            else if (fieldName == "quadratic")
            {
                DI[cellI] =
                    quadraticValue(C, minC, span, smoothAmplitude) + b;
                gradDI[cellI] =
                    quadraticGradient(C, minC, span, smoothAmplitude);
            }
            else if (fieldName == "checkerboard")
            {
                DI[cellI] =
                    b
                  + checkerboardValue
                    (
                        C,
                        minC,
                        span,
                        checkerboardCells,
                        checkerboardAmplitude
                    );
                gradDI[cellI] = tensor::zero;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown synthetic field " << fieldName << nl
                    << "Valid fields are translation, rigidRotation, affine, "
                    << "quadratic, smooth and checkerboard"
                    << abort(FatalError);
            }
        }

        D.correctBoundaryConditions();
        gradD.correctBoundaryConditions();

        volVectorField::Boundary& DB = D.boundaryFieldRef();
        volTensorField::Boundary& gradDB = gradD.boundaryFieldRef();

        forAll(DB, patchI)
        {
            if (DB[patchI].coupled())
            {
                continue;
            }

            const vectorField& Cf = mesh.Cf().boundaryField()[patchI];

            forAll(DB[patchI], faceI)
            {
                const vector& C = Cf[faceI];

                if (fieldName == "translation")
                {
                    DB[patchI][faceI] = b;
                    gradDB[patchI][faceI] = tensor::zero;
                }
                else if (fieldName == "rigidRotation" || fieldName == "affine")
                {
                    DB[patchI][faceI] = (C & affineGrad) + b;
                    gradDB[patchI][faceI] = affineGrad;
                }
                else if (fieldName == "smooth")
                {
                    DB[patchI][faceI] =
                        smoothValue(C, minC, span, smoothAmplitude, smoothFrequency);
                    gradDB[patchI][faceI] =
                        smoothGradient(C, minC, span, smoothAmplitude, smoothFrequency);
                }
                else if (fieldName == "quadratic")
                {
                    DB[patchI][faceI] =
                        quadraticValue(C, minC, span, smoothAmplitude) + b;
                    gradDB[patchI][faceI] =
                        quadraticGradient(C, minC, span, smoothAmplitude);
                }
                else if (fieldName == "checkerboard")
                {
                    DB[patchI][faceI] =
                        b
                      + checkerboardValue
                        (
                            C,
                            minC,
                            span,
                            checkerboardCells,
                            checkerboardAmplitude
                        );
                    gradDB[patchI][faceI] = tensor::zero;
                }
            }
        }
    }


    void report
    (
        const word& modelType,
        const word& fieldName,
        const word& gradientMode,
        const stabilisationModel& stab,
        const surfaceVectorField& faceVector,
        const volVectorField& cellVector
    )
    {
        const Norms internalNorms = vectorNorms(faceVector.primitiveField());
        const Norms boundaryNorms = boundaryVectorNorms(faceVector);
        const Norms cellNorms = vectorNorms(primitiveField(cellVector));

        vector netCellForce = vector::zero;
        scalar sumMagCellForce = 0.0;

        forAll(cellVector, cellI)
        {
            const vector extensive =
                cellVector[cellI]*cellVector.mesh().V()[cellI];

            netCellForce += extensive;
            sumMagCellForce += mag(extensive);
        }

        reduce(netCellForce, sumOp<vector>());
        reduce(sumMagCellForce, sumOp<scalar>());

        scalar highMax = -1.0;
        scalar highRms = -1.0;
        label highPassFallbacks = -1;
        scalar highPassReproductionMax = -1.0;

        if (isA<RhieChowHighPassStab>(stab))
        {
            const RhieChowHighPassStab& highPass =
                refCast<const RhieChowHighPassStab>(stab);

            const Norms highNorms =
                vectorNorms(primitiveField(highPass.highField()));

            highMax = highNorms.maxValue;
            highRms = highNorms.rmsValue;
            highPassFallbacks = highPass.nOrderFallbacks();
            highPassReproductionMax = highPass.maxReproductionError();
        }

        Info<< "auditResult"
            << " model=" << modelType
            << " field=" << fieldName
            << " gradient=" << gradientMode
            << " internalMax=" << internalNorms.maxValue
            << " internalRms=" << internalNorms.rmsValue
            << " boundaryMax=" << boundaryNorms.maxValue
            << " boundaryRms=" << boundaryNorms.rmsValue
            << " cellMax=" << cellNorms.maxValue
            << " cellRms=" << cellNorms.rmsValue
            << " netCellForce=" << netCellForce
            << " sumMagCellForce=" << sumMagCellForce
            << " highMax=" << highMax
            << " highRms=" << highRms
            << " highPassFallbacks=" << highPassFallbacks
            << " highPassReproductionMax=" << highPassReproductionMax
            << " internalCount=" << internalNorms.count
            << " boundaryCount=" << boundaryNorms.count
            << " cellCount=" << cellNorms.count
            << nl << endl;
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Audit runtime-selected solids4foam stabilisation models on synthetic "
        "displacement fields."
    );

    argList::addOption
    (
        "model",
        "word",
        "Override model.type in system/stabilisationModelAuditDict"
    );

    argList::addOption
    (
        "field",
        "word",
        "Synthetic field: translation, rigidRotation, affine, quadratic, smooth, checkerboard"
    );

    argList::addOption
    (
        "gradient",
        "word",
        "Gradient mode: exact or production"
    );

    argList::addOption
    (
        "scaleFactor",
        "scalar",
        "Override model.scaleFactor"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary auditDict
    (
        IOobject
        (
            "stabilisationModelAuditDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary modelDict(auditDict.subDict("model"));

    word modelType(modelDict.lookup("type"));
    word fieldName(auditDict.lookupOrDefault<word>("field", "translation"));
    word gradientMode(auditDict.lookupOrDefault<word>("gradient", "exact"));
    scalar scaleFactor(modelDict.lookupOrDefault<scalar>("scaleFactor", 1.0));
    const Switch writeFields(auditDict.lookupOrDefault<Switch>("writeFields", false));

    if (optionReadIfPresent(args, "model", modelType))
    {
        modelDict.set("type", modelType);
    }

    optionReadIfPresent(args, "field", fieldName);
    optionReadIfPresent(args, "gradient", gradientMode);

    if (optionReadIfPresent(args, "scaleFactor", scaleFactor))
    {
        modelDict.set("scaleFactor", scaleFactor);
    }

    volVectorField D
    (
        IOobject
        (
            "D_audit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            writeFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero),
        auditPatchTypes
        (
            mesh,
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    volTensorField gradD
    (
        IOobject
        (
            "gradD_audit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            writeFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero),
        auditPatchTypes
        (
            mesh,
            extrapolatedCalculatedFvPatchTensorField::typeName
        )
    );

    setSyntheticFields(mesh, auditDict, fieldName, D, gradD);

    if (gradientMode == "production")
    {
        gradD = fvc::grad(D);
        gradD.correctBoundaryConditions();
    }
    else if (gradientMode != "exact")
    {
        FatalErrorInFunction
            << "Unknown gradient mode " << gradientMode << nl
            << "Valid modes are exact and production" << abort(FatalError);
    }

    autoPtr<stabilisationModel> stab =
        stabilisationModel::New(mesh, modelDict, dimless);

    stab->updateVector(D, &gradD);

    const volVectorField& cellVector = stab->cellVector(nullptr, true);

    report
    (
        modelType,
        fieldName,
        gradientMode,
        stab(),
        stab->faceVector(),
        cellVector
    );

    if (writeFields)
    {
        D.write();
        gradD.write();
        stab->faceVector().write();
        cellVector.write();
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
