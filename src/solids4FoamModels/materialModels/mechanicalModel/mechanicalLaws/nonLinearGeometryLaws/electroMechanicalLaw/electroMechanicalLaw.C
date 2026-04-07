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
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    //activeSigmaPtr_(nullptr)
    
    activeSigma_
    (
     IOobject
     (
      "activeSigma",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
     (
      "zero",
      dimPressure,
      symmTensor::zero
      )
     )
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
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // if (!activeSigmaPtr_.valid())
    //   {
    // 	activeSigmaPtr_.reset
    // 	  (
    // 	   new volSymmTensorField
    // 	   (
    //         IOobject
    //         (
    // 	     "activeSigma",
    // 	     mesh().time().timeName(),
    // 	     mesh(),
    // 	     IOobject::NO_READ,
    // 	     IOobject::AUTO_WRITE
    // 	     ),
    //         mesh(),
    //         dimensionedSymmTensor
    //         (
    // 	     "zero",
    // 	     dimPressure,
    // 	     symmTensor::zero
    // 	     )
    // 	    )
    // 	   );
    //   }
    // Lookup the fibre directions
    const volVectorField& f0 = mesh().lookupObject<volVectorField>("f0");

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();
    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Build the active Cauchy stress from the current unit fibre direction.
    const volVectorField currentFibre(F & f0);
    const volScalarField currentFibreMag
    (
        max
        (
            mag(currentFibre),
            dimensionedScalar("small", dimless, SMALL)
        )
    );
    const volVectorField currentUnitFibre(currentFibre/currentFibreMag);
    activeSigma_ = currentTa*sqr(currentUnitFibre);

    if (mesh().foundObject<volScalarField>("p"))
    {
        activeSigma_ = dev(activeSigma_);
    }

    sigma += activeSigma_;
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Lookup the fibre directions
    const surfaceVectorField f0f
    (
        fvc::interpolate(mesh().lookupObject<volVectorField>("f0"))
    );

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const surfaceTensorField& F = this->Ff();
    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Build the active Cauchy stress from the current unit fibre direction.
    const surfaceVectorField currentFibre(F & f0f);
    const surfaceScalarField currentFibreMag
    (
        max
        (
            mag(currentFibre),
            dimensionedScalar("small", dimless, SMALL)
        )
    );
    const surfaceVectorField currentUnitFibre(currentFibre/currentFibreMag);
    sigma += currentTa*sqr(currentUnitFibre);
}


// ************************************************************************* //
