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

#include "rotatedSpatialDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatedSpatialDisplacementFvPatchVectorField::rotatedSpatialDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    thetha_(0.0),
    axis_(vector(0, 0, 1)),
    epsilon_(0.0)
{}


rotatedSpatialDisplacementFvPatchVectorField::rotatedSpatialDisplacementFvPatchVectorField
(
    const rotatedSpatialDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    thetha_(ptf.thetha_),
    axis_(ptf.axis_),
    epsilon_(ptf.epsilon_)
{}


rotatedSpatialDisplacementFvPatchVectorField::rotatedSpatialDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    thetha_(readScalar(dict.lookup("thetha"))),
    axis_(dict.lookup("axis")),
    epsilon_(readScalar(dict.lookup("epsilon")))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}

#ifndef OPENFOAM_ORG
rotatedSpatialDisplacementFvPatchVectorField::rotatedSpatialDisplacementFvPatchVectorField
(
    const rotatedSpatialDisplacementFvPatchVectorField& ptf
)
:
    fixedDisplacementFvPatchVectorField(ptf),
    thetha_(ptf.thetha_),
    axis_(ptf.axis_),
    epsilon_(ptf.epsilon_)
{}
#endif

rotatedSpatialDisplacementFvPatchVectorField::rotatedSpatialDisplacementFvPatchVectorField
(
    const rotatedSpatialDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(ptf, iF),
    thetha_(ptf.thetha_),
    axis_(ptf.axis_),
    epsilon_(ptf.epsilon_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void rotatedSpatialDisplacementFvPatchVectorField::updateCoeffs()
  {
    if (this->updated())
      {
        return;
      }

    scalar rad = thetha_ * constant::mathematical::pi / 180.0;
    scalar epsilon = epsilon_;
    scalar c = cos(rad);
    scalar s = sin(rad);

    vector k = axis_;
    scalar magnitude = mag(k);

    if (magnitude == 0)
      {
        FatalErrorInFunction << "Rotation axis has zero magnitude" << exit(FatalError);
      }

    k /= magnitude;


    const fvPatch& patch = this->patch();

    // Directly get face centres on this patch
    const vectorField& faceCentres = patch.patch().faceCentres();

    vectorField newDisp(patch.size(), vector::zero);

    forAll(newDisp, i)
      {
	vector X = faceCentres[i];
	//the perturbation below causes memory issues because it is using original reference rather than a copy of it
	//X.x() = X.x() + epsilon;

	// Rodrigues' rotation formula
	vector v_rot =
	  X * c +
	  (k ^ X) * s +
	  k * (k & X) * (1 - c);

	v_rot.x() *= (1.0 + epsilon);
	newDisp[i] = v_rot - X;
      }
    

    totalDisp() = newDisp;

    fixedDisplacementFvPatchVectorField::updateCoeffs();
  }



void rotatedSpatialDisplacementFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("thetha")
        << thetha_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("epsilon")
        << epsilon_ << token::END_STATEMENT << nl;

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    rotatedSpatialDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
