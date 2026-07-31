/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    ArosticaHolzapfelOgdenElastic

Description
    Purely elastic non-volumetric Holzapfel-Ogden law for the Aróstica
    cardiac elastodynamics benchmark.

    The constitutive law reads the supplied f0, s0, n0 and face f0f, s0f, n0f
    fields directly. It returns the complete non-volumetric passive Cauchy
    stress. Mixed pressure is owned by the calling solid model.

SourceFiles
    ArosticaHolzapfelOgdenElastic.C

\*---------------------------------------------------------------------------*/

#ifndef ArosticaHolzapfelOgdenElastic_H
#define ArosticaHolzapfelOgdenElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ArosticaHolzapfelOgdenElastic Declaration
\*---------------------------------------------------------------------------*/

class ArosticaHolzapfelOgdenElastic
:
    public mechanicalLaw
{
protected:

    // Protected data

        //- Density
        const dimensionedScalar rho_;

        //- Passive exponential coefficients
        const dimensionedScalar a_;
        const dimensionedScalar b_;
        const dimensionedScalar af_;
        const dimensionedScalar bf_;
        const dimensionedScalar as_;
        const dimensionedScalar bs_;
        const dimensionedScalar afs_;
        const dimensionedScalar bfs_;

        //- Finite bulk modulus used by the mixed pressure equation
        const dimensionedScalar bulkModulus_;

        //- Scaling coefficient for numerical implicit stiffness
        const scalar impKcoeff_;

        //- Compression-switch name and steepness
        const word compressionSwitch_;
        const scalar compressionSwitchK_;

        //- Direction validation and tangent perturbation controls
        const scalar directionTolerance_;
        const scalar tangentEps_;

        //- Cell field names
        const word fibreFieldName_;
        const word sheetFieldName_;
        const word sheetNormalFieldName_;

        //- Face field names
        const word faceFibreFieldName_;
        const word faceSheetFieldName_;
        const word faceSheetNormalFieldName_;

        //- Supplied cell reference triad
        volVectorField f0_;
        volVectorField s0_;
        volVectorField n0_;

        //- Supplied face reference triad
        surfaceVectorField f0f_;
        surfaceVectorField s0f_;
        surfaceVectorField n0f_;

        //- Positive representative small-strain shear modulus
        dimensionedScalar muRepresentative_;


    // Protected Member Functions

        struct CompressionSwitchResult
        {
            scalar chi;
            scalar dChi;
        };

        struct ConstitutiveResult
        {
            scalar J;
            scalar energy;
            scalar I1bar;
            scalar I4f;
            scalar I4s;
            scalar I8fs;
            symmTensor S;
            symmTensor sigma;
        };

        //- Evaluate chi and its derivative for one fibre/sheet invariant
        CompressionSwitchResult compressionSwitchValue
        (
            const scalar I4
        ) const;

        //- Shared tensor-level constitutive kernel
        ConstitutiveResult evaluateConstitutive
        (
            const tensor& F,
            const vector& f0,
            const vector& s0,
            const vector& n0
        ) const;

        //- Validate material dimensions and scalar values
        void validateMaterialParameters() const;

        //- Validate supplied cell and face triads
        void validateReferenceFields() const;

        //- Calculate the representative local shear modulus
        void calculateRepresentativeShearModulus();

        //- Fill one numerical constitutive tangent
        void fillTangent
        (
            mat66& tangent,
            const tensor& F,
            const vector& f0,
            const vector& s0,
            const vector& n0
        ) const;

        //- Disallow default bitwise copy construct
        ArosticaHolzapfelOgdenElastic
        (
            const ArosticaHolzapfelOgdenElastic&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const ArosticaHolzapfelOgdenElastic&
        );


public:

    //- Runtime type information
    TypeName("ArosticaHolzapfelOgdenElastic");


    // Constructors

        //- Construct from dictionary
        ArosticaHolzapfelOgdenElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~ArosticaHolzapfelOgdenElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the numerical implicit stiffness
        virtual tmp<volScalarField> impK() const;

        //- Return the finite mixed bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the representative shear modulus
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate the complete passive cell Cauchy stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the complete passive face Cauchy stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Populate the numerical constitutive tangent
        virtual void materialTangentField(List<mat66>& matTan) const;

        //- Populate the passive numerical constitutive tangent explicitly
        virtual void passiveMaterialTangentField(List<mat66>& matTan) const;

        //- Populate the direct-face passive nominal-force tangent
        virtual bool passiveNominalTangentField(List<tensor>& tangent) const;

        //- Set restart field write options
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
