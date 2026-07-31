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

Class
    GuccioneElastic

Description
    Guccione hyper-elasticity law suitable for myocardium. This implementation
    is based on the approach described in E Garcia-Blanco, R Ortigosa, AJ Gil,
    CH Lee, J Bonet, A new computational framework for electro-activation in
    cardiac mechanics, Comput. Methods Appl. Mech. Engrg. 348 (2019) 796–845.

SourceFiles
    GuccioneElastic.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef GuccioneElastic_H
#define GuccioneElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class GuccioneElastic Declaration
\*---------------------------------------------------------------------------*/

class GuccioneElastic
:
    public mechanicalLaw
{
    // Private data

        //- Bulk modulus
        dimensionedScalar bulkModulus_;

        //- Guccione modulus parameter
        dimensionedScalar k_;

        //- Material parameters
        const scalar cf_;
        const scalar ct_;
        const scalar cfs_;

        //- Shear modulus in the limited of small strains
        dimensionedScalar mu_;

        //- Pressure-displacement formulation
        const Switch pressureDisplacement_;

        //- Evaluate passive Guccione law using Fbar = J^(-1/3) F
        const Switch useIsochoricSplit_;

        //- Enable optional apex regularisation experiment
        const Switch apexRegularisationEnabled_;

        //- Apex cellZone name
        const word apexRegularisationZoneName_;

        //- Use isotropic Guccione coefficients in the selected apex zone
        const Switch apexIsotropisePassive_;

        //- Isotropic Guccione b coefficient for the apex zone
        const scalar apexIsotropicB_;

        //- Preserve baseline implicit stiffness and momentum-stab coefficient
        const Switch preserveBaselineImpK_;

        //- Write apex diagnostic fields
        const Switch writeApexDiagnosticFields_;

        //- Internal apex mask used by Guccione coefficient selection
        mutable autoPtr<volScalarField> apexRegularisationMaskPtr_;

        //- Interpolated apex mask for surface helper paths
        mutable autoPtr<surfaceScalarField> apexRegularisationMaskfPtr_;

        //- Passive isotropic diagnostic mask
        mutable autoPtr<volScalarField> apexPassiveIsotropicMaskPtr_;

        //- Baseline implicit-stiffness preservation diagnostic mask
        mutable autoPtr<volScalarField> apexBaselineImpKMaskPtr_;

        //- Effective shear modulus for coupledPressureDisplacementSolid
        volScalarField muEff_;

        //- Flag to indicate a uniform fibre field
        //  If true, the fibre direction is read from mechanicalProperties
        //  If false, the fibre field is read from the starting time directory
        const Switch uniformFibreField_;

        //- Fibre directions in the reference configuration
        volVectorField f0_;
        surfaceVectorField f0f_;

        //- Sheet directions in the reference configuration
        //  f0 & s0 == 0
        volVectorField s0_;
        surfaceVectorField s0f_;

        //- Sheet normal directions in the reference configuration
        //  n0 & f0 == n0 & s0 == 0
        volVectorField n0_;
        surfaceVectorField n0f_;

        //- Rotation matrix from local fibre coordinate system to the global
        //  coordinate system
        volTensorField R_;
        surfaceTensorField Rf_;

        //- Outer product of the fibre directions with themselves
        volSymmTensorField f0f0_;
        surfaceSymmTensorField f0f0f_;

        //- 2nd Piola-Kirchhoff stress (without the hydrostatic term)
        volSymmTensorField S_;
        surfaceSymmTensorField Sf_;

        //- Exponential of Qf for coupledPressureDisplacementSolid
        surfaceScalarField expQf_;

        //- Scaling coefficient for coupledPressureDisplacementSolid impK
        const scalar impKcoeff_;


    // Private Member Functions

        //- Create f0 field
        tmp<volVectorField> makeF0
        (
            const Switch& uniformFibreField,
            const fvMesh& mesh,
            const dictionary& dict
        ) const;

        //- Create f0f field
        tmp<surfaceVectorField> makeF0f
        (
            const Switch& uniformFibreField,
            const fvMesh& mesh,
            const dictionary& dict,
            const volVectorField& f0
        ) const;

        //- Calculate stress given gradD and pressure
        void calculateStress
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& gradD
        );

        //- True if passive stress coefficients should be changed
        bool useApexPassiveIsotropisation() const;

        //- True if implicit/numerical stiffness should use apex coefficients
        bool useApexImplicitStiffness() const;

        //- Create apex regularisation masks
        void makeApexRegularisationMasks() const;

        //- Return the apex regularisation mask
        const volScalarField& apexRegularisationMask() const;

        //- Return the apex regularisation surface mask
        const surfaceScalarField& apexRegularisationMaskf() const;

        //- Write apex diagnostic fields
        void writeApexDiagnosticFields() const;

        //- Report physical and implicit apex coefficients
        void reportApexRegularisationCoefficients() const;

        //- Validate assembled tangent columns for apex debug tests
        void validateMaterialTangentColumns
        (
            const List<mat66>& matTan,
            const word& context
        ) const;

        //- Optional debug validation of apex tangent column assembly
        void validateApexRegularisationTangent() const;

        //- Return a cell Guccione coefficient for implicit stiffness paths
        scalar implicitCoefficient
        (
            const scalar baseline,
            const label cellI
        ) const;

        //- Calculate deviatoric Cauchy stress from a deformation gradient
        void calcDevCauchy
        (
            const tensor& F,
            const symmTensor& f0f0,
            const tensor& R,
            symmTensor& devCauchy
        ) const;

        //- Calculate deviatoric Cauchy stress with explicit switches
        void calcDevCauchy
        (
            const tensor& F,
            const symmTensor& f0f0,
            const tensor& R,
            const Switch& useIsochoricSplit,
            const Switch& useLocalCoordSys,
            const scalar cf,
            const scalar ct,
            const scalar cfs,
            symmTensor& devCauchy
        ) const;

        //- Optional debug validation of the isochoric split algebra
        void validateIsochoricSplitImplementation() const;

        //- Calculate initial shear modulus by F perturbation
        void calcInitialShearModulus();

        //- Calculate effective shear modulus
        void calcEffectiveShearModulus();

        //- Disallow default bitwise copy construct
        GuccioneElastic(const GuccioneElastic&);

        //- Disallow default bitwise assignment
        void operator=(const GuccioneElastic&);

public:

    //- Runtime type information
    TypeName("GuccioneElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        GuccioneElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~GuccioneElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Material tangent for Newton Raphson loop
        virtual void materialTangentField(List<mat66>& matTan) const;

        //- Return the bulk modulus field
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the shear modulus field
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();

        //- Update total accumulated fields
        virtual void updateTotalFields();

        //- Return the cell fibre dyad used by this law
        const volSymmTensorField& fibreDyad() const
        {
            return f0f0_;
        }

        //- Return the face fibre dyad used by this law
        const surfaceSymmTensorField& faceFibreDyad() const
        {
            return f0f0f_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
