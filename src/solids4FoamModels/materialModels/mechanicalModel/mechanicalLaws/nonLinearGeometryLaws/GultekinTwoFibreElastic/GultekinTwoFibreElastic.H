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
    GultekinTwoFibreElastic

Description
    Aligned two-family exponential anisotropic hyperelastic law created for
    the Gultekin-Dal-Holzapfel artery benchmark.

    The matrix energy is

        Wiso = mu/2 (I1bar - 3)

    and each enabled fibre family contributes

        Wa = k1/(2 k2) (exp(k2 (I4a - 1)^2) - 1).

    By default I4a is the unsplit invariant Ma.C.Ma. The optional split form
    uses I4bar_a = J^(-2/3) I4a and includes the complete derivative of the
    modified invariant. Fibre directions are read directly as reference cell
    and face fields selected by dictionary names.

    Both correct() overloads return passive Cauchy stress. Mixed pressure is
    owned and applied by the calling solid model. The required finite
    bulkModulus controls compressibility through 1/K in the mixed pressure
    equation; it is not inserted as an additional constitutive stress.

    Example:

    \verbatim
        type                 GultekinTwoFibreElastic;
        rho                  rho [1 -3 0 0 0 0 0] 1000;
        mu                   mu  [1 -1 -2 0 0 0 0] 1e4;
        k1                   k1  [1 -1 -2 0 0 0 0] 5e5;
        k2                   k2  [0 0 0 0 0 0 0] 2;
        bulkModulus          K
                             [1 -1 -2 0 0 0 0] 1.0e6;
        fibreField1          f0;
        fibreField2          f1;
        faceFibreField1      f0f;
        faceFibreField2      f1f;
        anisotropicSplit     false;
        fibresTensionOnly    false;
        useSecondFibreFamily true;
        implicitShearModulus implicitShearModulus
                             [1 -1 -2 0 0 0 0] 2.01e6;
        impKcoeff            1.0;
        fibreUnitTolerance   1e-6;
        clipExponent         false;
        exponentLimit        650;
        writeDiagnostics     true;
    \endverbatim

    This is an aligned two-family law. It does not implement fibre dispersion
    and is not a complete dispersed Holzapfel-Gasser-Ogden model.

SourceFiles
    GultekinTwoFibreElastic.C

\*---------------------------------------------------------------------------*/

#ifndef GultekinTwoFibreElastic_H
#define GultekinTwoFibreElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class GultekinTwoFibreElastic Declaration
\*---------------------------------------------------------------------------*/

class GultekinTwoFibreElastic
:
    public mechanicalLaw
{
    // Private data

        //- Initial density
        const dimensionedScalar rho_;

        //- Ground-matrix shear modulus
        dimensionedScalar mu_;

        //- Fibre stress parameter
        dimensionedScalar k1_;

        //- Dimensionless fibre exponential parameter
        dimensionedScalar k2_;

        //- Finite bulk modulus used by the mixed pressure equation
        const dimensionedScalar bulkModulus_;

        //- Use modified fibre invariants
        const Switch anisotropicSplit_;

        //- Suppress compressed fibre response
        const Switch fibresTensionOnly_;

        //- Enable the second fibre family
        const Switch useSecondFibreFamily_;

        //- Cell field name for fibre family 1
        const word fibreField1Name_;

        //- Cell field name for fibre family 2
        const word fibreField2Name_;

        //- Face field name for fibre family 1
        const word faceFibreField1Name_;

        //- Face field name for fibre family 2
        const word faceFibreField2Name_;

        //- Reference cell directions for fibre family 1
        volVectorField fibre1_;

        //- Optional reference cell directions for fibre family 2
        autoPtr<volVectorField> fibre2Ptr_;

        //- Reference face directions for fibre family 1
        surfaceVectorField faceFibre1_;

        //- Optional reference face directions for fibre family 2
        autoPtr<surfaceVectorField> faceFibre2Ptr_;

        //- Scale applied to the numerical implicit stiffness
        const scalar impKcoeff_;

        //- Whether the numerical stiffness was supplied explicitly
        const bool implicitShearModulusSpecified_;

        //- Unscaled numerical stiffness input or default estimate
        dimensionedScalar unscaledImplicitShearModulus_;

        //- Effective constant numerical stiffness used by the mixed solver
        dimensionedScalar implicitShearModulus_;

        //- Unit-vector validation tolerance
        const scalar fibreUnitTolerance_;

        //- Clip the exponential argument instead of terminating
        const Switch clipExponent_;

        //- Maximum allowed exponential argument
        const scalar exponentLimit_;

        //- Create and write optional diagnostic fields
        const Switch writeDiagnostics_;

        //- Jacobian diagnostic
        autoPtr<volScalarField> diagnosticJPtr_;

        //- Fibre-family-1 invariant diagnostic
        autoPtr<volScalarField> diagnosticI4Ptr_;

        //- Fibre-family-2 invariant diagnostic
        autoPtr<volScalarField> diagnosticI6Ptr_;

        //- Fibre-family-1 stretch diagnostic
        autoPtr<volScalarField> diagnosticLambda4Ptr_;

        //- Fibre-family-2 stretch diagnostic
        autoPtr<volScalarField> diagnosticLambda6Ptr_;

        //- Matrix Cauchy stress diagnostic
        autoPtr<volSymmTensorField> diagnosticSigmaIsoPtr_;

        //- Fibre-family-1 Cauchy stress diagnostic
        autoPtr<volSymmTensorField> diagnosticSigmaFibre1Ptr_;

        //- Fibre-family-2 Cauchy stress diagnostic
        autoPtr<volSymmTensorField> diagnosticSigmaFibre2Ptr_;

        //- Complete passive Cauchy stress diagnostic
        autoPtr<volSymmTensorField> diagnosticSigmaPassivePtr_;


    // Private Member Functions

        struct ConstitutiveResult
        {
            scalar J;
            scalar I4;
            scalar I6;
            scalar lambda4;
            scalar lambda6;
            symmTensor sigmaIso;
            symmTensor sigmaFibre1;
            symmTensor sigmaFibre2;
            symmTensor sigmaPassive;
            bool exponentClipped;
            scalar maximumUnclippedExponent;
        };

        //- Evaluate one cell or face using the shared constitutive equations
        ConstitutiveResult evaluateConstitutive
        (
            const tensor& F,
            const vector& fibre1,
            const vector& fibre2
        ) const;

        //- Validate material coefficients and numerical controls
        void validateMaterialParameters() const;

        //- Validate cell and face reference fibre fields
        void validateFibreFields() const;

        //- Allocate optional diagnostic fields
        void makeDiagnostics();

        //- Report exponential clipping statistics
        void reportExponentClipping
        (
            const word& location,
            label clippedCount,
            scalar maximumUnclippedExponent
        ) const;

        //- Disallow default bitwise copy construct
        GultekinTwoFibreElastic(const GultekinTwoFibreElastic&);

        //- Disallow default bitwise assignment
        void operator=(const GultekinTwoFibreElastic&);

public:

    //- Runtime type information
    TypeName("GultekinTwoFibreElastic");


    // Constructors

        //- Construct from dictionary
        GultekinTwoFibreElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~GultekinTwoFibreElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the numerical implicit stiffness
        virtual tmp<volScalarField> impK() const;

        //- Return the finite mixed bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the reference numerical shear stiffness
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate passive cell Cauchy stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate passive face Cauchy stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Split the passive cell stress into the part the calling solid
        //  model may deviatorise and the part that must be preserved
        virtual void correctStressComponents
        (
            volSymmTensorField& sigmaToProject,
            volSymmTensorField& sigmaPreserved
        );

        //- Split the passive face stress into the part the calling solid
        //  model may deviatorise and the part that must be preserved
        virtual void correctStressComponents
        (
            surfaceSymmTensorField& sigmaToProject,
            surfaceSymmTensorField& sigmaPreserved
        );

        //- Set restart field write options
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
