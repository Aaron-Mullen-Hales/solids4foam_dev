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
    ArosticaActivation

Description
    Electro-coupled mechanical law implementing the time-dependent active
    stress model of Bestel, Clement and Sorine, as used for the Aróstica
    et al. cardiac elastodynamics benchmark ("A software benchmark for
    cardiac elastodynamics", Comput. Methods Appl. Mech. Engrg. 435 (2025)
    117485).

    The stress is composed of a run-time selectable passive hyperelastic
    law plus an active fibre-direction component sigma_active = tau(t) f0 x
    f0, where tau(t) solves the scalar evolution equation (Eq. 5 in the
    reference paper)

        dtau/dt = -|a(t)| tau(t) + sigma0 |a(t)|_+

    with activation function (Eq. 6)

        a(t)   = alphaMax f(t) + alphaMin (1 - f(t))
        f(t)   = Splus(t - tSys) Sminus(t - tDias)
        S+-(x) = 0.5 (1 +- tanh(x/gamma))

    tau(t) is spatially uniform (it depends only on time, not on the
    deformation), so it is integrated once per accepted time step using a
    fixed-substep RK4 scheme fine enough to resolve the gamma time scale,
    independent of the outer solver time step. The result is cached against
    the mesh time index and written out as a volScalarField "Ta" for
    visualisation.

    This law only reproduces the active-stress part of the benchmark
    (Case A / Step 1 activation). It does not implement the companion
    time-dependent "follower" endocardial pressure model (Eqs. 7-8), which
    remains a separate concern.

SourceFiles
    ArosticaActivation.C

\*---------------------------------------------------------------------------*/

#ifndef ArosticaActivation_H
#define ArosticaActivation_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ArosticaActivation Declaration
\*---------------------------------------------------------------------------*/

class ArosticaActivation
:
    public mechanicalLaw
{
    // Private Data

        //- Run-time selectable passive mechanical law
        autoPtr<mechanicalLaw> passiveMechLawPtr_;

        //- Fibre direction field (read from 0/solid/f0)
        volVectorField f0_;

        //- Fibre direction outer product f0 x f0
        volSymmTensorField f0f0_;

        //- Surface fibre direction field (read from 0/solid/f0f)
        surfaceVectorField f0f_;

        //- Surface fibre direction outer product
        surfaceSymmTensorField f0f0f_;

        //- Contractility, sigma0 in Eq. (5)
        const dimensionedScalar sigma0_;

        //- Activation function transition sharpness, gamma in Eq. (6)
        const dimensionedScalar gamma_;

        //- Minimum activation rate, alphaMin in Eq. (6)
        const dimensionedScalar alphaMin_;

        //- Maximum activation rate, alphaMax in Eq. (6)
        const dimensionedScalar alphaMax_;

        //- Systole onset time, tSys in Eq. (6)
        const dimensionedScalar tSys_;

        //- Diastole onset time, tDias in Eq. (6)
        const dimensionedScalar tDias_;

        //- Active tension field, written for post-processing
        mutable volScalarField Ta_;

        //- Mesh time index at which tau_ was last integrated
        mutable label cachedTimeIndex_;

        //- Active tension tau(t) at the cached time index
        mutable scalar tau_;

    // Private Member Functions

        //- Validate the activation model parameters
        void validateActivationParameters() const;

        //- S+-(deltaT) = 0.5*(1 +- tanh(deltaT/gamma_))
        scalar sPlus(const scalar deltaT) const;
        scalar sMinus(const scalar deltaT) const;

        //- Activation function a(t), Eq. (6)
        scalar activationFunction(const scalar t) const;

        //- Time derivative of tau, Eq. (5)
        scalar tauDerivative(const scalar t, const scalar tau) const;

        //- Integrate tau(t) from 0 to tEnd using fixed-substep RK4
        scalar integrateTau(const scalar tEnd) const;

        //- Return tau(t) at the current time, integrating and caching it
        //  once per accepted time step
        scalar currentTau() const;

        //- Disallow default bitwise copy construct
        ArosticaActivation(const ArosticaActivation&);

        //- Disallow default bitwise assignment
        void operator=(const ArosticaActivation&);

public:

    //- Runtime type information
    TypeName("ArosticaActivation");


    // Constructors

        //- Construct from dictionary
        ArosticaActivation
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~ArosticaActivation();


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

        //- True if this law currently provides a non-zero active stress
        bool hasActiveStress() const;

        //- Return the active Cauchy stress
        tmp<volSymmTensorField> activeCauchyStress() const;

        //- Return the active Cauchy stress surface field
        tmp<surfaceSymmTensorField> activeCauchyStressf() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate passive/projected and active/preserved cell stress
        virtual void correctStressComponents
        (
            volSymmTensorField& sigmaToProject,
            volSymmTensorField& sigmaPreserved
        );

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Calculate passive/projected and active/preserved face stress
        virtual void correctStressComponents
        (
            surfaceSymmTensorField& sigmaToProject,
            surfaceSymmTensorField& sigmaPreserved
        );

        //- Forward accumulated-field updates to the passive law
        //- Forward to the passive law: the face deformation gradient that
        //  the face constitutive routines and activeCauchyStressf() read is
        //  the PASSIVE law's Ff(), not this wrapper's own
        virtual void updateFf();

        virtual void updateTotalFields();

        //- Forward restart setup to the passive law
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
