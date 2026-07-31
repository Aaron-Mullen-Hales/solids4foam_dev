/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    ArosticaHolzapfelOgdenViscoelastic

Description
    Aróstica Holzapfel-Ogden passive law with exact material viscosity.

    The viscous contribution is Sviscous = eta*Edot, with
    E = 0.5*(F.T()*F - I), and is returned as the complete Cauchy stress
    symm(F&Sviscous&F.T())/J.  The law contains no pressure or volumetric
    penalty stress.

SourceFiles
    ArosticaHolzapfelOgdenViscoelastic.C

\*---------------------------------------------------------------------------*/

#ifndef ArosticaHolzapfelOgdenViscoelastic_H
#define ArosticaHolzapfelOgdenViscoelastic_H

#include "ArosticaHolzapfelOgdenElastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class ArosticaHolzapfelOgdenViscoelastic Declaration
\*---------------------------------------------------------------------------*/

class ArosticaHolzapfelOgdenViscoelastic
:
    public ArosticaHolzapfelOgdenElastic
{
    // Private data

        //- Material viscosity, with dimensions pressure*time
        const dimensionedScalar eta_;

        //- Optional diagnostic output switch
        const Switch writeDiagnostics_;

        //- Configured temporal derivative scheme
        const word ddtScheme_;

        //- Accepted current and previous strain histories
        volSymmTensorField E_;
        volSymmTensorField EOld_;
        volSymmTensorField EOldOld_;
        surfaceSymmTensorField Ef_;
        surfaceSymmTensorField EfOld_;
        surfaceSymmTensorField EfOldOld_;

        //- Trial derivative and viscous stress diagnostics
        volSymmTensorField Edot_;
        volSymmTensorField Sviscous_;
        volSymmTensorField sigmaViscous_;
        surfaceSymmTensorField Edotf_;
        surfaceSymmTensorField Sviscousf_;
        surfaceSymmTensorField sigmaViscousf_;

        //- Explicit histories are available after restart or first commit
        bool haveHistory_;
        bool haveSecondHistory_;

        //- Time index at which accepted histories were most recently shifted
        label lastHistoryUpdateTimeIndex_;


    // Private Member Functions

        //- Return a deterministic field name for this law instance
        word fieldName(const word& prefix) const;

        //- Check whether a typed field exists in the current time directory
        template<class FieldType>
        bool fieldPresent(const word& name) const;

        //- Read and validate the configured ddt scheme
        word readDdtScheme() const;

        //- Return the temporal coefficients used by OpenFOAM v2312 backward
        void temporalCoefficients
        (
            scalar& current,
            scalar& old,
            scalar& oldOld
        ) const;

        //- Calculate Green-Lagrange strain from a deformation gradient
        static symmTensor strain(const tensor& F);

        //- Seed accepted history only during initialisation/restart recovery
        void initialiseHistoryIfRequired();

        //- Refresh cell and face viscous fields for the current trial F/Ff
        void updateViscousFields();

        //- Return viscous stress for a tensor-level tangent evaluation
        symmTensor viscousStress
        (
            const tensor& F,
            const symmTensor& EOld,
            const symmTensor& EOldOld,
            const scalar current,
            const scalar old,
            const scalar oldOld
        ) const;

        //- Fill a numerical tangent containing passive and viscous stress
        void fillViscoelasticTangent
        (
            mat66& tangent,
            const tensor& F,
            const vector& f0,
            const vector& s0,
            const vector& n0,
            const symmTensor& EOld,
            const symmTensor& EOldOld,
            const scalar current,
            const scalar old,
            const scalar oldOld
        ) const;

        //- Disallow default bitwise copy construct
        ArosticaHolzapfelOgdenViscoelastic
        (
            const ArosticaHolzapfelOgdenViscoelastic&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const ArosticaHolzapfelOgdenViscoelastic&
        );


public:

    //- Runtime type information
    TypeName("ArosticaHolzapfelOgdenViscoelastic");


    // Constructors

        //- Construct from dictionary
        ArosticaHolzapfelOgdenViscoelastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~ArosticaHolzapfelOgdenViscoelastic();


    // Member Functions

        //- Calculate cell stress including material viscosity
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate direct face stress including material viscosity
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Populate the current-trial numerical tangent
        virtual void materialTangentField(List<mat66>& matTan) const;

        //- Populate d(F*Sviscous)/dF for each direct face. Nine output
        //  tensors are returned per face, one for each flattened F entry.
        virtual bool viscousNominalTangentField
        (
            List<tensor>& tangent
        ) const;

        //- Advance accepted strain histories once per accepted timestep
        virtual void updateTotalFields();

        //- Set restart write options for all required history fields
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
