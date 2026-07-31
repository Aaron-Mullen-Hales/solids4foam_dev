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

#include "nonLinGeomTotalLagTotalDispGultekinSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispGultekinSolid, 0);
addToRunTimeSelectionTable
(
    solidModel,
    nonLinGeomTotalLagTotalDispGultekinSolid,
    dictionary
);


tmp<volScalarField>
nonLinGeomTotalLagTotalDispGultekinSolid::mixedVolumetricConstraint
(
    const volScalarField& J
) const
{
    return (J - 1.0)/J;
}


scalar
nonLinGeomTotalLagTotalDispGultekinSolid::mixedVolumetricConstraintValue
(
    const scalar J
) const
{
    return (J - 1.0)/J;
}


scalar
nonLinGeomTotalLagTotalDispGultekinSolid::
mixedVolumetricConstraintDerivative
(
    const scalar J
) const
{
    return 1.0/sqr(J);
}


string
nonLinGeomTotalLagTotalDispGultekinSolid::
mixedVolumetricConstraintDescription() const
{
    return "(J - 1)/J";
}


nonLinGeomTotalLagTotalDispGultekinSolid::
nonLinGeomTotalLagTotalDispGultekinSolid
(
    Time& runTime,
    const word& region
)
:
    nonLinGeomTotalLagTotalDispSolid(typeName, runTime, region),
    writeGultekinMixedDiagnostics_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writeGultekinMixedDiagnostics",
            false
        )
    )
{
    if (!solvePressure())
    {
        FatalErrorInFunction
            << typeName << " requires solvePressure true" << nl
            << exit(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = mechanical();

    if
    (
        laws.size() != 1
     || laws[0].type() != word("GultekinTwoFibreElastic")
    )
    {
        FatalErrorInFunction
            << typeName << " requires exactly one mechanical law of type "
            << "GultekinTwoFibreElastic." << nl
            << "Selected law count = " << laws.size() << nl;

        if (laws.size())
        {
            FatalError
                << "First selected law type = " << laws[0].type() << nl;
        }

        FatalError << exit(FatalError);
    }

    const dictionary& lawDict = laws[0].dict();
    const Switch anisotropicSplit
    (
        lawDict.lookupOrDefault<Switch>("anisotropicSplit", false)
    );
    const Switch fibresTensionOnly
    (
        lawDict.lookupOrDefault<Switch>("fibresTensionOnly", false)
    );
    const Switch useSecondFibreFamily
    (
        lawDict.lookupOrDefault<Switch>("useSecondFibreFamily", true)
    );

    if (anisotropicSplit || fibresTensionOnly || !useSecondFibreFamily)
    {
        FatalErrorInFunction
            << typeName << " targets the original two-family unsplit WAS "
            << "continuum model and therefore requires:" << nl
            << "    anisotropicSplit false;" << nl
            << "    fibresTensionOnly false;" << nl
            << "    useSecondFibreFamily true;" << nl
            << "Selected values are " << anisotropicSplit << ", "
            << fibresTensionOnly << ", " << useSecondFibreFamily
            << exit(FatalError);
    }

    Info<< nl
        << "Gultekin mixed continuum mode active" << nl << nl
        << "    total Cauchy stress:" << nl
        << "        sigma = sigmaPassive - p I" << nl << nl
        << "    pressure constraint:" << nl
        << "        -p/K -(J - 1)/J + S_p = 0" << nl << nl
        << "    volumetric Cauchy stress when S_p = 0:" << nl
        << "        -p I = K*(J - 1)/J I" << nl << nl
        << "    finite-volume model; no Q1P0 finite elements are assembled"
        << nl
        << "    anisotropicSplit = false" << nl
        << "    fibresTensionOnly = false" << nl
        << "    useSecondFibreFamily = true" << nl
        << "    writeGultekinMixedDiagnostics = "
        << writeGultekinMixedDiagnostics_ << endl;
}


void nonLinGeomTotalLagTotalDispGultekinSolid::writeFields
(
    const Time& runTime
)
{
    solidModel::writeFields(runTime);

    if (!writeGultekinMixedDiagnostics_ || !runTime.writeTime())
    {
        return;
    }

    // Refresh all three terms using the same methods as the production
    // pressure residual. These fields are diagnostic only.
    refreshMixedPressureConstraintDiagnostics();

    const volSymmTensorField sigmaPassiveFull
    (
        IOobject
        (
            "sigmaPassiveFull",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sigma() + p()*I
    );

    const volScalarField fibreMeanStress
    (
        IOobject
        (
            "fibreMeanStress",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tr(sigmaPassiveFull)/3.0
    );

    const volSymmTensorField sigmaVolumetricFromP
    (
        IOobject
        (
            "sigmaVolumetricFromP",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -p()*dimensionedSymmTensor("I", dimless, I)
    );

    const volScalarField gultekinConstraintPhysical
    (
        IOobject
        (
            "gultekinConstraintPhysical",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureCompressibilityDiagnostic()
      + volumetricConstraintDiagnostic()
    );

    const volScalarField pressureStabilisationContribution
    (
        IOobject
        (
            "pressureStabilisationContribution",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pressureStabilisationDiagnostic()
    );

    const volScalarField completePressureResidual
    (
        IOobject
        (
            "completePressureResidual",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        completePressureResidualDiagnostic()
    );

    sigmaPassiveFull.write();
    fibreMeanStress.write();
    sigmaVolumetricFromP.write();
    gultekinConstraintPhysical.write();
    pressureStabilisationContribution.write();
    completePressureResidual.write();
}

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
