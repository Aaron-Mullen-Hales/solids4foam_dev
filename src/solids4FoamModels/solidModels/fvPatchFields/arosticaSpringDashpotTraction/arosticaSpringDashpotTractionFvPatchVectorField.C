/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "arosticaSpringDashpotTractionFvPatchVectorField.H"
#include "fvc.H"
#include "IStringStream.H"

#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

dimensionedScalar
arosticaSpringDashpotTractionFvPatchVectorField::readCoefficient
(
    const dictionary& dict,
    const word& dimensionedName,
    const word& legacyName,
    const dimensionSet& expectedDimensions
)
{
    dimensionedScalar coefficient
    (
        dimensionedName,
        expectedDimensions,
        0.0
    );

    if (dict.found(dimensionedName))
    {
        const dimensionedScalar input(dict.lookup(dimensionedName));

        // Check before assigning into the dimensioned result. OpenFOAM's
        // dimensioned assignment otherwise aborts generically before the
        // boundary condition can report the offending dictionary entry.
        if (input.dimensions() != expectedDimensions)
        {
            FatalIOErrorInFunction(dict)
                << dimensionedName << " has dimensions "
                << input.dimensions() << "; expected "
                << expectedDimensions << exit(FatalIOError);
        }

        coefficient = dimensionedScalar
        (
            dimensionedName,
            expectedDimensions,
            input.value()
        );
    }
    else if (dict.found(legacyName))
    {
        coefficient.value() = readScalar(dict.lookup(legacyName));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry " << dimensionedName
            << " (legacy alias " << legacyName << ')'
            << exit(FatalIOError);
    }

    if (coefficient.dimensions() != expectedDimensions)
    {
        FatalIOErrorInFunction(dict)
            << dimensionedName << " has dimensions "
            << coefficient.dimensions() << "; expected "
            << expectedDimensions << exit(FatalIOError);
    }

    if
    (
        !std::isfinite(coefficient.value())
     || coefficient.value() < 0.0
    )
    {
        FatalIOErrorInFunction(dict)
            << dimensionedName << " must be finite and non-negative; value = "
            << coefficient.value() << exit(FatalIOError);
    }

    return coefficient;
}


Switch
arosticaSpringDashpotTractionFvPatchVectorField::requireReferenceArea
(
    const dictionary& dict
)
{
    if
    (
        dict.found("useUndeformedArea")
     && !dict.lookupOrDefault<Switch>("useUndeformedArea", false)
    )
    {
        FatalIOErrorInFunction(dict)
            << "Aróstica spring-dashpot traction requires "
            << "useUndeformedArea true"
            << exit(FatalIOError);
    }

    return true;
}


dictionary
arosticaSpringDashpotTractionFvPatchVectorField::withTractionDefaults
(
    const dictionary& input
)
{
    dictionary result(input);

    if (!result.found("traction"))
    {
        IStringStream stream("traction uniform (0 0 0);");
        entry::New(result, stream);
    }

    if (!result.found("pressure"))
    {
        IStringStream stream("pressure uniform 0;");
        entry::New(result, stream);
    }

    return result;
}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    springCoefficient_
    (
        "springCoefficient",
        dimPressure/dimLength,
        0.0
    ),
    dashpotCoefficient_
    (
        "dashpotCoefficient",
        dimPressure*dimTime/dimLength,
        0.0
    ),
    forceReferenceArea_(true),
    writeDiagnostics_(false)
{}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, withTractionDefaults(dict)),
    springCoefficient_
    (
        readCoefficient
        (
            dict,
            "springCoefficient",
            "alpha",
            dimPressure/dimLength
        )
    ),
    dashpotCoefficient_
    (
        readCoefficient
        (
            dict,
            "dashpotCoefficient",
            "beta",
            dimPressure*dimTime/dimLength
        )
    ),
    forceReferenceArea_(requireReferenceArea(dict)),
    writeDiagnostics_
    (
        dict.lookupOrDefault<Switch>("writeDiagnostics", false)
    )
{
    Info<< "Creating " << type() << " boundary condition" << nl
        << "    springCoefficient = " << springCoefficient_ << nl
        << "    dashpotCoefficient = " << dashpotCoefficient_ << nl
        << "    useUndeformedArea = " << forceReferenceArea_ << endl;
}


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(pvf, p, iF, mapper),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    forceReferenceArea_(pvf.forceReferenceArea_),
    writeDiagnostics_(pvf.writeDiagnostics_)
{}


#ifndef OPENFOAM_ORG
arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf
)
:
    solidTractionFvPatchVectorField(pvf),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    forceReferenceArea_(pvf.forceReferenceArea_),
    writeDiagnostics_(pvf.writeDiagnostics_)
{}
#endif


arosticaSpringDashpotTractionFvPatchVectorField::
arosticaSpringDashpotTractionFvPatchVectorField
(
    const arosticaSpringDashpotTractionFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(pvf, iF),
    springCoefficient_(pvf.springCoefficient_),
    dashpotCoefficient_(pvf.dashpotCoefficient_),
    forceReferenceArea_(pvf.forceReferenceArea_),
    writeDiagnostics_(pvf.writeDiagnostics_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void arosticaSpringDashpotTractionFvPatchVectorField::reportDiagnostics
(
    const vectorField& displacement,
    const vectorField& velocity,
    const vectorField& spring,
    const vectorField& dashpot,
    const vectorField& total
) const
{
    if (!writeDiagnostics_)
    {
        return;
    }

    scalar maxDisplacement = 0.0;
    scalar maxVelocity = 0.0;
    scalar maxSpring = 0.0;
    scalar maxDashpot = 0.0;
    scalar maxTraction = 0.0;

    forAll(displacement, faceI)
    {
        maxDisplacement = max(maxDisplacement, mag(displacement[faceI]));
        maxVelocity = max(maxVelocity, mag(velocity[faceI]));
        maxSpring = max(maxSpring, mag(spring[faceI]));
        maxDashpot = max(maxDashpot, mag(dashpot[faceI]));
        maxTraction = max(maxTraction, mag(total[faceI]));
    }

    Info<< "Aróstica spring-dashpot diagnostics, patch "
        << patch().name() << ':' << nl
        << "    max |D| = " << maxDisplacement << nl
        << "    max |Ddot| = " << maxVelocity << nl
        << "    max |spring traction| = " << maxSpring << nl
        << "    max |dashpot traction| = " << maxDashpot << nl
        << "    max |traction| = " << maxTraction << endl;
}


void arosticaSpringDashpotTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& D =
        db().lookupObject<volVectorField>("D");

    const label patchI = patch().index();

    // solidTraction::evaluate() calls updateCoeffs() before fixedGradient
    // reconstructs the patch value.  The stored patch value can therefore
    // belong to the preceding PETSc trial.  Use the current trial cell value
    // as the stateless patch predictor before evaluating ddt(D) and the
    // support laws.  The base class subsequently reconstructs the final
    // fixed-gradient patch value from the newly calculated traction.
    Field<vector>::operator=(patchInternalField());

    const vectorField displacement
    (
        D.boundaryField()[patchI]
    );

    // fvc::ddt uses the same configured ddt(D) scheme as the volume
    // momentum equation. It reads current trial D and accepted old-time
    // fields without committing any state.
    const tmp<volVectorField> tDdot(fvc::ddt(D));
    const vectorField velocity
    (
        tDdot().boundaryField()[patchI]
    );

    vectorField spring(patch().size(), vector::zero);
    vectorField dashpot(patch().size(), vector::zero);

    calculateContributions(spring, dashpot, displacement, velocity);
    traction() = spring + dashpot;
    pressure() = scalar(0.0);

    reportDiagnostics
    (
        displacement,
        velocity,
        spring,
        dashpot,
        traction()
    );

    solidTractionFvPatchVectorField::updateCoeffs();
}


void arosticaSpringDashpotTractionFvPatchVectorField::write
(
    Ostream& os
) const
{
    solidTractionFvPatchVectorField::write(os);

    // dimensionedScalar::writeEntry omits the dimensioned value name when it
    // matches the dictionary keyword.  The constructor reads the entry as a
    // complete dimensionedScalar and therefore requires that name token to be
    // present in a written restart dictionary.
    dimensionedScalar springCoefficientForWrite
    (
        "springCoefficientValue",
        springCoefficient_.dimensions(),
        springCoefficient_.value()
    );
    dimensionedScalar dashpotCoefficientForWrite
    (
        "dashpotCoefficientValue",
        dashpotCoefficient_.dimensions(),
        dashpotCoefficient_.value()
    );

    springCoefficientForWrite.writeEntry("springCoefficient", os);
    dashpotCoefficientForWrite.writeEntry("dashpotCoefficient", os);
    os.writeKeyword("writeDiagnostics")
        << writeDiagnostics_ << token::END_STATEMENT << nl;
}

} // End namespace Foam

// ************************************************************************* //
