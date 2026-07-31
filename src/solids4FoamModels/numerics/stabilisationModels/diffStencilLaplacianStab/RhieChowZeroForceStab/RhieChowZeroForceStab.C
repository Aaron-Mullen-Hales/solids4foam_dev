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

#include "RhieChowZeroForceStab.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RhieChowZeroForceStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, RhieChowZeroForceStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::RhieChowZeroForceStab::RhieChowZeroForceStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    RhieChowStab(mesh, dict, dims),
    reportZeroForceDiagnostics_
    (
        dict.lookupOrDefault<Switch>("reportZeroForceDiagnostics", false)
    ),
    zeroForceAbsTolerance_
    (
        dict.lookupOrDefault<scalar>("zeroForceAbsTolerance", 1e-12)
    ),
    zeroForceRelTolerance_
    (
        dict.lookupOrDefault<scalar>("zeroForceRelTolerance", 1e-10)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RhieChowZeroForceStab::projectExtensiveVectorForce
(
    vectorField& stabilisationForce
) const
{
    const vector netForceBefore = gSum(stabilisationForce);
    const scalar sumMagBefore = gSum(mag(stabilisationForce));
    const scalar totalVolume = gSum(mesh().V());

    if (totalVolume <= VSMALL)
    {
        FatalErrorInFunction
            << "Cannot apply zero-force projection: total mesh volume is "
            << totalVolume << abort(FatalError);
    }

    const vector correctionDensity = netForceBefore/totalVolume;
    const scalarField& V = mesh().V();

    forAll(stabilisationForce, cellI)
    {
        stabilisationForce[cellI] -= V[cellI]*correctionDensity;
    }

    const vector netForceAfter = gSum(stabilisationForce);
    const scalar sumMagAfter = gSum(mag(stabilisationForce));

    const scalar ratioBefore =
        mag(netForceBefore)/(sumMagBefore + VSMALL);
    const scalar ratioAfter =
        mag(netForceAfter)/(sumMagAfter + VSMALL);

    if (reportZeroForceDiagnostics_)
    {
        Info<< "RhieChowZeroForceStab diagnostics:" << nl
            << "    net stabilisation force before projection = "
            << netForceBefore << nl
            << "    net stabilisation force after projection = "
            << netForceAfter << nl
            << "    sum of local stabilisation-force magnitudes before "
            << "projection = " << sumMagBefore << nl
            << "    sum of local stabilisation-force magnitudes after "
            << "projection = " << sumMagAfter << nl
            << "    ratio of net force to sum of magnitudes before "
            << "projection = " << ratioBefore << nl
            << "    ratio of net force to sum of magnitudes after "
            << "projection = " << ratioAfter << nl
            << "    correction force density = " << correctionDensity
            << endl;

        if (ratioBefore <= 100.0*SMALL)
        {
            Info<< "Existing Rhie-Chow term is already globally "
                << "self-equilibrated;" << nl
                << "zero-force projection is expected to have negligible "
                << "effect." << endl;
        }
    }

    if
    (
        mag(netForceAfter) > zeroForceAbsTolerance_
     && ratioAfter > zeroForceRelTolerance_
    )
    {
        FatalErrorInFunction
            << "Zero-force projection failed: net stabilisation force after "
            << "projection = " << netForceAfter
            << ", magnitude = " << mag(netForceAfter)
            << ", relative ratio = " << ratioAfter
            << ", absolute tolerance = " << zeroForceAbsTolerance_
            << ", relative tolerance = " << zeroForceRelTolerance_
            << abort(FatalError);
    }
}


// ************************************************************************* //
