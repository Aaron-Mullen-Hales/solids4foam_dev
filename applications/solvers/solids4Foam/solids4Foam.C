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

Application
    solids4Foam

Description
    General solver where the solved mathematical model (fluid, solid or
    fluid-solid) is chosen at run-time.

Author
    Philip Cardiff, UCD.
    Zeljko Tukovic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "solids4FoamWriteHeader.H"

    // Create the general physics class
    autoPtr<physicsModel> physics = physicsModel::New(runTime);

    // Set when an attempted time step could not be solved and must not be
    // accepted. See physicsModel::solutionFailed()
    bool solutionFailed = false;

    while (runTime.run())
    {
        // Update deltaT, if desired, before moving to the next step
        physics().setDeltaT(runTime);

        runTime++;

        if (physics().printInfo())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        // Solve the mathematical model
        // Note: the bool returned by evolve() is deliberately not used as the
        // success test. Its meaning is model-specific (e.g.
        // unsNonLinGeomTotalLagSolid returns false when enforced linearity is
        // engaged, which is not a failure) and it has always been discarded
        // here and by every fluid-solid interface. solutionFailed() is the
        // unambiguous query and defaults to false, so models that never report
        // failure are unaffected
        physics().evolve();

        if (physics().solutionFailed())
        {
            // The attempted time step could not be solved. It must not be
            // accepted: do not advance the accumulated history, do not write
            // it as a result, and do not go on to apply later loads as though
            // this interval had been solved. runTime has already been
            // incremented, so this time is an *attempted* time, not an
            // accepted solution; no attempt is made to rewind the Time object.
            //
            // There is no time-step retry or adaptive-deltaT mechanism in
            // solids4foam, so the safest minimal policy is to stop here rather
            // than silently skip the interval.
            solutionFailed = true;

            Info<< nl
                << "--------------------------------------------------" << nl
                << "SOLUTION FAILED at time = " << runTime.timeName() << nl
                << "--------------------------------------------------" << nl
                << "The physics model reported that this time step could not"
                << " be solved." << nl
                << "It has NOT been accepted: the accumulated fields were not"
                << " advanced and" << nl
                << "no output was written for it. The run stops here rather"
                << " than continuing" << nl
                << "to later loads as though this interval had converged."
                << nl << endl;

            break;
        }

        // Let the physics model know the end of the time-step has been reached
        physics().updateTotalFields();

        if (runTime.outputTime())
        {
            physics().writeFields(runTime);
        }

        if (physics().printInfo())
        {
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }

    physics().end();

    if (solutionFailed)
    {
        // A truthful non-zero exit status. This is a normal return, not an
        // OpenFOAM FatalError/abort: the failure is a physical/numerical
        // outcome that the user asked to be non-fatal, but the run did not
        // complete and must not report success to a calling script
        Info<< nl << "End (solution failed; simulation incomplete)" << nl
            << endl;

        return 1;
    }

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
