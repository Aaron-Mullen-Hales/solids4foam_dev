/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solidModel.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    autoPtr<solidModel> solid =
        solidModel::New(runTime, dynamicFvMesh::defaultRegion);
    volVectorField& D = solid->D();
    D.correctBoundaryConditions();

    const label nSamples = 20;
    scalar totalCpuTime = 0.0;

    for (label sampleI = 0; sampleI < nSamples; ++sampleI)
    {
        const scalar start = runTime.elapsedCpuTime();
        const tmp<volVectorField> tDdot0(fvc::ddt(D));
        const tmp<volVectorField> tDdot1(fvc::ddt(D));
        totalCpuTime += runTime.elapsedCpuTime() - start;
    }

    Info<< "PASS: fvc::ddt(D) cost probe" << nl
        << "    full-volume constructions per sample = 2" << nl
        << "    samples = " << nSamples << nl
        << "    total construction CPU time = " << totalCpuTime << " s" << nl
        << "    mean pair CPU time = " << totalCpuTime/nSamples << " s"
        << endl;

    return 0;
}

// ************************************************************************* //
