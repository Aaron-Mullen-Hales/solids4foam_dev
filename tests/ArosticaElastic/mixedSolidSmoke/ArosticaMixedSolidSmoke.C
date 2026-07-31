/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Application
    ArosticaMixedSolidSmoke

Description
    Instantiates the minimal Aróstica mixed solid model through the
    solidModel runtime-selection table.

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

    if
    (
        solid->type()
     != "arosticaNonLinearGeometryTotalLagrangianTotalDisplacement"
    )
    {
        FatalErrorInFunction
            << "Unexpected solid-model runtime type " << solid->type()
            << abort(FatalError);
    }

    Info<< "PASS: Aróstica mixed solid runtime selection" << nl
        << "    type = " << solid->type() << endl;

    return 0;
}


// ************************************************************************* //
