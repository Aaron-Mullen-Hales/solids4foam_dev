/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

Application
    ArosticaFailingFaceRegression

Description
    Validates the exact Case-B boundary-face state retained for the diagnostic
    lifecycle regression.  The invalid tensor is deliberately not evaluated:
    the test proves that its sheet exponential is the first overflowing
    operation.  It then evaluates the restored central-FD minus tensor and
    checks the finite passive-plus-viscous Cauchy stress recorded in the data.

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "IOstreams.H"
#include "dictionary.H"
#include "fileName.H"
#include "symmTensor.H"
#include "tensor.H"
#include "vector.H"

#include <cmath>
#include <limits>

using namespace Foam;


namespace
{

bool finiteTensor(const tensor& value)
{
    for (direction componentI = 0; componentI < tensor::nComponents; ++componentI)
    {
        if (!std::isfinite(value[componentI]))
        {
            return false;
        }
    }
    return true;
}


bool finiteSymmTensor(const symmTensor& value)
{
    for
    (
        direction componentI = 0;
        componentI < symmTensor::nComponents;
        ++componentI
    )
    {
        if (!std::isfinite(value[componentI]))
        {
            return false;
        }
    }
    return true;
}


void requireClose
(
    const word& name,
    const scalar actual,
    const scalar expected,
    const scalar absoluteTolerance,
    const scalar relativeTolerance
)
{
    const scalar scale = max(scalar(1.0), max(mag(actual), mag(expected)));
    const scalar tolerance = absoluteTolerance + relativeTolerance*scale;
    if (mag(actual - expected) > tolerance)
    {
        FatalErrorInFunction
            << name << " differs: actual=" << actual
            << ", expected=" << expected
            << ", tolerance=" << tolerance
            << abort(FatalError);
    }
}


struct SwitchResult
{
    scalar value;
    scalar derivative;
};


SwitchResult logisticSwitch(const scalar invariant, const scalar steepness)
{
    const scalar argument = -steepness*(invariant - 1.0);
    const scalar value = 1.0/(1.0 + std::exp(argument));
    return SwitchResult{value, steepness*value*(1.0 - value)};
}


symmTensor passiveStress
(
    const tensor& F,
    const vector& f0,
    const vector& s0,
    const dictionary& parameters
)
{
    const scalar a(readScalar(parameters.lookup("a")));
    const scalar b(readScalar(parameters.lookup("b")));
    const scalar af(readScalar(parameters.lookup("af")));
    const scalar bf(readScalar(parameters.lookup("bf")));
    const scalar as(readScalar(parameters.lookup("as")));
    const scalar bs(readScalar(parameters.lookup("bs")));
    const scalar afs(readScalar(parameters.lookup("afs")));
    const scalar bfs(readScalar(parameters.lookup("bfs")));
    const scalar switchK
    (
        readScalar(parameters.lookup("compressionSwitchK"))
    );

    const scalar J = det(F);
    const symmTensor C(symm(F.T() & F));
    const scalar traceC = tr(C);
    const scalar Jm23 = std::pow(J, -2.0/3.0);
    const scalar I1bar = Jm23*traceC;
    const scalar I4f = f0 & (C & f0);
    const scalar I4s = s0 & (C & s0);
    const scalar I8fs = f0 & (C & s0);
    const SwitchResult fibreSwitch = logisticSwitch(I4f, switchK);
    const SwitchResult sheetSwitch = logisticSwitch(I4s, switchK);
    const scalar fibreExponent = std::exp(bf*sqr(I4f - 1.0));
    const scalar sheetExponent = std::exp(bs*sqr(I4s - 1.0));
    const scalar fibreSheetExponent = std::exp(bfs*sqr(I8fs));

    const symmTensor dI1bar
    (
        Jm23*(symmTensor::I - (traceC/3.0)*symm(inv(C)))
    );
    symmTensor S(a*std::exp(b*(I1bar - 3.0))*dI1bar);
    const scalar qf = I4f - 1.0;
    const scalar qs = I4s - 1.0;
    const scalar dWf =
        af/(2.0*bf)
       *(
            fibreSwitch.derivative*(fibreExponent - 1.0)
          + fibreSwitch.value*fibreExponent*(2.0*bf*qf)
        );
    const scalar dWs =
        as/(2.0*bs)
       *(
            sheetSwitch.derivative*(sheetExponent - 1.0)
          + sheetSwitch.value*sheetExponent*(2.0*bs*qs)
        );
    S += 2.0*dWf*symm(f0*f0);
    S += 2.0*dWs*symm(s0*s0);
    S += 2.0*afs*I8fs*fibreSheetExponent*symm(f0*s0);
    return symm(F & S & F.T())/J;
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        FatalErrorInFunction
            << "Usage: ArosticaFailingFaceRegression <caseB-failing-face.dict>"
            << abort(FatalError);
    }

    const fileName dataFile(argv[1]);
    IFstream input(dataFile);
    if (!input.good())
    {
        FatalErrorInFunction
            << "Cannot read " << dataFile << abort(FatalError);
    }
    const dictionary data(input);
    const dictionary& failure = data.subDict("failureState");
    const dictionary& restored = data.subDict("restoredMinusState");
    const dictionary& directions = data.subDict("directions");
    const dictionary& histories = data.subDict("histories");
    const dictionary& parameters = data.subDict("materialParameters");
    const dictionary& expected = data.subDict("expected");

    const tensor lawOwnedF(failure.lookup("lawOwnedF"));
    const tensor independentlyReconstructedF
    (
        failure.lookup("independentlyReconstructedF")
    );
    const vector f0(directions.lookup("f0f"));
    const vector s0(directions.lookup("s0f"));
    const vector n0(directions.lookup("n0f"));
    const scalar bs(readScalar(parameters.lookup("bs")));
    const scalar bf(readScalar(parameters.lookup("bf")));
    const scalar bfs(readScalar(parameters.lookup("bfs")));

    if
    (
        !finiteTensor(lawOwnedF)
     || !finiteTensor(independentlyReconstructedF)
     || mag(lawOwnedF - independentlyReconstructedF) != 0.0
    )
    {
        FatalErrorInFunction
            << "Retained law/reconstruction tensors are not finite and equal"
            << abort(FatalError);
    }

    const symmTensor failureC(symm(lawOwnedF.T() & lawOwnedF));
    const scalar failureI4f = f0 & (failureC & f0);
    const scalar failureI4s = s0 & (failureC & s0);
    const scalar failureI8fs = f0 & (failureC & s0);
    const scalar fibreArgument = bf*sqr(failureI4f - 1.0);
    const scalar sheetArgument = bs*sqr(failureI4s - 1.0);
    const scalar fibreSheetArgument = bfs*sqr(failureI8fs);
    const scalar maximumFiniteExpArgument =
        std::log(std::numeric_limits<scalar>::max());

    requireClose
    (
        "failure I4f",
        failureI4f,
        readScalar(failure.lookup("I4f")),
        1.0e-11,
        1.0e-12
    );
    requireClose
    (
        "failure I4s",
        failureI4s,
        readScalar(failure.lookup("I4s")),
        1.0e-10,
        1.0e-12
    );
    requireClose
    (
        "failure I8fs",
        failureI8fs,
        readScalar(failure.lookup("I8fs")),
        1.0e-10,
        1.0e-12
    );

    if
    (
        !std::isfinite(std::exp(fibreArgument))
     || sheetArgument <= maximumFiniteExpArgument
     || !std::isfinite(std::exp(fibreSheetArgument))
    )
    {
        FatalErrorInFunction
            << "The retained state no longer identifies sheet exp overflow "
            << "as the first invalid operation" << abort(FatalError);
    }

    const tensor restoredF(restored.lookup("lawOwnedF"));
    const tensor restoredReconstruction
    (
        restored.lookup("independentlyReconstructedF")
    );
    if
    (
        !finiteTensor(restoredF)
     || mag(restoredF - restoredReconstruction) != 0.0
     || det(restoredF) <= 0.0
    )
    {
        FatalErrorInFunction
            << "Restored central-FD minus tensor is invalid"
            << abort(FatalError);
    }

    const symmTensor restoredC(symm(restoredF.T() & restoredF));
    const scalar restoredI1bar =
        std::pow(det(restoredF), -2.0/3.0)*tr(restoredC);
    const scalar restoredI4f = f0 & (restoredC & f0);
    const scalar restoredI4s = s0 & (restoredC & s0);
    const scalar restoredI8fs = f0 & (restoredC & s0);
    requireClose
    (
        "restored J",
        det(restoredF),
        readScalar(restored.lookup("J")),
        1.0e-12,
        1.0e-12
    );
    requireClose
    (
        "restored I1bar",
        restoredI1bar,
        readScalar(restored.lookup("I1bar")),
        1.0e-12,
        1.0e-12
    );
    requireClose
    (
        "restored I4f",
        restoredI4f,
        readScalar(restored.lookup("I4f")),
        1.0e-12,
        1.0e-12
    );
    requireClose
    (
        "restored I4s",
        restoredI4s,
        readScalar(restored.lookup("I4s")),
        1.0e-12,
        1.0e-12
    );
    requireClose
    (
        "restored I8fs",
        restoredI8fs,
        readScalar(restored.lookup("I8fs")),
        1.0e-12,
        1.0e-12
    );

    const symmTensor EOld(histories.lookup("EOld"));
    const symmTensor EOldOld(histories.lookup("EOldOld"));
    const scalar current
    (
        readScalar(histories.lookup("currentCoefficient"))
    );
    const scalar old(readScalar(histories.lookup("oldCoefficient")));
    const scalar oldOld
    (
        readScalar(histories.lookup("oldOldCoefficient"))
    );
    const scalar eta(readScalar(parameters.lookup("eta")));
    const symmTensor E(0.5*(restoredC - symmTensor::I));
    const symmTensor Edot(current*E + old*EOld + oldOld*EOldOld);
    const symmTensor Sviscous(eta*Edot);
    const symmTensor sigmaViscous
    (
        symm(restoredF & Sviscous & restoredF.T())/det(restoredF)
    );
    const symmTensor sigma
    (
        passiveStress(restoredF, f0, s0, parameters) + sigmaViscous
    );
    const symmTensor expectedSigma(expected.lookup("finiteCauchyStress"));
    const scalar stressError = mag(sigma - expectedSigma);
    const scalar stressScale = max
    (
        scalar(1.0),
        max(mag(sigma), mag(expectedSigma))
    );
    const scalar stressTolerance = 1.0e-8 + 1.0e-10*stressScale;

    if
    (
        !finiteSymmTensor(sigma)
     || stressError > stressTolerance
     || mag(f0) <= 0.0
     || mag(s0) <= 0.0
     || mag(n0) <= 0.0
    )
    {
        FatalErrorInFunction
            << "Restored finite stress check failed: error="
            << stressError << ", tolerance=" << stressTolerance
            << abort(FatalError);
    }

    Info<< "FAILING_FACE_FIRST_INVALID_OPERATION"
        << " sheetArgument=" << sheetArgument
        << " maxFiniteExpArgument=" << maximumFiniteExpArgument
        << " source="
        << word(expected.lookup("firstInvalidSourceLocation")) << nl
        << "RESTORED_FACE_FINITE_RESULT"
        << " J=" << det(restoredF)
        << " stressError=" << stressError
        << " tolerance=" << stressTolerance << nl
        << "PASS: exact Case-B failing-face data regression" << endl;
    return 0;
}


// ************************************************************************* //
