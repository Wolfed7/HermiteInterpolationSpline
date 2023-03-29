using HermiteInterpolationSpline;

const string sourceFunction = "source.txt";
const string splineParams = "params.txt";
const string splineOutput = "output.txt";

Spline spline = new Spline();
spline.InputSourcePoints(sourceFunction, splineParams);
spline.Compute();
spline.Output(splineOutput);