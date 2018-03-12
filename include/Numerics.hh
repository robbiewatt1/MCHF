
#include "Vector.hh"
#include "Array3D.hh"

namespace Numerics
{
// Integrates the vector integrand between lower and upper points. requires an odd number
// of points in the vectroe
double SimpsonsRule(const Vector<double> &integrand, double lowerBound, double upperBound);

double SimpsonsRule3D(const Vector <double> xAxis, const Vector <double> yAxis,
                      const Vector <double> zAxis, const Array3D <double> data);

// 1D linear interpolation. if the quere points are outside of the range of know points then
// a linear extrapolation method will be used
Vector<double> Interpolate1D(const Vector<double> &samplePoints, const Vector<double> &sampleData,
                             const Vector<double> &queryPoints);

}