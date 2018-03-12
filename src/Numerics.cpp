#include <cmath>
#include <iostream>
#include <cassert>

#include "Numerics.hh"
#include "Vector.hh"
#include "Array3D.hh"

double Numerics::SimpsonsRule(const Vector<double> &integrand, double lowerBound, double upperBound)
{
	assert(integrand.Length() % 2 != 0);

	int steps = integrand.Length() - 1;
	double deltaX = (upperBound - lowerBound) / steps;
	double integral = integrand(0);

	for (int i = 1; i < steps; i += 2)
	{
		integral += 4 * integrand(i);
	}
	for (int i = 2; i < steps - 1; i += 2)
	{
		integral += 2 * integrand(i);
	}
	integral += integrand(steps);

	return (integral * deltaX / 3);
}

double Numerics::SimpsonsRule3D(const Vector <double> xAxis, const Vector <double> yAxis,
                                const Vector <double> zAxis, const Array3D <double> data)
{
	double integralX(0), integralY(0), integralZ(0);
	double integralX0(0), integralY0(0), integralZ0(0);
	double deltaX = (xAxis.End() - xAxis.Start()) / xAxis.Length();
	double deltaY = (yAxis.End() - yAxis.Start()) / yAxis.Length();
	double deltaZ = (zAxis.End() - zAxis.Start()) / zAxis.Length();

	for (int k = 1; k < zAxis.Length(); k++)
	{
		integralZ0 += deltaZ * (data[0][0][k - 1] + data[0][0][k]) / 2;
	}
	for (int j = 1; j < yAxis.Length(); j++)
	{
		integralZ = 0;
		for (int k = 1; k < zAxis.Length(); k++)
		{
			integralZ += deltaZ * (data[0][j][k - 1] + data[0][j][k]) / 2;
		}
		integralY0 += deltaY * (integralZ0 + integralZ) / 2;
		integralZ0 = integralZ;
	}

	for (int i = 1; i < xAxis.Length(); i++)
	{
		integralZ0 = 0;
		integralY = 0;
		for (int k = 1; k < zAxis.Length(); k++)
		{
			integralZ0 += deltaZ * (data[i][0][k - 1] + data[i][0][k]) / 2;
		}
		for (int j = 1; j < yAxis.Length(); j++)
		{
			integralZ = 0;
			for (int k = 1; k < zAxis.Length(); k++)
			{
				integralZ += deltaZ * (data[i][j][k - 1] + data[i][j][k]) / 2;
			}
			integralY += deltaY * (integralZ0 + integralZ) / 2;
			integralZ0 = integralZ;
		}
		integralX += deltaX * (integralY0 + integralY) / 2;
		integralY0 = integralY;
	}
	return integralX;
}

Vector<double> Numerics::Interpolate1D(const Vector<double> &samplePoints, const Vector<double> &sampleData, const Vector<double> &queryPoints)
{
	Vector<double> queryData(queryPoints.Length());
	for (int i = 0; i < queryPoints.Length(); i++)
	{
		if (queryPoints(i) < samplePoints(0))
		{
			// value is not within the range so we need to use extrapolation
			queryData[i] = sampleData(0) + ((queryPoints(i) - samplePoints(0)) / (samplePoints(1) - samplePoints(0))) * (sampleData(1) - sampleData(0));
		} else if (queryPoints(i) > samplePoints.End())
		{
			// Value is again not within the range so we need to extrapolate
			int sampleEnd = samplePoints.Length() - 1;
			queryData[i] = sampleData(sampleEnd - 1) + ((queryPoints(i) - samplePoints(sampleEnd - 1)) * (sampleData(sampleEnd) - sampleData(sampleEnd - 1))
			               / (samplePoints(sampleEnd) - samplePoints(sampleEnd - 1)));
			std::cout << samplePoints(sampleEnd) << std::endl;
		} else
		{
			int lowerIndex, higherIndex;
			for (int j = 0; j < samplePoints.Length(); ++j)
			{
				if (queryPoints(i) < samplePoints(j))
				{
					lowerIndex  = j - 1;
					higherIndex = j;
					break;
				}
			}
			queryData[i] = sampleData(lowerIndex) + (queryPoints(i) - samplePoints(lowerIndex)) * (sampleData(higherIndex) - sampleData(lowerIndex))
			               / (samplePoints(higherIndex) - samplePoints(lowerIndex));
			std::cout << (samplePoints(higherIndex) - samplePoints(lowerIndex)) << std::endl;
		}
	}
	return queryData;
}