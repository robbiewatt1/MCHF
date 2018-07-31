#include <cassert>
#include <cmath>
#include <iostream>

#include <omp.h>
#include "Functions.hh"
#include "Constants.hh"
#include "Matrix.hh"
#include "Vector.hh"

double Functions::Factorial(double n)
{
	assert(n >= 0);
	double factorial = 1;
	if (n > 1)
	{
		for (int i = 0; i < n; ++i)
		{
			factorial *= (n - i);
		}
	}
	return factorial;
}

double Functions::SemiFactorial(double n)
{
	double semiFactorial = 1;
	if (n > 1)
	{
		for (int i = n; i > 1; i -= 2)
		{
			semiFactorial *= i;
		}
	}
	return semiFactorial;
}

unsigned int Functions::Gamma(int n)
{
	return (Factorial(n - 1));
}

double Functions::GammaPluseHalf(int n)
{
	double gamma  = std::pow(Constants::pi, 0.5)
	                * SemiFactorial((2 * n) - 1.0)
	                / (std::pow(2.0, n));
	return gamma;
}

int Functions::BinomialCoefficient(int n, int k)
{
	int binCoef = Factorial(n) / (Factorial(k) * Factorial(n - k));
	return binCoef;
}

double Functions::ErrorFunction(double x)
{
	double parity(1.0);
	if (x < 0)
	{
		parity = -1.0;
		x = std::abs(x);
	}
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p  = 0.3275911;
	double t = 1.0/(1.0 + (p * x));
	double erf = parity * (1.0 - ((a1 * t) + (a2 * std::pow(t, 2.0)) + (a3 * std::pow(t, 3.0)) + (a4 * std::pow(t, 4.0)) 
				 + (a5 * std::pow(t, 5.0))) * std::exp(- (x * x)));
	return erf;
}

Matrix<double> Functions::BoysGenerator(int maxV, double maxU, int deltaU)
{
	Matrix<double> boysResult(maxV, deltaU);
	Vector<double> uAxis(deltaU);

	for (int i = 0; i < deltaU; i++)
	{
		uAxis[i] =  maxU * (double)i / deltaU; 
	}
	for(int i = 0; i < deltaU; i++)
	{
		double result_p;
		if ( uAxis[i] > 1.0)
		{
			double sum = 0;
			for(int j = 0; j < 20; j++)
			{
				sum += SemiFactorial(2 * maxV - 1) * std::pow(2 * uAxis[i], j) / SemiFactorial(2 * maxV + 2 * j + 1);
			}
			result_p = std::exp( -1.0 * uAxis[i]) * sum;
		} else
		{
			result_p = (1.0 / (2.0 * maxV + 1.0))
					 - (uAxis[i] / (2.0 * maxV + 3.0))
					 + (uAxis[i] * uAxis[i] / (2.0 * (2.0 * maxV + 5.0)))
					 - (uAxis[i] * uAxis[i] * uAxis[i] / (6.0 * (2.0 * maxV + 7.0)))
					 + (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (24.0 * (2.0 * maxV + 11.0)))
					 - (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (120.0 * (2.0 * maxV + 13.0)))
					 + (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (720.0 * (2.0 * maxV + 15.0)));
		}
		boysResult[0][i] = result_p;
		for(int k = 1; k < maxV; k++)
		{
			boysResult[k][i] = (2 * uAxis[i] * boysResult[k-1][i] + std::exp(-uAxis[i]))
						 / (2 * (maxV - k - 1) + 1);
		}
	}
	return boysResult;
}

double Functions::BoysFunction(int v, double u)
{
	double result(0);
	if(u > 1.5) // Check this number
	{
		double factor1 = Factorial(2.0 * v) / (2.0 *  Factorial(v));
		double factor2 = std::sqrt(Constants::pi) * ErrorFunction(std::sqrt(u))
	                 / (std::pow(4.0, v) * std::pow(u, v + 0.5));
		double sum(0);
		for (int i = 0; i < v; i++)
		{
			sum += Factorial(1.0 * v - i) / (std::pow(4.0, i) * (Factorial(2.0 * v - 2.0 * i))
		   	                                   * std::pow(u, i + 1.0));
		}
		result = factor1 * ( factor2 - std::exp(-1.0 * u) * sum);
	} else
	{
		result = (1.0 / (2.0 * v + 1.0)) - (u / (2.0 * v + 3.0)) + (u * u/(2.0 * (2.0 * v + 5.0))) 
				 - (u * u * u/(6.0 * (2.0 * v + 7.0))) + (u * u * u * u/(24 * (2.0 * v + 11.0)))
				 - (u * u * u * u/(120.0 * (2.0 * v + 11.0)));
	}
	if (result < 0.0 || result > 1.0)
	{
		std::cerr << u << " EERRROORR " << result << std::endl;
	}
	return result;
}

