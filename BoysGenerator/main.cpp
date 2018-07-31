#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <omp.h>


static double pi = 3.14159265359;

double SemiFactorial(double n)
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

double Factorial(double n)
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

double ErrorFunction(double x)
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

void BoysFunction(int *vAxis, int vMax, double *uAxis, int uMax, double **result)
{
	#pragma omp parallel for
	for(int i = 0; i < uMax; i++)
	{
		double result_p;
//		if (uAxis[i] > 25.0)  //check this number
//		{
//			result_p = SemiFactorial(2.0 * vMax - 1.0) / (std::pow(2.0, vMax + 1.0))
//							* std::sqrt(pi / (std::pow(uAxis[i], 2.0 * vMax + 1.0)));
//			result_p = SemiFactorial(2 * vMax - 1) * std::sqrt(pi / std::pow(uAxis[i], 2 * vMax + 1)) / std::pow(2.0, vMax + 1);
		if ( uAxis[i] > 1.0)
		{
			double sum = 0;
			for(int j = 0; j < 20; j++)
			{
				sum += SemiFactorial(2 * vMax - 1) * std::pow(2 * uAxis[i], j) / SemiFactorial(2 * vMax + 2 * j + 1);
			}
			result_p = std::exp( -1.0 * uAxis[i]) * sum;
			if(result_p != result_p)
			{
				std::cout << SemiFactorial(2 * vMax + 2 * 200 + 1) << std::endl;
			}
/*
			double factor1 = Factorial(2.0 * vMax) / (2.0 *  Factorial(vMax));
                        double factor2 = std::sqrt(pi) * ErrorFunction(std::sqrt(uAxis[i]))
                                 / (std::pow(4.0, vMax) * std::pow(uAxis[i], vMax + 0.5));

                        float sum = 0;
                        for (int j = 0; j < vMax; j++)
                        {
                                sum += Factorial(1.0 * vMax - j) / (std::pow(4.0, j) * (Factorial(2.0 * vMax - 2.0 * j))
                                                                                 * std::pow(uAxis[i], j + 1.0));
                        }
                        result_p = factor1 * ( factor2 - std::exp(-1.0 * uAxis[i]) * sum);
*/
		} else
		{
			result_p = (1.0 / (2.0 * vMax + 1.0))
				  - (uAxis[i] / (2.0 * vMax + 3.0))
		 		  + (uAxis[i] * uAxis[i] / (2.0 * (2.0 * vMax + 5.0)))
				  - (uAxis[i] * uAxis[i] * uAxis[i] / (6.0 * (2.0 * vMax + 7.0)))
				  + (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (24.0 * (2.0 * vMax + 11.0)))
				  - (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (120.0 * (2.0 * vMax + 13.0)))
				  + (uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] * uAxis[i] / (720.0 * (2.0 * vMax + 15.0)));
		}
		result[i][0] = result_p;
		for(int k = 1; k < vMax; k++)
		{
			result[i][k] = (2 * uAxis[i] * result[i][k-1] + std::exp(-uAxis[i]))
					/ (2 * (vMax - k - 1) + 1);
		}
	}
}
		
int main(int argc, char *argv[])
{

	std::ofstream file("./BoysData");

	// Set the v axis
	int vAx[30];
	double *uAx = new double[1000];
	double **data =  new double *[1000];

	for(int i = 0; i < 1000; i++)
	{
		data[i] = new double [30];
	}

	for(int i = 0; i < 30; i++)
	{
		vAx[i] = i;
	}

	for(int i = 0; i < 1000; i++)
	{
		uAx[i] = 30 * (double)i/1000.0;
	}

	//BoysFunction(vAx,30,uAx,1000,data);
/*
	for(int i = 0; i < 1000; i++)
	{
		for(int j = 0; j < 30; j++)
		{
			file << data[i][j] << "\t";
		}
		file << "\n";
	}
*/
	return 0;
}
