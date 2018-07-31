#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include "Matrix.hh"


namespace Functions
{
	// Returns the factoriaol of int n
	double Factorial(double n);

	// Returns the seminfactorial of int n
	double SemiFactorial(double n);

	unsigned int Gamma(int n);
	// Returns the gamma function for an interger

	double GammaPluseHalf(int n);
	// Return the gamma function of (n + 0.5)

	int BinomialCoefficient(int n, int k);
	// Return the binomial coefficent given by n!/(k!(n-k)!)

	double ErrorFunction(double x);
	//Returns the error function of x to a minimum error of 1.5e-7

	Matrix<double> BoysGenerator(int maxV, double maxU, int deltaU);
	// Function to generate tabulated boys function and returned as a matrix.

	double BoysFunction(int v, double u);
}
#endif