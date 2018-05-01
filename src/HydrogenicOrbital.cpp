/*
#include <cassert>

#include "HydrogenicOrbital.hh"
#include "Matrix.hh"

HydrogenicOrbital::HydrogenicOrbital()
{
}

HydrogenicOrbital::HydrogenicOrbital(int maxL):
	m_maxL(maxL)
{
	m_coefficents = new double *[maxL + 1];
	for (int i = 0; i < maxL + 1; i++)
	{
		m_coefficents[i] = new double[(i + 1) * (i + 2) / 2];
	}
}

HydrogenicOrbital::~HydrogenicOrbital()
{
	delete[] m_coefficents;
}

void HydrogenicOrbital::SetCoefficent(int l, int ml, double value)
{
	assert(ml <= 2 * l + 1);
	m_coefficents[l][ml] = value;
}

int HydrogenicOrbital::GetMaxL()
{
	return m_maxL;
}

void HydrogenicOrbital::MinimiseEnergy()
{
	// First need to form the H matrix
	// Thej need to form the S matrix
	// then need to do eiganvalue problem
}

void HydrogenicOrbital::CalculateDataCartesian(const Vector<double> &xAxis,
	                            			   const Vector<double> &yAxis,
	                            			   const Vector<double> &zAxis)
{
}

void HydrogenicOrbital::CalculateDataSpherical(const Vector<double> &rAxis,
			   	                               const Vector<double> &thetaAxis,
	                            			   const Vector<double> &phiAxis)
{
}

Matrix<double> HydrogenicOrbital::CalculateSMatrix()
{
	// Should be symetic matrix
	int totalOrbits = NumberOfOrbitals();
	Matrix<double> SMaxtrix(totalOrbits,totalOrbits);
	for (int i = 0; i < totalOrbits; i++)
	{
		for (int j = 0; j <= i; j++)
		{


//			S[i][j] = value;
		}	
	}
	return SMaxtrix;
	// fill in the other half of the matrix
}

int HydrogenicOrbital::NumberOfOrbitals()
{

	int totalOrbits = 0;
	for (int i = 0; i < m_maxL + 1; ++i)
	{
		totalOrbits += (i + 1) * (i + 2) / 2;
	}
	return totalOrbits;
}
*/