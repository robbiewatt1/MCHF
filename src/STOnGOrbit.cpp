#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "STOnGOrbit.hh"
#include "GaussianOrbital.hh"

STOnGOrbit::STOnGOrbit():
	m_k(0), m_m(0), m_n(0), m_centrePositionX(0), m_centrePositionY(0), m_centrePositionZ(0)
{
}

STOnGOrbit::STOnGOrbit(std::string dataFile, int k, int m, int n, double centrePositionX,
                       double centrePositionY, double centrePositionZ):
	m_k(k), m_m(m), m_n(n), m_centrePositionX(centrePositionX), m_centrePositionY(centrePositionY),
	m_centrePositionZ(centrePositionZ)
{
	OpenDataFile(dataFile);
	AllocateBaseOrbitals();
	Normalise();				// Normalise here, although should already be normed
}

STOnGOrbit::~STOnGOrbit()
{
}

double STOnGOrbit::GetNormaliseConstant()
{
	return m_normaliseConstant;
}

GaussianOrbital STOnGOrbit::GetBaseOribtal(int i) const
{
	return m_baseOrbitalVector[i];
}

double STOnGOrbit::GetExponent(int i) const
{
	return m_exponents[i];
}

double STOnGOrbit::GetCoefficient(int i) const
{
	return m_coefficeints[i];
}

double STOnGOrbit::Overlap(const STOnGOrbit &orbit)
{
	double overlap(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			overlap += m_baseOrbitalVector[i].Overlap(orbit.GetBaseOribtal(j)) 
						* m_coefficeints[i] * orbit.m_coefficeints[j] 
						* m_baseOrbitalVector[i].GetNormaliseConstant()
						* orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return overlap;
}

double STOnGOrbit::KineticOverlap(const STOnGOrbit &orbit)
{
	double knieticEnergy(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			knieticEnergy += m_baseOrbitalVector[i].KineticOverlap(orbit.GetBaseOribtal(j)) 
							 * m_coefficeints[i] * orbit.m_coefficeints[j]
							 * m_baseOrbitalVector[i].GetNormaliseConstant()
							 * orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return knieticEnergy;
}

double STOnGOrbit::NuclearOverlap(const STOnGOrbit &orbit, int nuclearCharge, double nuclearX,
								  double nuclearY, double nuclearZ)
{
	double potentialEnergy(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			potentialEnergy += m_baseOrbitalVector[i].NuclearOverlap(orbit.GetBaseOribtal(j), 
									nuclearCharge, nuclearX, nuclearY, nuclearZ)
							 * m_coefficeints[i] * orbit.m_coefficeints[j]
							 * m_baseOrbitalVector[i].GetNormaliseConstant()
							 * orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return potentialEnergy;
}

void STOnGOrbit::CalculateDataCartesian(const Vector<double> &xAxis,
	    const Vector<double> &yAxis,
        const Vector<double> &zAxis)
{

}

void STOnGOrbit::CalculateDataSpherical(const Vector<double> &rAxis,
	    const Vector<double> &thetaAxis,
        const Vector<double> &phiAxis)
{
	
}

void STOnGOrbit::Normalise()
{
	m_normaliseConstant = std::sqrt(1.0 / Overlap(*this));
}

void STOnGOrbit::OpenDataFile(std::string fileName)
{
	std::ifstream file(fileName);

	if (!file)
	{
		std::cerr << "File " << fileName << " not found" << std::endl;
		std::exit(1);
	}
	
	double exponent, coefficeint;
	m_gaussianNumber = 0;
	while (file >> exponent >> coefficeint)
	{
		m_exponents.Append(exponent);
		m_coefficeints.Append(coefficeint);
		++m_gaussianNumber;
	}
}

void STOnGOrbit::AllocateBaseOrbitals()
{
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		m_baseOrbitalVector.Append(GaussianOrbital(m_k, m_m, m_n, m_exponents[i], m_centrePositionX,
		                                       m_centrePositionY, m_centrePositionZ));
	}
}