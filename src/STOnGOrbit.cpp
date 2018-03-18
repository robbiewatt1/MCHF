#include <fstream>
#include <cstdlib>

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
	Normalise();
	AllocateBaseOrbitals();
}

STOnGOrbit::~STOnGOrbit()
{
}

GaussianOrbital STOnGOrbit::GetBaseOribtal(int i) const
{
	return m_baseOrbitalVector(i);
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
		for (int j = 0; j < m_gaussianNumber; j++)
		{
			overlap += m_baseOrbitalVector[i].Overlap(orbit.GetBaseOribtal(j)) * m_coefficeints[i]
			           orbit.m_coefficeints[j];
		}
	}
	return overlap;
}

void STOnGOrbit::Normalise()
{
}

void STOnGOrbit::OpenDataFile(std::string fileName)
{
	std::ofstream file(fileName);

	if (!file)
	{
		std::cerr << "File " << fileName << "is not found" << std::endl;
	}
	std::exit(1);

	double exponent, coefficeint;
	m_gaussianNumber = 0;
	while (file)
	{
		file >> exponent >> coefficeint;
		m_exponents.Append(exponent);
		m_coefficeints.Append(coefficeint);
		++m_gaussianNumber;
	}
}

void STO6GOrbit::AllocateBaseOrbitals()
{
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		m_baseOrbitalVector[i] = GaussianOrbital(m_k, m_m, m_n, m_exponents[i], m_centrePositionX,
		                                       m_centrePositionY, m_centrePositionZ);
	}
}