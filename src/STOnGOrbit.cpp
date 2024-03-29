#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "STOnGOrbit.hh"
#include "GaussianOrbital.hh"

STOnGOrbit::STOnGOrbit():
	m_k(0), m_m(0), m_n(0), m_orbitPosition(3)
{
}

STOnGOrbit::STOnGOrbit(std::string dataFile, int k, int m, int n, Vector<double> orbitPosition):
	m_k(k), m_m(m), m_n(n)
{
	m_orbitPosition = orbitPosition;
	OpenDataFile(dataFile);
	AllocateBaseOrbitals();
	Normalise();
}

STOnGOrbit::~STOnGOrbit()
{
}

int STOnGOrbit::GetK() const
{
	return m_k;
}

int STOnGOrbit::GetM() const
{
	return m_m;
}

int STOnGOrbit::GetN() const
{
	return m_n;
}

double STOnGOrbit::GetNormaliseConstant() const
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

double STOnGOrbit::Overlap(const STOnGOrbit &orbit) const
{
	double overlap(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			overlap += m_baseOrbitalVector[i].Overlap(orbit.GetBaseOribtal(j)) 
						* m_coefficeints[i] * orbit.m_coefficeints[j]
						* m_normaliseConstant * orbit.m_normaliseConstant
						* m_baseOrbitalVector[i].GetNormaliseConstant()
						* orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return overlap;
}

double STOnGOrbit::KineticOverlap(const STOnGOrbit &orbit) const
{
	double knieticEnergy(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			knieticEnergy += m_baseOrbitalVector[i].KineticOverlap(orbit.GetBaseOribtal(j)) 
							 * m_coefficeints[i] * orbit.m_coefficeints[j]
							 * m_normaliseConstant * orbit.m_normaliseConstant
							 * m_baseOrbitalVector[i].GetNormaliseConstant()
							 * orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return knieticEnergy;
}

double STOnGOrbit::NuclearOverlap(const STOnGOrbit &orbit, int nuclearCharge, 
								  const Vector<double> &nuclearPosition, const BoysFunction &boyFn) const
{
	double potentialEnergy(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			potentialEnergy += m_baseOrbitalVector[i].NuclearOverlap(orbit.GetBaseOribtal(j), 
									nuclearCharge, nuclearPosition, boyFn)
							 * m_coefficeints[i] * orbit.m_coefficeints[j]
							 * m_normaliseConstant * orbit.m_normaliseConstant
							 * m_baseOrbitalVector[i].GetNormaliseConstant()
							 * orbit.GetBaseOribtal(j).GetNormaliseConstant();
		}
	}
	return potentialEnergy;
}

Vector<double> STOnGOrbit::MatrixElement(const STOnGOrbit &orbit) const
{
	Vector<double> matrixElements(3);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < orbit.m_gaussianNumber; j++)
		{
			matrixElements = matrixElements + m_baseOrbitalVector[i].MatrixElement(orbit.GetBaseOribtal(j))
						 	* m_coefficeints[i] * orbit.m_coefficeints[j]
						 	* m_normaliseConstant * orbit.m_normaliseConstant
							* m_baseOrbitalVector[i].GetNormaliseConstant()
							* orbit.GetBaseOribtal(j).GetNormaliseConstant();						 	
		}
	}
	return matrixElements;
}

Array3D<double> STOnGOrbit::CalculateDataCartesian(const Vector<double> &xAxis,
	    const Vector<double> &yAxis,
        const Vector<double> &zAxis)
{
	SetXAxis(xAxis);
	SetYAxis(yAxis);
	SetZAxis(zAxis);
	Array3D<double> data(xAxis.Length(), yAxis.Length(), zAxis.Length());

	for (int i = 0; i < m_gaussianNumber; i++)
	{
		// Calculate the data for each primitive guassian
		Array3D<double> orbitData = m_normaliseConstant * m_baseOrbitalVector[i].CalculateDataCartesian(xAxis, yAxis, zAxis) ;
		data = data + (m_coefficeints[i] * orbitData);
	}
	return data;
}

Array3D<double> STOnGOrbit::CalculateDataSpherical(const Vector<double> &rAxis,
	    const Vector<double> &thetaAxis,
        const Vector<double> &phiAxis)
{
	
}

void STOnGOrbit::Normalise()
{
	double overlap(0);
	for (int i = 0; i < m_gaussianNumber; i++)
	{
		for (int j = 0; j < m_gaussianNumber; j++)
		{
			overlap += m_baseOrbitalVector[i].Overlap(m_baseOrbitalVector[j]) 
						* m_coefficeints[i] * m_coefficeints[j]
						* m_baseOrbitalVector[i].GetNormaliseConstant()
						* m_baseOrbitalVector[j].GetNormaliseConstant();
		}
	}
	m_normaliseConstant = std::sqrt(1.0 / overlap);
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
		m_baseOrbitalVector.Append(GaussianOrbital(m_k, m_m, m_n, m_exponents[i], m_orbitPosition));
	}
}