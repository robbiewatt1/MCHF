#include "STO6GOrbit.hh"
#include "GaussianOrbital.hh"

STO6GOrbit::STO6GOrbit():
	m_k(0), m_m(0), m_n(0), centrePositionX(0), m_centrePositionY(0), m_centrePositionZ(0),
	centrePositionZ(0)
{
	Normalise();
	AllocateBaseOrbitals();
}

STO6GOrbit::STO6GOrbit(int k, int m, int n, double centrePositionX, double centrePositionY,
			   			double centrePositionZ):
	m_k(k), m_m(m), m_n(n), centrePositionX(centrePositionX), m_centrePositionY(centrePositionY),
	m_centrePositionZ(centrePositionZ)
{
	Normalise();
	AllocateBaseOrbitals();
}

STO6GOrbit::~STO6GOrbit()
{
}

int STO6GOrbit::GetK()
{
	return m_k;
}

int STO6GOrbit::GetM()
{
	return m_m;
}

int STO6GOrbit::GetN()
{
	return m_n;
}

double STO6GOrbit::GetAlpha()
{
	return m_alpha;
}

double STO6GOrbit::GetPositionX()
{
	return m_centrePositionX;
}

double STO6GOrbit::GetNormaliseConstant()
{
	return m_normaliseConstant;
}

GaussianOrbital STO6GOrbit::GetBaseOribtal(int i);
{
	return baseOrbitalVector[i];
}

double STO6GOrbit::GetCoefficient(int i)
{
	return m_coefficeints[i];
}

double STO6GOrbit::Overlap(const STO6GOrbit &orbit)
{
	double overlap(0);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			overlap += baseOrbitalVector[i].Overlap(orbit.GetBaseOribtal(j)) * m_coefficeints[i]
					   orbit.m_coefficeints[j];
		}
	}
	return overlap;
}


STO6GOrbit::AllocateBaseOrbitals()
{
	for (int i = 0; i < 6; i++)
	{
		baseOrbitalVector[i] = GaussianOrbital(m_k, m_m, m_n, m_exponents[i], m_centrePositionX, 
											   m_centrePositionY, m_centrePositionZ);
	}
}

