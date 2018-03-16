#include <cmath>
#include <algorithm>

#include "GaussianOrbital.hh"
#include "Functions.hh"
#include "Constants.hh"
#include "Array3D.hh"

GaussianOrbital::GaussianOrbital():
	m_k(1), m_m(1), m_n(1), m_alpha(1), m_centrePositionX(0), m_centrePositionY(0),
	m_centrePositionZ(0)
{
	Normalise();
}

GaussianOrbital::GaussianOrbital(int k, int m, int n, double alpha, double centrePositionX,
                                 double centrePositionY, double centrePositionZ):
	m_k(k), m_m(m), m_n(n), m_alpha(alpha), m_centrePositionX(centrePositionX),
	m_centrePositionY(centrePositionY), m_centrePositionZ(centrePositionZ)
{
	Normalise();
}

GaussianOrbital::~GaussianOrbital()
{
}

int GaussianOrbital::GetK()
{
	return m_k;
}

int GaussianOrbital::GetM()
{
	return m_m;
}

int GaussianOrbital::GetN()
{
	return m_n;
}

double GaussianOrbital::GetAlpha()
{
	return m_alpha;
}

double GaussianOrbital::GetPositionX()
{
	return m_centrePositionX;
}

double GaussianOrbital::GetNormaliseConstant()
{
	return m_normaliseConstant;
}

double GaussianOrbital::Overlap(const GaussianOrbital &orbit)
{
	double gamma  = m_alpha + orbit.m_alpha;

	double posPx  = (m_centrePositionX + orbit.m_centrePositionX) / 2;
	double posPy  = (m_centrePositionY + orbit.m_centrePositionY) / 2;
	double posPz  = (m_centrePositionZ + orbit.m_centrePositionZ) / 2;
	double posAB2 = std::pow(orbit.m_centrePositionX - m_centrePositionX , 2.0)
	                + std::pow(orbit.m_centrePositionY - m_centrePositionY , 2.0)
	                + std::pow(orbit.m_centrePositionZ -  m_centrePositionZ, 2.0);

	double overlapX = OverlapFunction(m_k, orbit.m_k, gamma, m_centrePositionX - posPx,
	                                  orbit.m_centrePositionX - posPx);
	double overlapY = OverlapFunction(m_m, orbit.m_m, gamma, m_centrePositionY - posPy,
	                                  orbit.m_centrePositionY - posPy);
	double overlapZ = OverlapFunction(m_n, orbit.m_n, gamma, m_centrePositionZ - posPz,
	                                  orbit.m_centrePositionZ - posPz);
	double overlap = std::exp(-1 * m_alpha * orbit.m_alpha * posAB2 / gamma)
	                 * overlapX * overlapY * overlapZ;
	return overlap;
}

double GaussianOrbital::KineticOverlap(const GaussianOrbital &orbit)
{
	double normFactor = m_normaliseConstant * orbit.m_normaliseConstant;
	double overlap = this->Overlap(orbit);
	double plus2k  = this->Overlap(orbit.ChangeK(2));
	double plus2m  = this->Overlap(orbit.ChangeM(2));
	double plus2n  = this->Overlap(orbit.ChangeN(2));
	double minus2K(0), minus2M(0), minus2N(0);

	if (orbit.m_k > 1)
	{
		minus2K = this->Overlap(orbit.ChangeK(-2));
	}
	if (orbit.m_m > 1)
	{
		minus2M = this->Overlap(orbit.ChangeM(-2));
	}
	if (orbit.m_n > 1)
	{
		minus2N = this->Overlap(orbit.ChangeN(-2));
	}

	double kinOverlap = (orbit.m_alpha * (2 * (orbit.m_k + orbit.m_m + orbit.m_n) + 3) * overlap)
	                    - (2 * std::pow(orbit.m_alpha, 2.0) * (plus2k + plus2m + plus2n))
	                    - (0.5 * (orbit.m_k * (orbit.m_k - 1) * minus2K
	                              + orbit.m_m * (orbit.m_m - 1) * minus2M
	                              + orbit.m_n * (orbit.m_n - 1) * minus2N));
	return kinOverlap * normFactor / std::pow(Constants::pi,2.5);
}

double GaussianOrbital::NuclearOverlap(const GaussianOrbital &orbit, int nuclearCharge, double nuclearX,
                                       double nuclearY, double nuclearZ)
{
	double gamma  = m_alpha + orbit.m_alpha;

	double posPx = (m_centrePositionX + orbit.m_centrePositionX) / 2;
	double posPy = (m_centrePositionY + orbit.m_centrePositionY) / 2;
	double posPz = (m_centrePositionZ + orbit.m_centrePositionZ) / 2;
	double posAB2 = std::pow(orbit.m_centrePositionX - m_centrePositionX , 2.0)
	                + std::pow(orbit.m_centrePositionY - m_centrePositionY , 2.0)
	                + std::pow(orbit.m_centrePositionZ -  m_centrePositionZ, 2.0);
	double posPC2 = std::pow(nuclearX - posPx, 2) + std::pow(nuclearY - posPy, 2)
	                + std::pow(nuclearZ - posPz, 2);

	double exponetialFactor = std::exp(-1 * m_alpha * orbit.m_alpha * posAB2 / gamma);
	double normFactor = m_normaliseConstant * orbit.m_normaliseConstant;

	double sum(0);
	int maxK = m_k + orbit.m_k;
	int maxM = m_m + orbit.m_m;
	int maxN = m_n + orbit.m_n;
	for (int k = 0; k <= maxK; k++)
	{
		int maxR = k / 2;
		for (int r = 0; r <= maxR ; r++)
		{
			int maxI = (k - 2 * r) / 2;
			for (int i = 0; i <= maxI; i++)
			{
				double Ax = NuclearFunction(k, r, i, m_k, orbit.m_k, gamma, m_centrePositionX,
				                            orbit.m_centrePositionX, nuclearX);
				for (int m = 0; m <= maxM; m++)
				{
					int maxS = m / 2;
					for (int s = 0; s <= maxS; s++)
					{
						int maxJ = (m - 2 * s) / 2;
						for (int j = 0; j <= maxJ; j++)
						{
							double Ay = NuclearFunction(m, s, j, m_m, orbit.m_m, gamma,
							                            m_centrePositionY, orbit.m_centrePositionY,
							                            nuclearY);
							for (int n = 0; n <= maxN; n++)
							{
								int maxT = n / 2;
								for (int t = 0; t <= maxT; t++)
								{
									int maxL = (n - 2 * t) / 2;
									for (int  l = 0; l <= maxL; l++)
									{
										double Az = NuclearFunction(n, t, l,
										                            m_n, orbit.m_n, gamma,
										                            m_centrePositionZ, orbit.m_centrePositionZ,
										                            nuclearZ);
										int boysIndex = (k + m + n) - (2 * (r + s + t)) - (i + j + l);
										sum += Ax * Ay * Az * Functions::BoysFunction(boysIndex, gamma * posPC2);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return (-2 * Constants::pi / gamma) * exponetialFactor * sum * normFactor;
}

void GaussianOrbital::CalculateDataCartesian(const Vector<double> &xAxis,
        const Vector<double> &yAxis,
        const Vector<double> &zAxis)
{
	// Set the axis data;
	SetXAxis(xAxis);
	SetYAxis(yAxis);
	SetZAxis(zAxis);
	SetData(xAxis, yAxis, zAxis);

	// Set the data array
	for (int i = 0; i < xAxis.Length(); i++)
	{
		for (int j = 0; j < yAxis.Length(); j++)
		{
			for (int k = 0; k < zAxis.Length(); k++)
			{
				m_data[i][j][k] =  m_normaliseConstant * std::pow(xAxis(i), m_k)
				                   * std::pow(yAxis(j), m_m) * std::pow(zAxis(k), m_n)
				                   * std::exp(-1 * m_alpha * (std::pow(xAxis(i), 2.0)
				                              + std::pow(yAxis(j), 2.0) + std::pow(zAxis(k), 2.0)));
			}
		}
	}
}

void GaussianOrbital::CalculateDataSpherical(const Vector<double> &rAxis,
        const Vector<double> &thetaAxis,
        const Vector<double> &phiAxis)
{
	SetXAxis(rAxis);
	SetYAxis(thetaAxis);
	SetZAxis(phiAxis);
	SetData(rAxis, thetaAxis, phiAxis);

	// Set the data array
	for (int i = 0; i < rAxis.Length(); i++)
	{
		for (int j = 0; j < thetaAxis.Length(); j++)
		{
			for (int k = 0; k < phiAxis.Length(); k++)
			{
				m_data[i][j][k] =  m_normaliseConstant * std::pow(rAxis(i), (m_k + m_n + m_m))
				                   * std::pow(std::sin(thetaAxis(j)), m_k + m_m)
				                   * std::pow(std::cos(thetaAxis(j)), m_n)
				                   * std::pow(std::cos(phiAxis(k)), m_k)
				                   * std::pow(std::sin(phiAxis(k)), m_m)
				                   * std::exp(-1 * m_alpha * std::pow(rAxis(i), 2.0));
			}
		}
	}
}

void GaussianOrbital::Normalise()
{
	m_normaliseConstant = std::sqrt(1.0 / Overlap(*this));
}

double GaussianOrbital::OverlapFunction(int l1, int l2, double gamma, double posA, double posB)
{
	int maxSum = (l1 + l2) / 2;
	double sum(0);
	for (int i = 0; i <= maxSum; i++)
	{
		sum += std::sqrt(Constants::pi / gamma) * GaussianProduct(2 * i, l1, l2, posA, posB)
		       * Functions::SemiFactorial(2 * i - 1) / std::pow(2 * gamma, i);
	}
	return sum;
}

double GaussianOrbital::NuclearFunction(int l, int r, int i, int l1, int l2, double gamma,
                                        double posA, double posB, double posC)
{
	double epsilon = 1 / (4 * gamma);
	double posP = (posA + posB) / 2;
	double result = std::pow(-1, l) * GaussianProduct(l, l1, l2, posA - posP, posB - posP) * std::pow(-1, i)
	                * Functions::Factorial(l) * std::pow(posC - posP, l - 2 * r - 2 * i)
	                * std::pow(epsilon, r + i) / (Functions::Factorial(r) * Functions::Factorial(i)
	                        * Functions::Factorial(l - (2 * r) - (2 * i)));
	return result;
}

double GaussianOrbital::GaussianProduct(int k, int l1, int l2, double pos1, double pos2)
{
	int maxSum = std::min(k, 2 * l1 - k);
	int minSum = std::max(-k, k - 2 * l2);
	double sum(0);
	for (int i = minSum; i <= maxSum; i = i + 2)
	{
		int q = (k + i) / 2;
		int p = (k - i) / 2;
		sum += Functions::BinomialCoefficient(l1, q) * Functions::BinomialCoefficient(l2, p)
		       * std::pow(pos1, l1 - q) * std::pow(pos2, l2 - p);
	}
	return sum;
}

double GaussianOrbital::GaussianProduct2(int k, int l1, int l2, double pos1, double pos2)
{
	double sum(0);
	for (int i = 0; i <= l1; i++)
	{
		for (int j = k - i; j <= l2; j++)
		{
			sum += Functions::BinomialCoefficient(l1, i) * Functions::BinomialCoefficient(l2, j)
			       * std::pow(pos1, l1 - i) * std::pow(pos2, l2 - j);
		}
	}
	return sum;
}

GaussianOrbital GaussianOrbital::ChangeK(int k) const
{
	GaussianOrbital orbit(m_k + k, m_m, m_n, m_alpha, m_centrePositionX, m_centrePositionY,
	                      m_centrePositionZ);
	return orbit;
}

GaussianOrbital GaussianOrbital::ChangeM(int m) const
{
	GaussianOrbital orbit(m_k, m_m + m, m_n, m_alpha, m_centrePositionX, m_centrePositionY,
	                      m_centrePositionZ);
	return orbit;
}

GaussianOrbital GaussianOrbital::ChangeN(int n) const
{
	GaussianOrbital orbit(m_k, m_m, m_n + n, m_alpha,
	                      m_centrePositionX, m_centrePositionY, m_centrePositionZ);
	return orbit;
}