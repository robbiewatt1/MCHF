#include <cmath>
#include <algorithm>

#include "GaussianOrbital.hh"
#include "Functions.hh"
#include "Constants.hh"
#include "Array3D.hh"

GaussianOrbital::GaussianOrbital():
	m_k(0), m_m(0), m_n(0), m_alpha(1), m_orbitPosition(3)
{
	Normalise();
}

GaussianOrbital::GaussianOrbital(int k, int m, int n, double alpha, Vector<double> orbitPosition):
	m_k(k), m_m(m), m_n(n), m_alpha(alpha)
{
	m_orbitPosition = orbitPosition;
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
	return m_orbitPosition[0];
}

double GaussianOrbital::GetNormaliseConstant()
{
	return m_normaliseConstant;
}

double GaussianOrbital::Overlap(const GaussianOrbital &orbit)
{
	double gamma  = m_alpha + orbit.m_alpha;
	Vector<double> positionP(3);
	positionP[0] = ((m_orbitPosition[0] * m_alpha) + (orbit.m_orbitPosition[0] * orbit.m_alpha))
	               / gamma;
	positionP[1] = ((m_orbitPosition[1] * m_alpha) + (orbit.m_orbitPosition[1] * orbit.m_alpha))
	               / gamma;
	positionP[2] = ((m_orbitPosition[2] * m_alpha) + (orbit.m_orbitPosition[2] * orbit.m_alpha))
	               / gamma;
	double posAB2 = std::pow(orbit.m_orbitPosition[0] - m_orbitPosition[0], 2.0)
	                + std::pow(orbit.m_orbitPosition[1] - m_orbitPosition[1], 2.0)
	                + std::pow(orbit.m_orbitPosition[2] - m_orbitPosition[2], 2.0);

	double overlapX = OverlapFunction(m_k, orbit.m_k, gamma, m_orbitPosition[0] - positionP[0],
	                                  orbit.m_orbitPosition[0] - positionP[0]);
	double overlapY = OverlapFunction(m_m, orbit.m_m, gamma, m_orbitPosition[1] - positionP[1],
	                                  orbit.m_orbitPosition[1] - positionP[1]);
	double overlapZ = OverlapFunction(m_n, orbit.m_n, gamma, m_orbitPosition[2] - positionP[2],
	                                  orbit.m_orbitPosition[2] - positionP[2]);
	double overlap = std::exp(-1 * m_alpha * orbit.m_alpha * posAB2 / gamma)
	                 * overlapX * overlapY * overlapZ;
	return overlap * m_normaliseConstant * orbit.m_normaliseConstant;
}

double GaussianOrbital::KineticOverlap(const GaussianOrbital &orbit)
{
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

	double kinOverlap = (orbit.m_alpha * (2.0 * (orbit.m_k + orbit.m_m + orbit.m_n) + 3.0) * overlap)
	                    - (2.0 * std::pow(orbit.m_alpha, 2.0) * (plus2k + plus2m + plus2n))
	                    - (0.5 * (orbit.m_k * (orbit.m_k - 1.0) * minus2K
	                              + orbit.m_m * (orbit.m_m - 1.0) * minus2M
	                              + orbit.m_n * (orbit.m_n - 1.0) * minus2N));
	return kinOverlap;
}

double GaussianOrbital::NuclearOverlap(const GaussianOrbital &orbit, int nuclearCharge,
                                       Vector<double> nuclearPosition)
{
	double gamma  = m_alpha + orbit.m_alpha;
	Vector<double> positionP(3);
	positionP[0] = ((m_orbitPosition[0] * m_alpha) + (orbit.m_orbitPosition[0] * orbit.m_alpha))
	               / gamma;
	positionP[1] = ((m_orbitPosition[1] * m_alpha) + (orbit.m_orbitPosition[1] * orbit.m_alpha))
	               / gamma;
	positionP[2] = ((m_orbitPosition[2] * m_alpha) + (orbit.m_orbitPosition[2] * orbit.m_alpha))
	               / gamma;
	double posAB2 = std::pow(orbit.m_orbitPosition[0] - m_orbitPosition[0], 2.0)
	                + std::pow(orbit.m_orbitPosition[1] - m_orbitPosition[1], 2.0)
	                + std::pow(orbit.m_orbitPosition[2] - m_orbitPosition[2], 2.0);
	double posPC2 = std::pow(nuclearPosition[0] - positionP[0], 2.0)
	                + std::pow(nuclearPosition[1] - positionP[1], 2.0)
	                + std::pow(nuclearPosition[2] - positionP[2], 2.0);

	double exponetialFactor = std::exp(-1 * m_alpha * orbit.m_alpha * posAB2 / gamma);

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
				double Ax = NuclearFunction(k, r, i, m_k, orbit.m_k, m_alpha, orbit.m_alpha,
				                            m_orbitPosition[0], orbit.m_orbitPosition[0],
				                            nuclearPosition[0]);
				for (int m = 0; m <= maxM; m++)
				{
					int maxS = m / 2;
					for (int s = 0; s <= maxS; s++)
					{
						int maxJ = (m - 2 * s) / 2;
						for (int j = 0; j <= maxJ; j++)
						{
							double Ay = NuclearFunction(m, s, j, m_m, orbit.m_m, m_alpha,
							                            orbit.m_alpha, m_orbitPosition[1],
							                            orbit.m_orbitPosition[1],
							                            nuclearPosition[1]);
							for (int n = 0; n <= maxN; n++)
							{
								int maxT = n / 2;
								for (int t = 0; t <= maxT; t++)
								{
									int maxL = (n - 2 * t) / 2;
									for (int  l = 0; l <= maxL; l++)
									{
										double Az = NuclearFunction(n, t, l,
										                            m_n, orbit.m_n, m_alpha,
										                            orbit.m_alpha,
										                            m_orbitPosition[2],
										                            orbit.m_orbitPosition[2],
										                            nuclearPosition[2]);
										int boysIndex = (k + m + n) - (2 * (r + s + t))
										                - (i + j + l);
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
	return (-2.0 * Constants::pi / gamma) * exponetialFactor * sum * nuclearCharge 
		   * m_normaliseConstant * orbit.m_normaliseConstant;
}

Array3D<double> GaussianOrbital::CalculateDataCartesian(const Vector<double> &xAxis,
        const Vector<double> &yAxis,
        const Vector<double> &zAxis)
{
	SetXAxis(xAxis);
	SetYAxis(yAxis);
	SetZAxis(zAxis);
	Array3D<double> data(xAxis.Length(), yAxis.Length(), zAxis.Length());

	for (int i = 0; i < xAxis.Length(); i++)
	{
		for (int j = 0; j < yAxis.Length(); j++)
		{
			for (int k = 0; k < zAxis.Length(); k++)
			{
				data[i][j][k] =  m_normaliseConstant * std::pow(xAxis(i) - m_orbitPosition[0], m_k)
				                   * std::pow(yAxis(j) - m_orbitPosition[1], m_m) 
				                   * std::pow(zAxis(k) - m_orbitPosition[2], m_n)
				                   * std::exp(-1 * m_alpha * (std::pow(xAxis(i) - m_orbitPosition[0], 2.0)
				                              + std::pow(yAxis(j) - m_orbitPosition[1], 2.0) 
				                              + std::pow(zAxis(k) - m_orbitPosition[2], 2.0)));
			}
		}
	}
	return data;
}

Array3D<double> GaussianOrbital::CalculateDataSpherical(const Vector<double> &rAxis,
        const Vector<double> &thetaAxis,
        const Vector<double> &phiAxis)
{
	SetXAxis(rAxis);
	SetYAxis(thetaAxis);
	SetZAxis(phiAxis);
	Array3D<double> data(rAxis.Length(), thetaAxis.Length(), phiAxis.Length());

	// Set the data array
	for (int i = 0; i < rAxis.Length(); i++)
	{
		for (int j = 0; j < thetaAxis.Length(); j++)
		{
			for (int k = 0; k < phiAxis.Length(); k++)
			{
				data[i][j][k] =  m_normaliseConstant * std::pow(rAxis(i), (m_k + m_n + m_m))
				                   * std::pow(std::sin(thetaAxis(j)), m_k + m_m)
				                   * std::pow(std::cos(thetaAxis(j)), m_n)
				                   * std::pow(std::cos(phiAxis(k)), m_k)
				                   * std::pow(std::sin(phiAxis(k)), m_m)
				                   * std::exp(-1.0 * m_alpha * std::pow(rAxis(i), 2.0));
			}
		}
	}
	return data;
}

void GaussianOrbital::Normalise()
{
	double gamma  = 2.0 * m_alpha;
	double overlapX = OverlapFunction(m_k, m_k, gamma, 0, 0);
	double overlapY = OverlapFunction(m_m, m_m, gamma, 0, 0);
	double overlapZ = OverlapFunction(m_n, m_n, gamma, 0, 0);
	double overlap = overlapX * overlapY * overlapZ;
	m_normaliseConstant = std::sqrt(1.0 / overlap);
}

double GaussianOrbital::OverlapFunction(int l1, int l2, double gamma, double posA, double posB)
{
	int maxSum = (l1 + l2) / 2;
	double sum(0);
	for (int i = 0; i <= maxSum; i++)
	{
		sum += GaussianProduct(2 * i, l1, l2, posA, posB) * Functions::SemiFactorial(2 * i - 1)
		       / std::pow(2.0 * gamma, i);
	}
	return std::sqrt(Constants::pi / gamma) * sum;
}

double GaussianOrbital::NuclearFunction(int l, int r, int i, int l1, int l2, double alpha,
                                        double beta, double posA, double posB, double posC)
{
	double gamma = alpha + beta;
	double epsilon = 1.0 / (4.0 * gamma);
	double posP = (posA * alpha + posB * beta) / gamma;
	double result = std::pow(-1.0, l) * GaussianProduct(l, l1, l2, posA - posP, posB - posP) * std::pow(-1.0, i)
	                * Functions::Factorial(l) * std::pow(posC - posP, l - 2.0 * r - 2.0 * i)
	                * std::pow(epsilon, r + i) / (Functions::Factorial(r) * Functions::Factorial(i)
	                        * Functions::Factorial(l - (2.0 * r) - (2 * i)));
	return result;
}

double GaussianOrbital::GaussianProduct3(int k, int l1, int l2, double pos1, double pos2)
{
	int maxSum = std::min(k, 2 * l1 - k);
	int minSum = std::max(-k, k - 2 * l2);
	double sum(0);
	for (int q = minSum; q <= maxSum; q = q + 2)
	{
		int i = (k + q) / 2;
		int j = (k - q) / 2;
		sum += Functions::BinomialCoefficient(l1, i) * Functions::BinomialCoefficient(l2, j)
		       * std::pow(pos1, l1 - i) * std::pow(pos2, l2 - j);
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

double GaussianOrbital::GaussianProduct(int k, int l1, int l2, double pos1, double pos2)
{
	double sum(0);
	int minSum = std::max(0, k - l2);
	int maxSum = std::min(k, l1);
	for (int i = minSum; i <= maxSum ; i++)
	{
		sum += Functions::BinomialCoefficient(l1, i) * Functions::BinomialCoefficient(l2, k - i)
		       * std::pow(pos1, l1 - i) * std::pow(pos2, l2 + i - k);
	}
	return sum;
}

GaussianOrbital GaussianOrbital::ChangeK(int k) const
{
	GaussianOrbital orbit(m_k + k, m_m, m_n, m_alpha, m_orbitPosition);
	return orbit;
}

GaussianOrbital GaussianOrbital::ChangeM(int m) const
{
	GaussianOrbital orbit(m_k, m_m + m, m_n, m_alpha, m_orbitPosition);
	return orbit;
}

GaussianOrbital GaussianOrbital::ChangeN(int n) const
{
	GaussianOrbital orbit(m_k, m_m, m_n + n, m_alpha, m_orbitPosition);
	return orbit;
}
