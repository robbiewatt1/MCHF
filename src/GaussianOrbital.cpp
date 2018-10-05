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

int GaussianOrbital::GetK() const
{
	return m_k;
}

int GaussianOrbital::GetM() const
{
	return m_m;
}

int GaussianOrbital::GetN() const
{
	return m_n;
}

double GaussianOrbital::GetAlpha() const
{
	return m_alpha;
}

double GaussianOrbital::GetPositionX() const
{
	return m_orbitPosition[0];
}

double GaussianOrbital::GetNormaliseConstant() const
{
	return m_normaliseConstant;
}

double GaussianOrbital::Overlap(const GaussianOrbital &orbit) const
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
	return overlap;
}

double GaussianOrbital::KineticOverlap(const GaussianOrbital &orbit) const
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
                                       const Vector<double> &nuclearPosition, const BoysFunction &boyFn) const
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

	double exponetialFactor = std::exp(-1.0 * m_alpha * orbit.m_alpha * posAB2 / gamma);

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
										int boysIndex = (k + m + n) - ( 2 * (r + s + t))
										                - (i + j + l);
										sum += Ax * Ay * Az * boyFn.Interpolate(boysIndex, gamma * posPC2);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return (-2.0 * Constants::pi / gamma) * exponetialFactor * sum * nuclearCharge;
}

double GaussianOrbital::ElectronRepulsion(const GaussianOrbital &orbit1, const GaussianOrbital &orbit2, 
										  const GaussianOrbital &orbit3, const BoysFunction &boyFn) const
{
	double gamma1 = m_alpha + orbit1.m_alpha;
	double gamma2 = orbit2.m_alpha + orbit3.m_alpha;
	double delta  = (1.0 / (4.0 * gamma1)) + (1.0 / (4.0 * gamma2));
	Vector<double> positionP(3);
	Vector<double> positionQ(3);
	positionP[0] = ((m_orbitPosition[0] * m_alpha) + (orbit1.m_orbitPosition[0] * orbit1.m_alpha))
	               / gamma1;
	positionP[1] = ((m_orbitPosition[1] * m_alpha) + (orbit1.m_orbitPosition[1] * orbit1.m_alpha))
	               / gamma1;
	positionP[2] = ((m_orbitPosition[2] * m_alpha) + (orbit1.m_orbitPosition[2] * orbit1.m_alpha))
	               / gamma1;
	positionQ[0] = ((orbit2.m_orbitPosition[0] * orbit2.m_alpha) + (orbit3.m_orbitPosition[0] * orbit3.m_alpha))
	               / gamma2;
	positionQ[1] = ((orbit2.m_orbitPosition[1] * orbit2.m_alpha) + (orbit3.m_orbitPosition[1] * orbit3.m_alpha))
	               / gamma2;
	positionQ[2] = ((orbit2.m_orbitPosition[2] * orbit2.m_alpha) + (orbit3.m_orbitPosition[2] * orbit3.m_alpha))
	               / gamma2;

	double posAB2 = std::pow(orbit1.m_orbitPosition[0] - m_orbitPosition[0], 2.0)
	              + std::pow(orbit1.m_orbitPosition[1] - m_orbitPosition[1], 2.0)
	              + std::pow(orbit1.m_orbitPosition[2] - m_orbitPosition[2], 2.0);
	double posCD2 = std::pow(orbit3.m_orbitPosition[0] - orbit2.m_orbitPosition[0], 2.0)
	              + std::pow(orbit3.m_orbitPosition[1] - orbit2.m_orbitPosition[1], 2.0)
	              + std::pow(orbit3.m_orbitPosition[2] - orbit2.m_orbitPosition[2], 2.0);
	double posQP2 = std::pow(positionP[0] - positionQ[0], 2.0)
				  + std::pow(positionP[1] - positionQ[1], 2.0)
				  + std::pow(positionP[2] - positionQ[2], 2.0);

	double gaussProduct = 2.0 * std::pow(Constants::pi, 2.0) / (gamma1 * gamma2) * std::sqrt(Constants::pi / (gamma1 + gamma2))
						* std::exp(-1.0 * (m_alpha * orbit1.m_alpha * posAB2 / gamma1) 
								+ (-1.0 * (orbit2.m_alpha * orbit3.m_alpha * posCD2 / gamma2)));
	if (gaussProduct < 1e-8)
	{
		return 0;
	}
	double sum(0);
	int maxK   = m_k + orbit1.m_k;
	int maxK_p = orbit2.m_k + orbit3.m_k;
	int maxM   = m_m + orbit1.m_m;
	int maxM_P = orbit2.m_m + orbit3.m_m;
	int maxN   = m_n + orbit1.m_n;
	int maxN_p = orbit2.m_n + orbit3.m_n;
	for (int k = 0; k <= maxK; k++)
	{
		int maxR = k / 2;
		for (int r = 0; r <= maxR; r++)
		{
			for (int k_p = 0; k_p <= maxK_p; k_p++)
			{
				int maxR_p = k_p / 2;
				for (int r_p = 0; r_p <= maxR_p; r_p++)
				{
					int maxI = (k + k_p - 2 * (r + r_p)) / 2;
					for (int i = 0; i <= maxI; i++)
					{			
						double Bx = ElectronFunction(k, k_p, r, r_p, i, m_k, orbit1.m_k, m_orbitPosition[0], orbit1.m_orbitPosition[0], positionP[0], gamma1, 
													 orbit2.m_k, orbit3.m_k, orbit2.m_orbitPosition[0], orbit3.m_orbitPosition[0], positionQ[0], gamma2);
						if (Bx == 0)
						{
							continue;
						}
						for (int m = 0; m <= maxM; m++)
						{
							int maxS = m / 2;
							for (int s = 0; s <= maxS; s++)
							{
								for (int m_p = 0; m_p <= maxM_P; m_p++)
								{
									int maxS_p = m_p / 2;
									for (int s_p = 0; s_p <= maxS_p; s_p++)
									{
										int maxJ = (m + m_p - 2 * (s + s_p)) / 2;
										for (int j = 0; j <= maxJ; j++)
										{
											double By = ElectronFunction(m, m_p, s, s_p, j, m_m, orbit1.m_m, m_orbitPosition[1], orbit1.m_orbitPosition[1], positionP[1], gamma1,
													 					 orbit2.m_m, orbit3.m_m, orbit2.m_orbitPosition[1], orbit3.m_orbitPosition[1], positionQ[1], gamma2);
											if (By == 0)
											{
												continue;
											}
											for (int n = 0; n <= maxN; n++)
											{
												int maxT = n / 2;
												for (int t = 0; t <= maxT; t++)
												{
													for (int n_p = 0; n_p <= maxN_p; n_p++)
													{
														int maxT_p = n_p / 2;
														for (int t_p = 0; t_p <= maxT_p; t_p++)
														{
															int maxL = (n + n_p - 2 * (t + t_p)) / 2;
															for (int l = 0; l <= maxL; l++)
															{
																double Bz = ElectronFunction(n, n_p, t, t_p, l, m_n, orbit1.m_n, m_orbitPosition[2], orbit1.m_orbitPosition[2], positionP[2],
																			gamma1, orbit2.m_n, orbit3.m_n, orbit2.m_orbitPosition[2], orbit3.m_orbitPosition[2], positionQ[2], gamma2);
																if (Bz == 0)
																{
																	continue;
																}
																int boysIndex = k + k_p + m + m_p + n + n_p - 2 * (r + r_p + s + s_p + t + t_p) - (i + j + l);
																sum += Bx * By * Bz * boyFn.Interpolate(boysIndex, posQP2 / (4.0 * delta));
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return gaussProduct * sum;
}

Vector<double> GaussianOrbital::MatrixElement(const GaussianOrbital &orbit) const
{
	
	Vector<double> matrixElements(3);
	double plusK  = this->Overlap(orbit.ChangeK(1));
	double plusM  = this->Overlap(orbit.ChangeM(1));
	double plusN  = this->Overlap(orbit.ChangeN(1));
	double minusK(0), minusM(0), minusN(0);
	if (orbit.m_k > 0)
	{
		minusK = this->Overlap(orbit.ChangeK(-1));
	}
	if (orbit.m_m > 0)
	{
		minusM = this->Overlap(orbit.ChangeM(-1));
	}
	if (orbit.m_n > 0)
	{
		minusN = this->Overlap(orbit.ChangeN(-1));
	}
	matrixElements[0] = orbit.m_k * minusK - 2 * orbit.m_alpha * plusK;
	matrixElements[1] = orbit.m_m * minusM - 2 * orbit.m_alpha * plusM;
	matrixElements[2] = orbit.m_n * minusN - 2 * orbit.m_alpha * plusN;
	return matrixElements;
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

double GaussianOrbital::OverlapFunction(int l1, int l2, double gamma, double posA, double posB) const
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
                                        double beta, double posA, double posB, double posC) const
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

double GaussianOrbital::ElectronFunction(int l, int l_p, int r, int r_p, int i, int l1, int l2, double posA,
										 double posB, double posP, double gamma1, int l3, int l4, double posC, 
										 double posD, double posQ, double gamma2) const
{
	double delta = (1.0 / (4.0 * gamma1)) + (1.0 / (4.0 * gamma2));
	double posR  = posP - posQ;
	double posPA = posA - posP;
	double posPB = posB - posP;
	double posQC = posC - posQ;
	double posQD = posD - posQ;
	double result = std::pow(-1.0, l_p) * ThetaFn(l, l1, l2, posPA, posPB, r, gamma1) 
				  * ThetaFn(l_p, l3, l4, posQC, posQD, r_p, gamma2) * std::pow(-1.0, i) 
				  * std::pow(2.0 * delta, 2.0 * (r + r_p)) * Functions::Factorial(l + l_p - 2.0 * (r + r_p))
				  * std::pow(delta, i) * std::pow(posR, l + l_p - 2.0 * (r + r_p + i))
				  / (std::pow(4.0 * delta, l + l_p) * Functions::Factorial(i) * 
				  	Functions::Factorial(l + l_p - 2.0 * (r + r_p + i)));
	return result;
}

double GaussianOrbital::ThetaFn(int l, int l1, int l2, double a, double b, double c, double d) const
{
	return GaussianProduct(l, l1, l2, a, b) * Functions::Factorial(l) * std::pow(d, c - l) 
		   / (Functions::Factorial(c) * Functions::Factorial(l - 2.0 * c));
}

double GaussianOrbital::GaussianProduct(int k, int l1, int l2, double pos1, double pos2) const
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
