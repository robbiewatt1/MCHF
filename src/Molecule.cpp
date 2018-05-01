#include <cmath>
#include <omp.h>
#include "Molecule.hh"
#include "LinearAlgebra.hh"
#include "STOnGOrbit.hh"

Molecule::Molecule()
{
}

Molecule::Molecule(const Vector<double> &nuclearCharges,
                   const Vector<Vector<double>> &nuclearPositions, int maxL):
	m_nuclearCharges(nuclearCharges), m_nuclearPositions(nuclearPositions), m_maxL(maxL)
{
	SetBasisSet();
}

Molecule::~Molecule()
{
}

Vector<double> Molecule::GetEnergyLevels()
{
	return m_energyLevels;
}

Matrix<double> Molecule::GetBasisCoefficients()
{
	return m_basisSetCoefficients;
}

void Molecule::CalculateEnergy()
{
	Matrix<double> energyMaxtrix(m_basisSet.Length(), m_basisSet.Length());
	Matrix<double> overlapMatrix(m_basisSet.Length(), m_basisSet.Length());

	double potential(0);
	#pragma omp parallel for
	for (int i = 0; i < m_nuclearCharges.Length(); i++)
	{
		for (int j = i + 1; j < m_nuclearCharges.Length(); j++)
		{
			potential += CalculatePotential(m_nuclearCharges[i], m_nuclearCharges[j],
			                                m_nuclearPositions[i], m_nuclearPositions[j]);
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = i; j < m_basisSet.Length(); j++)
		{
			//Check omp is working
			//std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
			// loop over ion sights to get nuclear attraction integral
			double nuclearPotential(0);
			for (int k = 0; k < m_nuclearCharges.Length(); k++)
			{
				nuclearPotential += m_basisSet[i].NuclearOverlap(m_basisSet[j],
				                    m_nuclearCharges[k],
				                    m_nuclearPositions[k])
				                    * m_basisSet[j].GetNormaliseConstant()
				                    * m_basisSet[i].GetNormaliseConstant();
			}
			double kineticEnergy = m_basisSet[i].KineticOverlap(m_basisSet[j])
			                       * m_basisSet[i].GetNormaliseConstant()
			                       * m_basisSet[j].GetNormaliseConstant();

			overlapMatrix[i][j] = m_basisSet[i].Overlap(m_basisSet[j])
			                      * m_basisSet[i].GetNormaliseConstant()
			                      * m_basisSet[j].GetNormaliseConstant();
			energyMaxtrix[i][j] = kineticEnergy + (potential * overlapMatrix[i][j]) + nuclearPotential;
			
			// Matrix is symetric
			overlapMatrix[j][i] = overlapMatrix[i][j];
			energyMaxtrix[j][i] = energyMaxtrix[i][j];
		}
	}
	m_energyLevels = Vector<double>(m_basisSet.Length());
	m_basisSetCoefficients = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	LinearAlgebra::GeneralisedEigenSolver(energyMaxtrix, overlapMatrix, m_basisSetCoefficients,
	                                      m_energyLevels);
}

void Molecule::SetBasisSet()
{
	// Loop over k, m and n such that the sum is less than the maximum L
	for (int k = 0; k <= m_maxL; k++)
	{
		for (int m = 0; m <= m_maxL; m++)
		{
			for (int n = 0; n <= m_maxL; n++)
			{
				if (n + m + k <= m_maxL)
				{
					// Now loop over all ion sights
					for (int i = 0; i < m_nuclearPositions.Length(); i++)
					{
						STOnGOrbit orbital;
						if (m_nuclearCharges[i] == 1)
						{
							orbital = STOnGOrbit("./OrbitalData/STO6test", k, m, n,
											     m_nuclearPositions[i]);
						} else if(m_nuclearCharges[i] == 2)
						{
							orbital = STOnGOrbit("./OrbitalData/STO6He", k, m, n,
												 m_nuclearPositions[i]);
						}
						m_basisSet.Append(orbital);
					}
				}
			}
		}
	}
}

double Molecule::CalculatePotential(int z1, int z2, Vector<double> ionLocation1,
                                    Vector<double> ionLocation2)
{
	double ionR = std::sqrt(std::pow(ionLocation1[0] - ionLocation2[0], 2.0)
	                        + std::pow(ionLocation1[1] - ionLocation2[1], 2.0)
	                        + std::pow(ionLocation1[2] - ionLocation2[2], 2.0));
	return z1 * z2 / ionR;
}
