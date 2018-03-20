#include <cmath>

#include "Molecule.hh"
#include "LinearAlgebra.hh"

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
	for (int i = 0; i < m_nuclearCharges.Length(); i++)
	{
		for (int j = i + 1; j < m_nuclearCharges.Length(); j++ )
		{
			potential += CalculatePotential(m_nuclearCharges[i], m_nuclearCharges[j],
			                                m_nuclearPositions[i], m_nuclearPositions[j]);
		}
	}

	int ionIndex(0);
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			energyMaxtrix[i][j] = (m_basisSet[i].KineticOverlap(m_basisSet[j])
			                       + m_basisSet[i].NuclearOverlap(m_basisSet[j],
			                               m_nuclearCharges[ionIndex],
			                               m_nuclearPositions[ionIndex][0],
			                               m_nuclearPositions[ionIndex][1],
			                               m_nuclearPositions[ionIndex][2]))
			                      * m_basisSet[i].GetNormaliseConstant()
			                      * m_basisSet[j].GetNormaliseConstant()
			                      + potential;

			overlapMatrix[i][j] = m_basisSet[i].Overlap(m_basisSet[j])
			                      * m_basisSet[i].GetNormaliseConstant()
			                      * m_basisSet[j].GetNormaliseConstant();
			ionIndex = (ionIndex + 1) % m_nuclearCharges.Length();
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
					for (int i = 0; i < m_nuclearPositions.Length() ; ++i)
					{
						STOnGOrbit orbital = STOnGOrbit("./OrbitalData/STO3Test", k, m, n,
						                                m_nuclearPositions[i][0],
						                                m_nuclearPositions[i][1],
						                                m_nuclearPositions[i][2]);
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
	                        + std::pow(ionLocation1[0] - ionLocation2[0], 2.0));
	return z1 * z2 / ionR;
}