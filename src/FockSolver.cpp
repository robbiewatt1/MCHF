#include <cmath>

#include "FockSolver.hh"
#include "LinearAlgebra.hh"

FockSolver::FockSolver(const Vector<double> &nuclearCharges,
					   const Vector<Vector<double>> &nuclearPositions,
					   const Vector<STOnGOrbit> &basisSet, int nElectrons,
					   const BoysFunction &boyFn):
m_nElectrons(nElectrons), m_basisSet(basisSet), m_nuclearPositions(nuclearPositions),
m_nuclearCharges(nuclearCharges)
{
	m_boyFn = boyFn;
}

FockSolver::~FockSolver()
{
}


void FockSolver::Solve()
{
	// First we need to set the m_coeff matrix. Not sure of the best initail guess so will start 
	// with everything 1
	m_coeffC = Matrix<double>(m_basisSet.Length(), m_nElectrons);
	InitialGuess(m_coeffC);
	
	// Solve the system 
	OneElectronSolver();
	CoulombSolver();
	ExchangeSolver();

	// Form the fock matrix
	Matrix<double> fockMaxtrix = m_oneElectronEnergy + 2.0 * m_coulombEnergy - m_exchangeEnergy;

	// Transform into the othogonal basis
//	Vector<double> overlapValues;
//	Matrix<double> overlapVector;
//	LinearAlgebra::EigenSolver(overlapMatrix, overlapVector, overlapValues);


}

void FockSolver::OneElectronSolver()
{
	m_basisOverlap = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	m_oneElectronEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());

	// calculate the ion repulsion potenital
	double ionPotential(0);
	#pragma omp parallel for
	for (int i = 0; i < m_nuclearCharges.Length(); i++)
	{
		for (int j = i + 1; j < m_nuclearCharges.Length(); j++)
		{
			ionPotential += IonPotential(m_nuclearCharges[i], m_nuclearCharges[j],
									  m_nuclearPositions[i], m_nuclearPositions[j]);
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = i; j < m_basisSet.Length(); j++)
		{
			double nuclearPotential(0);
			for (int k = 0; k < m_nuclearCharges.Length(); k++)
			{
				nuclearPotential += m_basisSet[i].NuclearOverlap(m_basisSet[j],
				                    m_nuclearCharges[k],
									m_nuclearPositions[k], m_boyFn);
			}
			double kineticEnergy = m_basisSet[i].KineticOverlap(m_basisSet[j]);
			m_basisOverlap[i][j] = m_basisSet[i].Overlap(m_basisSet[j]);
			m_oneElectronEnergy[i][j] = kineticEnergy + nuclearPotential + ionPotential * m_basisOverlap[i][j];
			m_basisOverlap[j][i] = m_basisOverlap[i][j];
			m_oneElectronEnergy[j][i] = m_oneElectronEnergy[i][j];
		}
	}
}

void FockSolver::CoulombSolver()
{
	m_coulombEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	Matrix<double> weights(m_basisSet.Length(), m_basisSet.Length());

	// Check this is in the right order
	weights = m_coeffC * LinearAlgebra::Transpose(m_coeffC);
	#pragma omp parallel for
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			double coulombEnergy(0);
			for (int k = 0; k < m_basisSet.Length(); k++)
			{
				for (int l = 0; l < m_basisSet.Length(); +l++)
				{
					coulombEnergy += weights[k][l] * m_basisSet[i].ElectronRepulsion(m_basisSet[j],
																					 m_basisSet[k],
																					 m_basisSet[l],
																					 m_boyFn);
				}
			}
			m_coulombEnergy[i][j] = coulombEnergy;
			// Think this should be symetirc but will need to check
			std::cout << i << " " << j << std::endl;
		}
	}
}

void FockSolver::ExchangeSolver()
{
	m_exchangeEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	Matrix<double> weights(m_basisSet.Length(), m_basisSet.Length());

	// Check this is in the right order
	weights = m_coeffC * LinearAlgebra::Transpose(m_coeffC);
	#pragma omp parallel for
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			double exchangeEnergy(0);
			for (int k = 0; k < m_basisSet.Length(); k++)
			{
				for (int l = 0; l < m_basisSet.Length(); +l++)
				{
					exchangeEnergy += weights[k][l] * m_basisSet[i].ElectronRepulsion(m_basisSet[l],
																					 m_basisSet[k],
																					 m_basisSet[j],
																					 m_boyFn);
				}
			}
			m_exchangeEnergy[i][j] = exchangeEnergy;
			// Think this should be symetirc but will need to check
		}
	}
}

double FockSolver::IonPotential(int z1, int z2, Vector<double> ionLocation1,
								Vector<double> ionLocation2)
{
	double ionR = std::sqrt(std::pow(ionLocation1[0] - ionLocation2[0], 2.0)
	                        + std::pow(ionLocation1[1] - ionLocation2[1], 2.0)
	                        + std::pow(ionLocation1[2] - ionLocation2[2], 2.0));
	return z1 * z2 / ionR;
}

void FockSolver::InitialGuess(Matrix<double> &coeffC)
{
	for (int i = 0; i < coeffC.GetRows(); ++i)
	{
		for (int j = 0; j < coeffC.GetColumns(); ++j)
		{
			if (i == j)
			{
				coeffC[i][j] = 1.0;
			} else
			{
				coeffC[i][j] = 0.0;
			}
		}
	}
}
