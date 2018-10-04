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
	InitialGuess();

	// Solve for the one electron matrix. This only needs to be done once
	OneElectronSolver();
	ElectronRepulsionSolver();
	for (int i = 0; i < 2; i++)
	{
		TwoElectronSolver();
		Matrix<double> fockMaxtrix = m_oneElectronEnergy + 2.0 * m_coulombEnergy - m_exchangeEnergy;
		//(2.0 * m_coulombEnergy - m_exchangeEnergy).Print();
		// Find matgrix used to transform to othogonal basis. Here I might need to remove bad matricies!
		Vector<double> othgTransVector;
		Matrix<double> othgTransMatrix;
		LinearAlgebra::EigenSolver(m_basisOverlap, othgTransMatrix, othgTransVector);
		Matrix<double> eigenMatrix(m_basisOverlap.GetRows(), m_basisOverlap.GetColumns());
		for (int j = 0; j < m_basisOverlap.GetRows(); j++)
		{
			for (int k = 0; k < m_basisOverlap.GetColumns(); k++)
			{
				if (k == j)
				{
					eigenMatrix[k][k] = 1.0 / std::sqrt(othgTransVector[k]);
				} else
				{
					eigenMatrix[j][k] = 0;
				}
			}
		}
		othgTransMatrix = othgTransMatrix * eigenMatrix * LinearAlgebra::Transpose(othgTransMatrix);
		fockMaxtrix = LinearAlgebra::Transpose(othgTransMatrix) * fockMaxtrix * othgTransMatrix;
		fockMaxtrix.Print();
		Matrix<double> orbitalCoeff;
		LinearAlgebra::EigenSolver(fockMaxtrix, orbitalCoeff, m_energyLevels);
		orbitalCoeff = othgTransMatrix * orbitalCoeff;
		// only occupided oribtals
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			for (int k = 0; k <  m_basisSet.Length(); k++)
			{
				double sum(0);
				for (int l = 0; l < m_nElectrons; l++)
				{
					sum += orbitalCoeff[j][l] * orbitalCoeff[k][l];
				}
				m_density[j][k] = sum;
			}
		}
		m_density.Print();
		double sum(0);
		for (int j = 0; j < m_basisSet.Length(); ++j)
		{
			for (int k = 0; k < m_basisSet.Length(); ++k)
			{
				sum += m_density[j][k] * (2.0 * m_oneElectronEnergy[j][k] + 2.0 * m_coulombEnergy[j][k] - m_exchangeEnergy[j][k]);
			}
		}
		std::cout << sum << std::endl;
		m_energyLevels.Print();
	//	m_density.Print();
	//	std::cout << LinearAlgebra::Trace(fockMaxtrix * m_density) << std::endl;
	//	m_energyLevels.Print();
	}

}

void FockSolver::OneElectronSolver()
{
	m_basisOverlap = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	m_oneElectronEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());

	// calculate the ion repulsion potenital
	double ionPotential(0);
//	#pragma omp parallel for
//	for (int i = 0; i < m_nuclearCharges.Length(); i++)
//	{
//		for (int j = i + 1; j < m_nuclearCharges.Length(); j++)
//		{
//			ionPotential += IonPotential(m_nuclearCharges[i], m_nuclearCharges[j],
//			                             m_nuclearPositions[i], m_nuclearPositions[j]);
//		}
//	}
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

void FockSolver::TwoElectronSolver()
{
	m_coulombEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	m_exchangeEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	#pragma omp parallel for	
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = i; j < m_basisSet.Length(); j++)
		{
			double exchangeEnergy(0);
			double coulombEnergy(0);
			for (int k = 0; k < m_basisSet.Length(); k++)
			{
				for (int l = 0; l < m_basisSet.Length(); l++)
				{
					exchangeEnergy += m_density[k][l] * m_eRepulsion[i][k][j][l];
					coulombEnergy  += m_density[k][l] * m_eRepulsion[i][j][k][l];
				}
			}
			m_exchangeEnergy[i][j] = exchangeEnergy;
			m_exchangeEnergy[j][i] = exchangeEnergy;
			m_coulombEnergy[i][j] = coulombEnergy;
			m_coulombEnergy[j][i] = coulombEnergy;
		}
	}
}

void FockSolver::ElectronRepulsionSolver()
{
	// initailsie the 4D array
	m_eRepulsion = Matrix<Matrix<double>>(m_basisSet.Length(),m_basisSet.Length());
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			m_eRepulsion[i][j] = Matrix<double>(m_basisSet.Length(),m_basisSet.Length());
		}
	}
	#pragma omp parallel for	
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j <= i; j++)
		{
			for (int k = 0; k <= i; k++)
			{
				int lMax;
				if (i == k)
				{
					lMax = j;
				} else
				{
					lMax = k;
				}
				for (int l = 0; l <= lMax; l++)
				{
					m_eRepulsion[i][j][k][l] = m_basisSet[i].ElectronRepulsion(m_basisSet[j],
										  			  m_basisSet[k], m_basisSet[l], m_boyFn);
					m_eRepulsion[i][j][l][k] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[j][i][k][l] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[j][i][l][k] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[k][l][i][j] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[k][l][j][i] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[l][k][i][j] = m_eRepulsion[i][j][k][l];
					m_eRepulsion[l][k][j][i] = m_eRepulsion[i][j][k][l];			
				}
			}
		}
	}
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j <= i; j++)
		{
			for (int k = 0; k <= j; k++)
			{
				for (int l = 0; l <= k; l++)
				{
					std::cout << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << " " << m_eRepulsion[i][j][k][l] << std::endl;
				}
			}
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

void FockSolver::InitialGuess()
{
	m_density = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());

	// This should be part of the matrix class
	for (int i = 0; i < m_basisSet.Length(); ++i)
	{
		for (int j = 0; j < m_basisSet.Length(); ++j)
		{
			m_density[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 1; ++i)
	{
		for (int j = 0; j < m_basisSet.Length(); ++j)
		{
			m_density[i+(j*m_basisSet.Length() / m_basisSet.Length())][j] = 0.0;
		}
	}
}
