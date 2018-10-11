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

double FockSolver::GetGroundEnergy() const
{
	return m_groundEnergy;
}

void FockSolver::GroundSolve(int itrMax, int nInterp)
{
	InitialGuess();
	OneElectronSolver();
	ElectronRepulsionSolver();

	for (int i = 0; i < itrMax; i++)
	{
		TwoElectronSolver();
		Matrix<double> fockMaxtrix = m_oneElectronEnergy + m_coulombEnergy - m_exchangeEnergy;
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

		LinearAlgebra::EigenSolver(fockMaxtrix, m_orbitalCoeff, m_fockyLevels);
		m_orbitalCoeff = othgTransMatrix * m_orbitalCoeff;
		
		// only occupided oribtals for the aubjua lowest energy
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			for (int k = 0; k <  m_basisSet.Length(); k++)
			{
				double sum(0);
				for (int l = 0; l < m_nElectrons; l++)
				{
					sum += m_orbitalCoeff[k][l] * m_orbitalCoeff[j][l];
				}
				m_density[j][k] = sum;
			}
		}
	}

	double sum(0);
	for (int j = 0; j < m_basisSet.Length(); ++j)
	{
		for (int k = 0; k < m_basisSet.Length(); ++k)
		{
			sum += m_density[j][k] * (m_oneElectronEnergy[j][k] + 0.5 * (m_coulombEnergy[j][k] - m_exchangeEnergy[j][k]));
		}
	}
	m_groundEnergy = sum + IonPotential();
}

void FockSolver::ExcitedSolver(int itrMax)
{
	// check that groun state has been solved.
	if (m_fockyLevels.Length() == 0)
	{
		std::cerr << "Error: Must find ground state first." << std::endl;
		exit(-1);
	}
	// Form the starting density matrix using the next highest virtual
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j <  m_basisSet.Length(); j++)
		{
			double sum(0);
			for (int k = 0; k < m_nElectrons; k++)
			{
				if (k == m_nElectrons+1)
				{
					continue;
				} else
				{
					sum += m_orbitalCoeff[j][k] * m_orbitalCoeff[i][k];
				}
			}
			m_density[i][j] = sum;
		}
	}

	// Enter the normal fock loop
	Matrix<double> orbitalCoeffNew;
	for (int i = 0; i < itrMax; i++)
	{
		TwoElectronSolver();
		Matrix<double> fockMaxtrix = m_oneElectronEnergy + m_coulombEnergy - m_exchangeEnergy;

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
		LinearAlgebra::EigenSolver(fockMaxtrix, orbitalCoeffNew, m_fockyLevels);
		orbitalCoeffNew = othgTransMatrix * orbitalCoeffNew;
		
		Matrix<double> maxOverlap = LinearAlgebra::Transpose(m_orbitalCoeff) * m_basisOverlap * orbitalCoeffNew;
		Vector<double> projections(maxOverlap.GetRows());
		for (int j = 0; j < maxOverlap.GetRows(); ++j)
		{
			for (int k = 0; k < maxOverlap.GetColumns(); k++)
			{
				projections[j] += (maxOverlap[j][k]);
			}
		}
		projections.Print();
		// find the n largest projections
		Vector<int> indexVector(m_nElectrons);
		for (int j = 0; j < m_nElectrons; j++)
		{
			int maxIndex = LinearAlgebra::MaxValueLocation(projections);
			indexVector[j] = j;
			projections[maxIndex] = 0;
		}
		// only occupided oribtals for the aubjua lowest energy
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			for (int k = 0; k <  m_basisSet.Length(); k++)
			{
				double sum(0);
				for (int l = 0; l < m_nElectrons; l++)
				{
					sum += orbitalCoeffNew[k][indexVector[l]] * orbitalCoeffNew[j][indexVector[l]];
				}
				m_density[j][k] = sum;
			}
		}
	}
	double sum(0);
	for (int j = 0; j < m_basisSet.Length(); ++j)
	{
		for (int k = 0; k < m_basisSet.Length(); ++k)
		{
			sum += m_density[j][k] * (m_oneElectronEnergy[j][k] + 0.5 * (m_coulombEnergy[j][k] - m_exchangeEnergy[j][k]));
		}
	}
	double excitedEnergy;
	m_groundEnergy = sum + IonPotential();
	std::cout << m_groundEnergy << std::endl;
}



void FockSolver::OneElectronSolver()
{
	m_basisOverlap = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	m_oneElectronEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());

	// calculate the ion repulsion potenital
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
			m_oneElectronEnergy[i][j] = kineticEnergy + nuclearPotential;
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
}

double FockSolver::IonPotential()
{
	double potential(0);
	for (int i = 0; i < m_nuclearCharges.Length(); i++)
	{
		for (int j = i + 1; j < m_nuclearCharges.Length(); j++)
		{
			potential += (m_nuclearCharges[i] * m_nuclearCharges[j]) / 
						 (std::sqrt(std::pow(m_nuclearPositions[i][0] - m_nuclearPositions[j][0], 2.0)
	                			  + std::pow(m_nuclearPositions[i][1] - m_nuclearPositions[j][1], 2.0)
	                			  + std::pow(m_nuclearPositions[i][2] - m_nuclearPositions[j][2], 2.0)));
		}
	}
	return potential; 
}

void FockSolver::InitialGuess()
{
	m_density = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());

	// This should be part of the matrix class
	for (int i = 0; i < m_basisSet.Length(); ++i)
	{
		for (int j = 0; j < m_basisSet.Length(); ++j)
		{
			if (i == j && i == 0)
			{
				m_density[i][j] = 1.0;
			} else
			{
				m_density[i][j] = 1.0;
			}
		}
	}
}
