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

void FockSolver::Solve(int itrMax, int nInterp)
{
	InitialGuess();
	OneElectronSolver();
	ElectronRepulsionSolver();

	Vector<double> fockCoeffs(nInterp);
	Vector<Matrix<double>> fockVector(nInterp);
	Vector<Matrix<double>> errorVector(nInterp);
	for (int i = 0; i < 6; i++)
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

		// apply the DIIS process to improve convergence
		Matrix<double> errorMatrix = LinearAlgebra::Transpose(othgTransMatrix) 
						            * (fockMaxtrix * m_density * m_basisOverlap
						            - m_basisOverlap * m_density * fockMaxtrix) * othgTransMatrix;

		// Fill the error and fock Vectors
		if (i < nInterp)
		{
			errorVector[i] = errorMatrix;
			fockVector[i]  = fockMaxtrix;
		} else
		{
			errorVector[i%nInterp] = errorMatrix;
			fockVector[i%nInterp] = fockMaxtrix;
		}

		// Form the Beta matrix and solution vector
		if (i = 0)
		{
			fockCoeffs[0] = 1;
			fockMaxtrix = LinearAlgebra::Transpose(othgTransMatrix) * fockMaxtrix * othgTransMatrix;
		}else
		{
			Matrix<double> beta(i+2,i+2);
			Vector<double> alpha(i+2);
			for (int j = 0; j < i+2; j++)
			{
				if (j == i + 1)
				{
					alpha[j] = -1;
				}else
				{
					alphaj[j] = 0;
				}
				for (int k = 0; k < i+2; k++)
				{
					if (j == i + 1 || k == i + 1)
					{
						beta[j][k] = -1;
					} else
					{
						beta[j][k] = LinearAlgebra::Trace(LinearAlgebra::Transpose(errorVector[j]) * errorVector[k]);
					}
				}
			}
		}

//		if (i < nErrors)
//		{
//			fockVector[i] = fockMaxtrix;
//		}

		Matrix<double> orbitalCoeff;
		LinearAlgebra::EigenSolver(fockMaxtrix, orbitalCoeff, m_energyLevels);
		orbitalCoeff = othgTransMatrix * orbitalCoeff;
		
		// only occupided oribtals for the aubjua lowest energy
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			for (int k = 0; k <  m_basisSet.Length(); k++)
			{
				double sum(0);
				for (int l = 0; l < m_nElectrons; l++)
				{
					sum += orbitalCoeff[k][l] * orbitalCoeff[j][l];
				}
				m_density[j][k] = sum;
			}
		}
	}


//	m_groundEnergy = LinearAlgebra::Trace(m_oneElectronEnergy * m_density) + 0.5 * LinearAlgebra::Trace((2.0 * m_coulombEnergy - m_exchangeEnergy) * m_density) + IonPotential(); 
//	sum + IonPotential();
	double sum(0);
	for (int j = 0; j < m_basisSet.Length(); ++j)
	{
		for (int k = 0; k < m_basisSet.Length(); ++k)
		{
			sum += m_density[j][k] * (m_oneElectronEnergy[j][k] + 0.5 * (m_coulombEnergy[j][k] - m_exchangeEnergy[j][k]));
		}
	}
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
