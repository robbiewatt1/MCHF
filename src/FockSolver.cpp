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
	m_coeffC.Print();
	getchar();

	// Solve for the one electron matrix. This only needs to be done once
	OneElectronSolver();

	for (int i = 0; i < 50; i++)
	{
		
		TwoElectronSolver();

		Matrix<double> fockMaxtrix = m_oneElectronEnergy + m_coulombEnergy - m_exchangeEnergy;
		// Find matgrix used to transform to othogonal basis. Here I might need to remove bad matricies!
		(LinearAlgebra::Transpose(m_coeffC) * m_exchangeEnergy * m_coeffC).Print();
		(LinearAlgebra::Transpose(m_coeffC) * m_coulombEnergy * m_coeffC).Print();
		Vector<double> othgTransVector;
		Matrix<double> othgTransMatrix;
		LinearAlgebra::EigenSolver(m_basisOverlap, othgTransMatrix, othgTransVector);
		fockMaxtrix.Print();
		for (int j = 0; j < m_basisOverlap.GetRows(); j++)
		{
			for (int k = 0; k < m_basisOverlap.GetColumns(); k++)
			{
				othgTransMatrix[j][k] /=  std::sqrt(othgTransVector[k]);
			}
		}

		// Transform to the othogonal coods
		fockMaxtrix = LinearAlgebra::Transpose(othgTransMatrix) * fockMaxtrix * othgTransMatrix;
		

		// Solve the eigenvalue HF equation
		Matrix<double> orbitalCoeff;
		LinearAlgebra::EigenSolver(fockMaxtrix, orbitalCoeff, m_energyLevels);

		// Reduce the matrix back to correct size
		orbitalCoeff = othgTransMatrix * orbitalCoeff;
	//	orbitalCoeff.Print();

		for (int j = 0; j < m_nElectrons; j++)
		{
			for (int k = 0; k < m_basisSet.Length(); k++)
			{
				m_coeffC[k][j] = orbitalCoeff[k][j];			
			}
		}
		m_coeffC.Print();
		m_energyLevels.Print();
		std::cout << "\n\n\n\n";
	}

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

void FockSolver::TwoElectronSolver()
{
	m_coulombEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	m_exchangeEnergy = Matrix<double>(m_basisSet.Length(), m_basisSet.Length());
	Matrix<double> weights(m_basisSet.Length(), m_basisSet.Length());
	// Check this is in the right order
	weights = m_coeffC * LinearAlgebra::Transpose(m_coeffC);
	#pragma omp parallel for	
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = i; j < m_basisSet.Length(); j++)
		{
			double exchangeEnergy(0);
			double coulombEnergy(0);
			for (int k = 0; k < m_basisSet.Length(); k++)
			{
				for (int l = k; l < m_basisSet.Length(); l++)
				{
					if (weights[k][l] < 1e-10)
					{
						continue;
					}else
					{
						exchangeEnergy += weights[k][l] * m_basisSet[i].ElectronRepulsion(m_basisSet[l],
										  m_basisSet[k], m_basisSet[j], m_boyFn);
						coulombEnergy += weights[k][l] * m_basisSet[i].ElectronRepulsion(m_basisSet[j],
										  m_basisSet[k], m_basisSet[l], m_boyFn);
					}
				}
			}
			m_exchangeEnergy[i][j] = exchangeEnergy;
			m_exchangeEnergy[j][i] = exchangeEnergy;
			m_coulombEnergy[i][j] = coulombEnergy;
			m_coulombEnergy[j][i] = coulombEnergy;
		}
	}
}

/*
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
				for (int l = 0; l < m_basisSet.Length(); l++)
				{
					if (weights[k][l] == 0)
					{
						continue;
					} else
					{
						exchangeEnergy += weights[k][l] * m_basisSet[i].ElectronRepulsion(m_basisSet[l],
					    	              m_basisSet[k], m_basisSet[j], m_boyFn);
					}
				}
			}
			m_exchangeEnergy[i][j] = exchangeEnergy;
	//		m_exchangeEnergy[j][i] = exchangeEnergy;
			// Think this should be symetirc but will need to check
		}
	}
}
*/
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
	m_coeffC = Matrix<double>(m_basisSet.Length(), m_nElectrons);

	// This should be part of the matrix class
	for (int i = 0; i < m_basisSet.Length(); ++i)
	{
		for (int j = 0; j < m_nElectrons; ++j)
		{
			m_coeffC[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 1; ++i)
	{
		for (int j = 0; j < m_nElectrons; ++j)
		{
			m_coeffC[i+(j*m_basisSet.Length() / m_nElectrons)][j] = 1.0;
		}
	}
}
