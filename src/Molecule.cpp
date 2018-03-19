#include "Molecule.hh"

Molecule::Molecule()
{
}

Molecule::Molecule(const Vector<double> &nuclearCharges, 
				   const Vector<Vector<double>> &nuclearPositions, int maxL):
	m_nuclearCharges(nuclearCharges), m_nuclearPositions(nuclearPositions), m_maxL(maxL)
{
}

Molecule::~Molecule()
{

}

void Molecule::CalculateEnergy()
{
	// 


}

void SetBasisSet(int max)
{
//	int baseSetNumber = m_maxL * (m_maxL + 1) * (m_maxL + 2) / 6.0
	// Loop over att orbits less than lMax
	for (int k = 0; k <= max; k++)
	{
		for (int m = 0; m <= max; m++)
		{
			for (int n = 0; n < max; n++)
			{
				std::cout << k << " " << m " " << n << std::endl;
			}
		}
	}
}
			
		

/*
		}
		//Loop over all ion locations
		for (int j = 0; j < m_nuclearCharges.Length(); j++)
		{
			

			STOnGOrbit orbital = STOnGOrbit("./OrbitalData/STO3Test", k, m, n,
											m_nuclearPositions[j][0],
											m_nuclearPositions[j][1],
											m_nuclearPositions[j][2])
			m_basisSet.Append(orbital);
		}
		
	}
}
*/