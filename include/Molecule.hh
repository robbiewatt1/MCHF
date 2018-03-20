#ifndef MOLECULE_HH
#define MOLECULE_HH

#include "Vector.hh"
#include "Matrix.hh"
#include "STOnGOrbit.hh"

class Molecule
{
public:

	Molecule();

	Molecule(const Vector<double> &nuclearCharges, const Vector<Vector<double>> &nuclearPositions,
	         int maxL);

	~Molecule();

	Vector<double> GetEnergyLevels();

	Matrix<double> GetBasisCoefficients();

	void CalculateEnergy();

private:

	void SetBasisSet();
	// Sets up the basis set vecotr. All orbitals l = n + m + k <= lMax * number of ions

	double CalculatePotential(int z1, int z2, Vector<double> ionLocation1,
	                    	Vector<double> ionLocation2);
	// Calculate the potential due to two ions

private:

	Vector<double> m_nuclearCharges;			// Vector containing the charge of each ion
	Vector<Vector<double>> m_nuclearPositions;	// Vector containing the position of each ion
	int m_maxL;									// The maximum angular momentum of each basis set

	Vector<STOnGOrbit> m_basisSet;				// Vector of the basis set functions
	Vector<double> m_energyLevels;				// Vector of the energy levels
	Matrix<double> m_basisSetCoefficients;		// maxtix with collumns containing coefficents of basis set
};


#endif