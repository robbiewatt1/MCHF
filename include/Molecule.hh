#ifndef MOLECULE_HH
#define MOLECULE_HH

#include <string>
#include "Vector.hh"
#include "Matrix.hh"
#include "Array3D.hh"
#include "STOnGOrbit.hh"
#include "BoysFunction.hh"

class Molecule
{
public:

	Molecule();

	Molecule(int nElectrons, const Vector<double> &nuclearCharges, 
			 const Vector<Vector<double>> &nuclearPositions, int maxL,
			 const BoysFunction &boyFn, std::string basisSetDir);

	~Molecule();

	Vector<double> GetEnergyLevels();

	Matrix<double> GetBasisCoefficients();

	double CalculateEnergy();

	Vector<double> MatrixElement(int level1, int level2);

	double OscilatorStrength(int level1, int level2);

	void SetXAxis(const Vector<double> &xAxis);

	void SetYAxis(const Vector<double> &yAxis);

	void SetZAxis(const Vector<double> &zAxis);

	void SetWavefunction(const Vector<double> &xAxis, const Vector<double> &yAxis,
						 const Vector<double> &zAxis);

	Array3D<double> CalculateWavefunction(int level);

	void OutputData(int level, std::string fileName);

private:

	void SetBasisSet(std::string basisSetDir);
	// Sets up the basis set vecotr. All orbitals l = n + m + k <= lMax * number of ions

	double CalculatePotential(int z1, int z2, Vector<double> ionLocation1,
	                    	Vector<double> ionLocation2);
	// Calculate the potential due to two ions

private:

	Vector<double> m_nuclearCharges;			// Vector containing the charge of each ion
	Vector<Vector<double>> m_nuclearPositions;	// Vector containing the position of each ion
	int m_maxL;									// The maximum angular momentum of each basis set
	int m_nElectrons;							// Number of electrons in molecule

	Vector<STOnGOrbit> m_basisSet;				// Vector of the basis set functions
	Vector<double> m_energyLevels;				// Vector of the energy levels
	Matrix<double> m_basisSetCoefficients;		// maxtix with collumns containing coefficents of basis set

	Vector<double> m_xAxis;						// Xasis of the wavefunction
	Vector<double> m_yAxis;						// Xasis of the wavefunction
	Vector<double> m_zAxis;						// Xasis of the wavefunction

	BoysFunction m_boyFn;						// The BoysFunction data  
};


#endif