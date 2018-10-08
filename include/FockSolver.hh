#ifndef FOCKSOLVER_HH
#define FOCKSOLVER_HH

#include "Matrix.hh"
#include "Vector.hh"
#include "STOnGOrbit.hh"
#include "BoysFunction.hh"

// For now this will follow a restircted method. 

class FockSolver
{
public:

	FockSolver(const Vector<double> &nuclearCharges,
			   const Vector<Vector<double>> &nuclearPositions,
			   const Vector<STOnGOrbit> &basisSet, int nElectrons,
			   const BoysFunction &boyFn);

	~FockSolver();

	double GetGroundEnergy() const;

	void Solve();

private:

	void OneElectronSolver();

	void TwoElectronSolver();

	void ElectronRepulsionSolver();

	double IonPotential();

	void InitialGuess();

private:
	int m_nElectrons;
	Vector<STOnGOrbit> m_basisSet;
	Vector<Vector<double>> m_nuclearPositions;
	Vector<double> m_nuclearCharges;
	BoysFunction m_boyFn;
	Matrix<Matrix<double>> m_eRepulsion;
	double m_groundEnergy;


	Matrix<double> m_density;
	Matrix<double> m_basisOverlap;
	Matrix<double> m_oneElectronEnergy;
	Matrix<double> m_exchangeEnergy;
	Matrix<double> m_coulombEnergy;

	Vector<double> m_energyLevels;

};
#endif