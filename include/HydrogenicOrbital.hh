#ifndef HYDROGENICORBITAL__HH
#define HYDROGENICORBITAL__HH

#include "Orbital.hh"
#include "Matrix.hh"

class HydrogenicOrbital : Orbital
{
public:
	HydrogenicOrbital();

	HydrogenicOrbital(int maxL);

	~HydrogenicOrbital();

	void SetCoefficent(int l, int ml, double value);

	int GetMaxL();

	void MinimiseEnergy();

	void CalculateDataCartesian(const Vector<double> &xAxis,
	                            const Vector<double> &yAxis,
	                            const Vector<double> &zAxis) override;

	void CalculateDataSpherical(const Vector<double> &rAxis,
	                            const Vector<double> &thetaAxis,
	                            const Vector<double> &phiAxis) override;


public:

	Matrix<double> CalculateHMatrix();

	Matrix<double> CalculateSMatrix();

	int NumberOfOrbitals();

private:
	int m_maxL;

	double m_energy;

	double **m_coefficents;
};
#endif