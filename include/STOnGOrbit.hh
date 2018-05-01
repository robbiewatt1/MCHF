#ifndef STONG_ORBIT_HH
#define STONG_ORBIT_HH

#include <string>

#include "Orbital.hh"
#include "GaussianOrbital.hh"
#include "Vector.hh"

class STOnGOrbit : public Orbital
{
public:

	STOnGOrbit();

	STOnGOrbit(std::string dataFile, int k, int m, int n, Vector<double> orbitPosition);

	~STOnGOrbit();

	int GetK();

	int GetM();

	int GetN();

	double GetNormaliseConstant();

	GaussianOrbital GetBaseOribtal(int i) const;

	double GetExponent(int i) const;

	double GetCoefficient(int i) const;

public:

	double Overlap(const STOnGOrbit &orbit);
	// Calculate the overlap intergral of the two orbirtsals

	double KineticOverlap(const STOnGOrbit &orbit);
	// Calculates the kinetic overlap of the two orbitals. This includes the -0.5

	double NuclearOverlap(const STOnGOrbit &orbit, int nuclearCharge, Vector<double> nuclearPosition);
	// Cacultes the nuclear overlap intergral of the orbit with another orbital. the nuclearXYZ
	// are the x y z position of the nuclus that the integral is being calculated around

	Array3D<double> CalculateDataCartesian(const Vector<double> &xAxis,
	                            const Vector<double> &yAxis,
	                            const Vector<double> &zAxis) override;
	// Calcultes the data for a cartesian grid

	Array3D<double> CalculateDataSpherical(const Vector<double> &rAxis,
	                            const Vector<double> &thetaAxis,
	                            const Vector<double> &phiAxis) override;

	// Calcultes the data for a spherical grid

private:

	void Normalise();
	// Orbit should be normilised by defult, but set m_normaliseConstant if not

	void OpenDataFile(std::string fileName);
	// Opens the data file and fills the

	void AllocateBaseOrbitals();
	// Sets the orbital vector

	unsigned int m_k;			// x quntum number of angular momentum
	unsigned int m_m;			// y quntum number of angular momentum
	unsigned int m_n;			// z quntum number of angular momentum
	Vector<double> m_orbitPosition;	// x position of centre of gaussian

	double m_normaliseConstant;	// Factor to normilise the wavefunction
	int m_gaussianNumber;

	Vector<GaussianOrbital> m_baseOrbitalVector;
	// Vector holding the instances of the base orbitals
	Vector<double> m_exponents;
	// Exponents of the gussian basis functions form basis set exchange
	Vector<double> m_coefficeints;
	// Coefficeintsof the gussian basis functions form basis set exchange
};
#endif