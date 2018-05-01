#ifndef GAUSSIANORBITAL_HH
#define GAUSSIANORBITAL_HH

#include "Orbital.hh"

class GaussianOrbital : public Orbital
{

public:

	GaussianOrbital();	// Default constructor

	GaussianOrbital(int k, int m, int n, double alpha, Vector<double> orbitPosition);

	~GaussianOrbital();

	int GetK();

	int GetM();

	int GetN();

	double GetAlpha();

	double GetPositionX();

	double GetPositionY();

	double GetPositionZ();

	double GetNormaliseConstant();

public:

	double Overlap(const GaussianOrbital &orbit);
	// Calculate the overlap intergral of the two orbirtsals

	double KineticOverlap(const GaussianOrbital &orbit);
	// Calculates the kinetic overlap of the two orbitals. This includes the -0.5

	double NuclearOverlap(const GaussianOrbital &orbit, int nuclearCharge,
	                      Vector<double> nuclearPosition);
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
	// Normalises the wavefunction, and sets the norm constant. uses the overlap of the orbit with itself

	double OverlapFunction(int l1, int l2, double gamma, double posA, double posB);

	double NuclearFunction(int l, int r, int i, int l1, int l2, double alpha, double beta,
	                       double pos1, double pos2, double pos3);

	double GaussianProduct(int k, int l1, int l2, double pos1, double pos2);

	double GaussianProduct2(int k, int l1, int l2, double pos1, double pos2);

	double GaussianProduct3(int k, int l1, int l2, double pos1, double pos2);

	GaussianOrbital ChangeK(int k) const;

	GaussianOrbital ChangeM(int m) const;

	GaussianOrbital ChangeN(int n) const;

private:
	unsigned int m_k;			// x quntum number of angular momentum
	unsigned int m_m;			// y quntum number of angular momentum
	unsigned int m_n;			// z quntum number of angular momentum
	double m_alpha;				// a factor in the exponetial (exp(ax^2))
	Vector<double> m_orbitPosition;	// position of centre of gaussian
	double m_normaliseConstant;	// Factor to normilise the wavefunction
};
#endif