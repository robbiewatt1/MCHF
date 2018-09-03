#ifndef GAUSSIANORBITAL_HH
#define GAUSSIANORBITAL_HH

#include "Orbital.hh"
#include "BoysFunction.hh"

class GaussianOrbital : public Orbital
{

public:

	GaussianOrbital();	// Default constructor

	GaussianOrbital(int k, int m, int n, double alpha, Vector<double> orbitPosition);

	~GaussianOrbital();

	int GetK() const;

	int GetM() const;

	int GetN() const;

	double GetAlpha() const;

	double GetPositionX() const;

	double GetPositionY() const;

	double GetPositionZ() const;

	double GetNormaliseConstant() const;

public:

	double Overlap(const GaussianOrbital &orbit) const;
	// Calculate the overlap intergral of the two orbirtsals

	double KineticOverlap(const GaussianOrbital &orbit) const;
	// Calculates the kinetic overlap of the two orbitals. This includes the -0.5

	double NuclearOverlap(const GaussianOrbital &orbit, int nuclearCharge,
	                      const Vector<double> &nuclearPosition, const BoysFunction &boyFn) const;
	// Cacultes the nuclear overlap intergral of the orbit with another orbital. the nuclearXYZ
	// are the x y z position of the nuclus that the integral is being calculated around

	double ElectronRepulsion(const GaussianOrbital &orbit1, const GaussianOrbital &orbit2, 
							 const GaussianOrbital &orbit3, const BoysFunction &boyFn) const;
	// Calculates the repulsion between gaussian orbitals. Only the hartree and fock terms 
	// are considered

	Vector<double> MatrixElement(const GaussianOrbital &orbit) const;

	Array3D<double> CalculateDataCartesian(const Vector<double> &xAxis,
	                            const Vector<double> &yAxis,
	                            const Vector<double> &zAxis) override;
	// Calcultes the data for a cartesian grid

	Array3D<double> CalculateDataSpherical(const Vector<double> &rAxis,
	                            const Vector<double> &thetaAxis,
	                            const Vector<double> &phiAxis) override;

	// Calcultes the data for a spherical grid

public:

	void Normalise();
	// Normalises the wavefunction, and sets the norm constant. uses the overlap of the orbit with itself

	double OverlapFunction(int l1, int l2, double gamma, double posA, double posB) const;

	double NuclearFunction(int l, int r, int i, int l1, int l2, double alpha, double beta,
	                       double pos1, double pos2, double pos3) const;

	double ElectronFunction(int l, int l_p, int r, int r_p, int i, int l1, int l2, double posA,
							double posB, double posP, double gamma1, int l3, int l4, double posC, 
							double posD, double posQ, double gamma2) const;

	double ThetaFn(int l, int l1, int l2, double a, double b, double c, double d) const;

	double GaussianProduct(int k, int l1, int l2, double pos1, double pos2) const;

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