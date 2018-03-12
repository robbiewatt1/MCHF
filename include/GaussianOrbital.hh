#ifndef GAUSSIANORBITAL_HH
#define GAUSSIANORBITAL_HH

#include "Orbital.hh"

class GaussianOrbital : public Orbital
{

public:

	GaussianOrbital();		// Default constructor

	GaussianOrbital(int k, int m, int n, double centrePositionX, double centrePositionY,
	                double centrePositionZ);

	~GaussianOrbital();

	int GetK();

	int GetM();

	int GetN();

	double GetPositionX();

	double GetPositionY();

	double GetPositionZ();

	double GetNormaliseConstant();

public:

	double Overlap(const GaussianOrbital &orbit);
	// Calculate the overlap intergral of the two orbirtsals

	double KineticOverlap(const GaussianOrbital &orbit);
	// Calculates the kinetic overlap of the two orbitals. This includes the -0.5

	double NuclearOverlap(const GaussianOrbital &orbit, int nuclearCharge, double nuclearX,
	                      double nuclearY, double nuclearZ);
	// Cacultes the nuclear overlap intergral of the orbit with another orbital. the nuclearXYZ
	// are the x y z position of the nuclus that the integral is being calculated around

	void CalculateDataCartesian(const Vector<double> &xAxis,
	                            const Vector<double> &yAxis,
	                            const Vector<double> &zAxis) override;
	// Calcultes the data for a cartesian grid

	void CalculateDataSpherical(const Vector<double> &rAxis,
	                            const Vector<double> &thetaAxis,
	                            const Vector<double> &phiAxis) override;

	// Calcultes the data for a spherical grid

private:

	void Normalise();

	double OverlapFunction(int l1, int l2, double posA, double posB);

	double NuclearFunction(int l, int r, int i, int l1, int l2, double pos1, double pos2, double pos3);

	double GaussianProduct(int k, int l1, int l2, double pos1, double pos2);

	double GaussianProduct2(int k, int l1, int l2, double pos1, double pos2);

	GaussianOrbital ChangeK(int k) const;

	GaussianOrbital ChangeM(int m) const;

	GaussianOrbital ChangeN(int n) const;

private:
	unsigned int m_k;	// x quntum number of angular momentum
	unsigned int m_m;	// y quntum number of angular momentum
	unsigned int m_n;	// z quntum number of angular momentum

	double m_centrePositionX;
	double m_centrePositionY;
	double m_centrePositionZ;
	double m_normaliseConstant;

};
#endif