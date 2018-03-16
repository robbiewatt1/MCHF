#ifndef STO6G_ORBIT_HH
#define STO6G_ORBIT_HH

#include "Orbital.hh"
#include "Vector.hh"

class STO6GOrbit : public Orbital
{
public:

	STO6GOrbit();

	STO6GOrbit(int k, int m, int n, double centrePositionX, double centrePositionY, 
			   double centrePositionZ);

	~STO6GOrbit();
	
	int GetK();

	int GetM();

	int GetN();

	double GetAlpha();

	double GetPositionX();

	double GetPositionY();

	double GetPositionZ();

	double GetNormaliseConstant();

public:

	double Overlap(const STO6GOrbit &orbit);
	// Calculate the overlap intergral of the two orbirtsals

	double KineticOverlap(const STO6GOrbit &orbit);
	// Calculates the kinetic overlap of the two orbitals. This includes the -0.5

	double NuclearOverlap(const STO6GOrbit &orbit, int nuclearCharge, double nuclearX,
	                      double nuclearY, double nuclearZ);
	// Cacultes the nuclear overlap intergral of the orbit with another orbital. the nuclearXYZ
	// are the x y z position of the nuclus that the integral is being calculated around

private:

	unsigned int m_k;			// x quntum number of angular momentum
	unsigned int m_m;			// y quntum number of angular momentum
	unsigned int m_n;			// z quntum number of angular momentum
	double m_centrePositionX;	// x position of centre of gaussian
	double m_centrePositionY;	// x position of centre of gaussian
	double m_centrePositionZ;	// x position of centre of gaussian
	double m_normaliseConstant;	// Factor to normilise the wavefunction

	// Exponents of the gussian basis functions form basis set exchange
	Vector<double> m_exponents(6); 
	m_exponents[0] = 0.100112428;
	m_exponents[1] = 0.243076747;
	m_exponents[2] = 0.625955266;
	m_exponents[3] = 1.822142904;
	m_exponents[4] = 6.513143725;
	m_exponents[5] = 6.513143725;

	// Coefficeintsof the gussian basis functions form basis set exchange
	Vector<double> m_coefficeints(6); 
	m_coefficeints[0] = 0.13033408410;
	m_coefficeints[1] = 0.41649152980;
	m_coefficeints[2] = 0.37056279970;
	m_coefficeints[3] = 0.16853830490;
	m_coefficeints[4] = 0.04936149294;
	m_coefficeints[5] = 0.00916359628;
};
#endif