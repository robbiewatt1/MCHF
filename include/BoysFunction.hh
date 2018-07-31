#ifndef BOYSFUNCTION_HH
#define BOYSFUNCTION_HH

#include "Matrix.hh"
#include "Vector.hh"

class BoysFunction
{
public:
	// This is a class that is used to calculate the boys function which is a standard integral. 
	// The maximum values for u and v should be found before and then passed to the class. 

	BoysFunction(int maxV, double maxU, int uLength);

	~BoysFunction();

	double Interpolate(int v, double u) const;

private:
	void BoysGenerator();

private:
	int m_uLength;
	int m_vLength;

	Matrix<double> m_boysData;
	Vector<double> m_uAxis;
	Vector<int>	   m_vAxis;	
};
#endif