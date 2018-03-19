#ifndef ORBITAL_HH
#define ORBITAL_HH

#include "Vector.hh"
#include "Array3D.hh"

class Orbital
{
public:

	Orbital();

	Orbital(const Vector<double> &xAxis, const Vector<double> &yAxis, 
			const Vector<double> &zAxis);

	Orbital(const Orbital &orbital);

	~Orbital();

	void SetXAxis(const Vector<double> &xAxis);

	void SetYAxis(const Vector<double> &yAxis);

	void SetZAxis(const Vector<double> &zAxis);	

	void SetData(const Vector<double> &xAxis, const Vector<double> &yAxis, 
				 const Vector<double> &zAxis);

	Vector<double> GetXAxis() const;

	Vector<double> GetYAxis() const;

	Vector<double> GetZAxis() const;

	Array3D<double> GetData() const;

	// Function that outputs the orbital data to an H5 file
	void OutputData();

	// Calculates the values at each element of data assuming a cartesian grid
	virtual void CalculateDataCartesian(const Vector<double> &xAxis,
	                           		   const Vector<double> &yAxis,
	                           		   const Vector<double> &zAxis) = 0;

	// Calculates the values at each element in data assuming a spherical grid
	virtual void CalculateDataSpherical(const Vector<double> &rAxis,
	                           		   const Vector<double> &thetaAxis,
	                                   const Vector<double> &phiAxis) = 0;


protected:
	Array3D<double> m_data;
	Vector<double> m_xAxis;
	Vector<double> m_yAxis;
	Vector<double> m_zAxis;
};
#endif