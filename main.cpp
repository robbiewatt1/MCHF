#include "Functions.hh"
#include "Orbital.hh"
#include "GaussianOrbital.hh"
#include "Constants.hh"
#include "Array3D.hh"
#include "Numerics.hh"
#include "HydrogenicOrbital.hh"
#include "LinearAlgebra.hh"
#include "Matrix.hh"
#include "STOnGOrbit.hh"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <climits>
#include <string>
#include "H5Cpp.h"
#include <cmath>
int main()
{
	
	STOnGOrbit a = STOnGOrbit("./OrbitalData/STO3",0,0,0,0,0,0);
	std::cout << (a.KineticOverlap(a) + a.NuclearOverlap(a,1,0,0,0)) * std::pow(a.GetNormaliseConstant(),2.0) << std::endl;
	std::cout << "kin " << (a.KineticOverlap(a) * std::pow(a.GetNormaliseConstant(),2.0)) << std::endl;
	std::cout << "pot " << (a.NuclearOverlap(a,1,0,0,0) * std::pow(a.GetNormaliseConstant(),2.0)) << std::endl;
	std::cout << std::pow(a.GetNormaliseConstant(),2.0) << std::endl;

	Vector<GaussianOrbital> orbitV(3);
	double alpha1 = 3.42525091;
	double alpha2 = 0.62391373;
	double alpha3 = 0.16885540;
	orbitV[0] = GaussianOrbital(0,0,0,alpha1,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(0,0,0,alpha2,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(0,0,0,alpha3,0.0,0.0,0.0);
	
	Matrix<double> SMatrix(3,3);
	Matrix<double> EMatrix(3,3);
	Matrix<double> eigenVectors(3,3);
	Vector<double> eigenvalues(3);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			EMatrix[i][j] = (( orbitV[i].KineticOverlap(orbitV[j])) + orbitV[i].NuclearOverlap(orbitV[j],1,0,0,0) )* orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
			SMatrix[i][j] = orbitV[i].Overlap(orbitV[j]) * orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
		}
	}
	std::cout << "M1" << std::endl;
	EMatrix.Print();
	std::cout << "M2" << std::endl;
	SMatrix.Print();

	LinearAlgebra::GeneralisedEigenSolver(EMatrix,SMatrix,eigenVectors,eigenvalues);
	std::cout << "RESULT " << std::endl;
	eigenVectors.Print();
	eigenvalues.Print();
return 0;
}
