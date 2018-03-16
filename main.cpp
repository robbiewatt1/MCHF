#include "Functions.hh"
#include "Orbital.hh"
#include "GaussianOrbital.hh"
#include "Constants.hh"
#include "Array3D.hh"
#include "Numerics.hh"
#include "HydrogenicOrbital.hh"
#include "LinearAlgebra.hh"
#include "Matrix.hh"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <climits>
#include <string>
#include "H5Cpp.h"
#include <cmath>
int main()
{
	double alpha = 0.003;
	double alpha2 = 0.034;
	double alpha3 = 0.8135;
	Vector<GaussianOrbital> orbitV(16);
	orbitV[0] = GaussianOrbital(1,1,1,0.3,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(2,1,1,alpha2,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(1,2,1,alpha3,0.0,0.0,0.0);
	orbitV[3] = GaussianOrbital(1,1,2,alpha2,0.0,0.0,0.0);
	orbitV[4] = GaussianOrbital(2,2,1,alpha2,0.0,0.0,0.0);
	orbitV[5] = GaussianOrbital(2,1,2,alpha2,0.0,0.0,0.0);
	orbitV[6] = GaussianOrbital(1,2,2,alpha2,0.0,0.0,0.0);
	orbitV[7] = GaussianOrbital(3,1,1,alpha2,0.0,0.0,0.0);
	orbitV[8] = GaussianOrbital(1,3,1,alpha2,0.0,0.0,0.0);
	orbitV[9] = GaussianOrbital(1,1,3,alpha2,0.0,0.0,0.0);
    orbitV[10] = GaussianOrbital(5,1,1,alpha2,0.0,0.0,0.0);
    orbitV[11] = GaussianOrbital(1,5,1,alpha2,0.0,0.0,0.0);
    orbitV[12] = GaussianOrbital(1,1,5,alpha2,0.0,0.0,0.0);
    orbitV[13] = GaussianOrbital(7,1,1,alpha2,0.0,0.0,0.0);
    orbitV[14] = GaussianOrbital(1,7,1,alpha2,0.0,0.0,0.0);
    orbitV[15] = GaussianOrbital(1,1,7,alpha2,0.0,0.0,0.0);

	Matrix<double> SMatrix(16,16);
	Matrix<double> EMatrix(16,16);
	Matrix<double> eigenVectors(16,16);
	Vector<double> eigenvalues(16);

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			EMatrix[i][j] = ( orbitV[i].KineticOverlap(orbitV[j])) + orbitV[i].NuclearOverlap(orbitV[j],1,0,0,0);
			SMatrix[i][j] = orbitV[i].Overlap(orbitV[j]) * orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
		}
	}
//	SMatrix.Print();
	LinearAlgebra::GeneralisedEigenSolver(EMatrix,SMatrix,eigenVectors,eigenvalues);
	eigenVectors.Print();
	eigenvalues.Print();
	
	std::cout << orbitV[0].KineticOverlap(orbitV[0]) + orbitV[0].NuclearOverlap(orbitV[0],1,0,0,0) << std::endl;

	return 0;
}
