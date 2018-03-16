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
	double alpha = 0.0334946;
	double alpha2 = 0.2347269;
	double alpha3 = 0.8137573;
	Vector<GaussianOrbital> orbitV(12);
	orbitV[0] = GaussianOrbital(0,0,0,alpha,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(0,0,0,alpha2,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(0,0,0,alpha3,0.0,0.0,0.0);
	orbitV[3] = GaussianOrbital(2,0,0,alpha,0.0,0.0,0.0);
	orbitV[4] = GaussianOrbital(2,0,0,alpha2,0.0,0.0,0.0);
	orbitV[5] = GaussianOrbital(2,0,0,alpha3,0.0,0.0,0.0);
	orbitV[6] = GaussianOrbital(0,2,0,alpha,0.0,0.0,0.0);
	orbitV[7] = GaussianOrbital(0,2,0,alpha2,0.0,0.0,0.0);
	orbitV[8] = GaussianOrbital(0,2,0,alpha3,0.0,0.0,0.0);
	orbitV[9] = GaussianOrbital(0,0,2,alpha,0.0,0.0,0.0);
	orbitV[10] = GaussianOrbital(0,0,2,alpha2,0.0,0.0,0.0);
	orbitV[11] = GaussianOrbital(0,0,2,alpha3,0.0,0.0,0.0);

/*

	orbitV[10] = GaussianOrbital(5,1,1,alpha2,0.0,0.0,0.0);
    orbitV[11] = GaussianOrbital(1,5,1,alpha2,0.0,0.0,0.0);
    orbitV[12] = GaussianOrbital(1,1,5,alpha2,0.0,0.0,0.0);
    orbitV[13] = GaussianOrbital(7,1,1,alpha2,0.0,0.0,0.0);
    orbitV[14] = GaussianOrbital(1,7,1,alpha2,0.0,0.0,0.0);
    orbitV[15] = GaussianOrbital(1,1,7,alpha2,0.0,0.0,0.0);
*/
	Matrix<double> SMatrix(12,12);
	Matrix<double> EMatrix(12,12);
	Matrix<double> eigenVectors(12,12);
	Vector<double> eigenvalues(12);

	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 12; j++)
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
