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
	double alpha1 = 35.52322122;// 0.00916359628;
	double alpha2 = 6.513143725;// 0.04936149294;
	double alpha3 = 1.822142904;// 0.16853830490;
	double alpha4 = 0.625955266;// 0.37056279970;
	double alpha5 = 0.243076747;// 0.41649152980;
	double alpha6 = 0.100112428;// 3.13033408410;
	Vector<GaussianOrbital> orbitV(24);
	orbitV[0] = GaussianOrbital(0,0,0,alpha1,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(0,0,0,alpha2,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(0,0,0,alpha3,0.0,0.0,0.0);
	orbitV[3] = GaussianOrbital(0,0,0,alpha4,0.0,0.0,0.0);
	orbitV[4] = GaussianOrbital(0,0,0,alpha5,0.0,0.0,0.0);
	orbitV[5] = GaussianOrbital(0,0,0,alpha6,0.0,0.0,0.0);

	orbitV[6] = GaussianOrbital(2,0,0,alpha1,0.0,0.0,0.0);
	orbitV[7] = GaussianOrbital(2,0,0,alpha2,0.0,0.0,0.0);
	orbitV[8] = GaussianOrbital(2,0,0,alpha3,0.0,0.0,0.0);
	orbitV[9] = GaussianOrbital(2,0,0,alpha4,0.0,0.0,0.0);
	orbitV[10] = GaussianOrbital(2,0,0,alpha5,0.0,0.0,0.0);
	orbitV[11] = GaussianOrbital(2,0,0,alpha6,0.0,0.0,0.0);

	orbitV[12] = GaussianOrbital(0,2,0,alpha1,0.0,0.0,0.0);
	orbitV[13] = GaussianOrbital(0,2,0,alpha2,0.0,0.0,0.0);
	orbitV[14] = GaussianOrbital(0,2,0,alpha3,0.0,0.0,0.0);
	orbitV[15] = GaussianOrbital(0,2,0,alpha4,0.0,0.0,0.0);
	orbitV[16] = GaussianOrbital(0,2,0,alpha5,0.0,0.0,0.0);
	orbitV[17] = GaussianOrbital(0,2,0,alpha6,0.0,0.0,0.0);

	orbitV[18] = GaussianOrbital(0,0,2,alpha1,0.0,0.0,0.0);
	orbitV[19] = GaussianOrbital(0,0,2,alpha2,0.0,0.0,0.0);
	orbitV[20] = GaussianOrbital(0,0,2,alpha3,0.0,0.0,0.0);
	orbitV[21] = GaussianOrbital(0,0,2,alpha4,0.0,0.0,0.0);
	orbitV[22] = GaussianOrbital(0,0,2,alpha5,0.0,0.0,0.0);
	orbitV[23] = GaussianOrbital(0,0,2,alpha6,0.0,0.0,0.0);


	Matrix<double> SMatrix(24,24);
	Matrix<double> EMatrix(24,24);
	Matrix<double> eigenVectors(24,24);
	Vector<double> eigenvalues(24);

	for (int i = 0; i < 24; i++)
	{
		for (int j = 0; j < 24; j++)
		{
			EMatrix[i][j] = ( orbitV[i].KineticOverlap(orbitV[j])) + orbitV[i].NuclearOverlap(orbitV[j],1,0,0,0);
			SMatrix[i][j] = orbitV[i].Overlap(orbitV[j]) * orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
		}
	}

	LinearAlgebra::GeneralisedEigenSolver(EMatrix,SMatrix,eigenVectors,eigenvalues);
	eigenVectors.Print();
	eigenvalues.Print();
	return 0;
}
