#include <cmath>
#include <iostream>
#include <fstream>

#include "Molecule.hh"
#include "Vector.hh"
#include "STOnGOrbit.hh"
#include "Numerics.hh"
#include "GaussianOrbital.hh"
#include "Functions.hh"
#include "BoysFunction.hh"

int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		std::cerr << "ERROR: Incorrect use of " << argv[0] << std::endl;
		std::cerr << "Correct use: " << argv[0] << " filename" << std::endl;
		return 1;
	}
	
	std::ifstream inputFile(argv[1]);
	if(!inputFile)
	{
		std::cerr << "Could not find file " << argv[1] << std::endl;
		return 1;
	}


	Vector<Vector<double>> positions;
	Vector<double> charges;
	double Z, posx, posy, posz;
	while(inputFile >> Z >> posx >> posy >> posz)
	{
		charges.Append(Z);

		Vector<double> position(3);
		position[0] = posx;
		position[1] = posy;
		position[2] = posz;
		positions.Append(position);
	}

	int maxL = 0;
	int maxU = std::pow(2,12);
	BoysFunction boyFn = BoysFunction(4 * maxL, maxU, 8*maxU);
	Molecule mol = Molecule(2, charges, positions, maxL, boyFn, "./OrbitalData/STO3/");
	Vector<double> axis(300);
	Vector<double> energy(300);

	std::ofstream outfile("./OutputData/H2l.dat");
	for(int i = 300; i > 1; i--)
	{
		axis[300-i] = (double)(i) / 20.0;
	}
//	energy = mol.CalculateGroundCurve(axis);
	energy = mol.CalculateExcitedCurve(axis, 1);
	for (int i = 0; i < 300; ++i)
	{
		outfile << axis[i] << "\t" << energy[i] << "\n";
	}
	std::cout << "end" << std::endl;
	
	return 0;
}
