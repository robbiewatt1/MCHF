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

	int maxL = 1;
	int maxU = std::pow(2,16);
	BoysFunction boyFn = BoysFunction(4 * maxL, maxU, 8*maxU);
	std::ofstream outfile("./OutputData/H2+.dat");
	for(int i = 100; i >= 98; i--)
	{
		std::cout << i << std::endl;
		positions[0][0] = (double)i / 5.0;
		Molecule mol = Molecule(1, charges, positions, maxL, boyFn, "./OrbitalData/6311GSS/");
		double energy = mol.CalculateEnergy();
		outfile << positions[0][0] << "\t" << energy << "\n";
	}
	return 0;
}
