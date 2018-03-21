#include "Molecule.hh"
#include "Vector.hh"
#include <cmath>
#include "STOnGOrbit.hh"
#include "Functions.hh"
#include "Matrix.hh"
#include "LinearAlgebra.hh"
#include <iostream>
#include <fstream>

int main()
{

	Vector<Vector<double>> positions(2);
	Vector<double> pos1(3), pos2(3);

	pos1[0] = 0.0;
    pos1[1] = 0.0;
    pos1[2] = 0.0;

    pos2[0] = 0.0;
    pos2[1] = 0.0;
	pos2[2] = 0.0;

	positions[0] = pos1;
	positions[1] = pos2;

	Vector<double> charges(2);
	charges[0] = 1;
	charges[1] = 1;
	std::ofstream outfile("./testData.dat");

	for(int i = 0; i < 100; i++)
	{
		positions[0][0] = 10.0 * (i /50.0);
		
		positions[0].Print();
		Molecule a = Molecule(charges, positions, 3);
		a.CalculateEnergy();
		outfile << positions[0][0] << "\t" << a.GetEnergyLevels()[0] << "\t" << a.GetEnergyLevels()[1] << "\n";
	}
	return 0;
}
