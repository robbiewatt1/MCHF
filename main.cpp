#include "Molecule.hh"
#include "Vector.hh"
#include <cmath>
#include <iostream>
#include <fstream>

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
	}

	while(inputFile)
	{
		std::cout << inputFile << std::endl;
	}


	Vector<Vector<double>> positions(2);
	Vector<double> pos1(3),  pos2(3);

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
	std::ofstream outfile("./H2Plus.dat");

	for(int i = 0; i < 500; i++)
	{
		positions[0][0] = 20.0 * ((double)i /500.0);
		std::cout << i << std::endl;
		Molecule a = Molecule(charges, positions, 6);
		a.CalculateEnergy();
		outfile << positions[0][0] << "\t" << a.GetEnergyLevels()[0] << "\t" << a.GetEnergyLevels()[1] << "\t" <<
				   a.GetEnergyLevels()[2] << "\t" << a.GetEnergyLevels()[3] << "\t" << a.GetEnergyLevels()[4]  <<
				   "\t" << a.GetEnergyLevels()[5] << "\t" << a.GetEnergyLevels()[6] << "\t" << a.GetEnergyLevels()[7] << "\n";
	}
	return 0;
}
