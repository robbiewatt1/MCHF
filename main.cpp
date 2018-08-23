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
//	Vector<double> pos(3);
//	pos[0] = 0.0;
//	pos[1] = 0.0;
//	pos[2] = 0.0;
//	GaussianOrbital test = GaussianOrbital(0,0,0,1.0,pos);
//	GaussianOrbital test2 = GaussianOrbital(1,0,0,1.0,pos);
//	std::cout << test.MatrixElement(test2)[0] << " " << test.MatrixElement(test2)[1] << " " << test.MatrixElement(test2)[2] << std::endl;
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
	int maxU = 1000000;
	BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 100*maxU);
	std::ofstream outfile("./OutputData/testMe.dat");
	for(int i = 10; i >= 1; i--)
	{
		srand((unsigned int) time(0));
		positions[0][0] = (double)i / 10.0;
		Molecule mol = Molecule(charges, positions, maxL, boyFn);
		mol.CalculateEnergy();
		Vector<double> me = mol.MatrixElement(0,1);
		Vector<double> me2 = mol.MatrixElement(0,2);
		Vector<double> me3 = mol.MatrixElement(0,3);
		mol.GetEnergyLevels().Print();
		me.Print();
		me2.Print();
		me3.Print();
		outfile << positions[0][0] << "\t" << me[0] << "\t" << me[1] << "\t" << me[2] << std::endl;
		//outfile << positions[0][0] << "\t" << mol.GetEnergyLevels()[0] << "\t" << mol.GetEnergyLevels()[1] << "\t" << std::endl;
			//	<<  mol.GetEnergyLevels()[2] << "\t" << mol.GetEnergyLevels()[3] << "\t" << mol.GetEnergyLevels()[4]  
				//<< "\t" << mol.GetEnergyLevels()[5] << "\t" << mol.GetEnergyLevels()[6] 
				//<< "\t" << mol.GetEnergyLevels()[7] << "\n";
	//	std::cout <<positions[0][0] << " " << mol.GetEnergyLevels()[0] << " " << mol.GetEnergyLevels()[1] << " " << mol.GetEnergyLevels()[2] << mol.GetEnergyLevels()[3] << " " << mol.GetEnergyLevels()[4] << " " << mol.GetEnergyLevels()[5] << std::endl;

	}



/*
    int maxL = 2;
    int maxU = 1000000;
    BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 100*maxU);

    std::ofstream file("./test");
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 1000; ++j)
		{
			file << boyFn.Interpolate(i,j) << "\t";
		}
		file << "\n";
	}
*/	

	return 0;
}
