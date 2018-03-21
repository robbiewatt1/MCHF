#include "Molecule.hh"
#include "Vector.hh"
#include <cmath>
#include "STOnGOrbit.hh"
#include "Functions.hh"
#include "Matrix.hh"
#include "LinearAlgebra.hh"

int main()
{

	Vector<Vector<double>> positions(2);

	Vector<double> pos1(3), pos2(3);

	pos1[0] = 0;
    pos1[1] = 0.0;
    pos1[2] = 0.0;

    pos2[0] = -0;
    pos2[1] = 2.05;
	pos2[2] = 0.0;

	positions[0] = pos1;
	positions[1] = pos2;

	Vector<double> charges(2);
	charges[0] = 1;
	charges[1] = 1;
	
	Molecule a = Molecule(charges, positions, 3);

	a.CalculateEnergy();

	std::cout << a.GetEnergyLevels()[0] << std::endl;
	std::cout << a.GetEnergyLevels()[1] << std::endl;


	return 0;
}
