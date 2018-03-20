#include "Molecule.hh"
#include "Vector.hh"


int main()
{
	Vector<Vector<double>> positions(2);

	Vector<double> pos1(3), pos2(3);

	pos1[0] = -1;
        pos1[1] = 0;
        pos1[2] = 0;

        pos2[0] = 1;
        pos2[1] = 0;
        pos2[2] = 0;

	positions[0] = pos1;
	positions[1] = pos2;

	Vector<double> charges(2);
	charges[0] = 0;
	charges[1] = 0;
	
	Molecule a = Molecule(charges, positions, 1);

	a.CalculateEnergy();

	a.GetEnergyLevels().Print();

	return 0;
}
