#include <cmath>
#include <omp.h>
#include <string>
#include "H5Cpp.h"
#include "boost/filesystem.hpp"

#include "Molecule.hh"
#include "LinearAlgebra.hh"
#include "STOnGOrbit.hh"
#include "Numerics.hh"
#include "FockSolver.hh"

Molecule::Molecule()
{
}

Molecule::Molecule(int nElectrons, const Vector<double> &nuclearCharges,
				   const Vector<Vector<double>> &nuclearPositions, int maxL,
				   const BoysFunction &boyFn, std::string basisSetDir):
m_nuclearCharges(nuclearCharges), m_nuclearPositions(nuclearPositions), m_maxL(maxL),
m_nElectrons(nElectrons)
{
	SetBasisSet(basisSetDir);
	m_boyFn = boyFn;
}

Molecule::~Molecule()
{
}

Vector<double> Molecule::GetEnergyLevels()
{
	return m_energyLevels;
}

Matrix<double> Molecule::GetBasisCoefficients()
{
	return m_basisSetCoefficients;
}

void Molecule::CalculateEnergy()
{

	FockSolver solver = FockSolver(m_nuclearCharges, m_nuclearPositions, m_basisSet,
								   m_nElectrons, m_boyFn);
	solver.Solve();
	/*
	Matrix<double> energyMaxtrix(m_basisSet.Length(), m_basisSet.Length());
	Matrix<double> overlapMatrix(m_basisSet.Length(), m_basisSet.Length());

	double potential(0);
	#pragma omp parallel for
	for (int i = 0; i < m_nuclearCharges.Length(); i++)
	{
		for (int j = i + 1; j < m_nuclearCharges.Length(); j++)
		{
			potential += CalculatePotential(m_nuclearCharges[i], m_nuclearCharges[j],
			                                m_nuclearPositions[i], m_nuclearPositions[j]);
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = i; j < m_basisSet.Length(); j++)
		{
			//Check omp is working
			//std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
			// loop over ion sights to get nuclear attraction integral
			double nuclearPotential(0);
			for (int k = 0; k < m_nuclearCharges.Length(); k++)
			{
				nuclearPotential += m_basisSet[i].NuclearOverlap(m_basisSet[j],
				                    m_nuclearCharges[k],
									m_nuclearPositions[k], m_boyFn);
			}
			double kineticEnergy = m_basisSet[i].KineticOverlap(m_basisSet[j]);
			overlapMatrix[i][j] = m_basisSet[i].Overlap(m_basisSet[j]);
			energyMaxtrix[i][j] = kineticEnergy + nuclearPotential + potential * overlapMatrix[i][j];

			// Matrix is symetric
			overlapMatrix[j][i] = overlapMatrix[i][j];
			energyMaxtrix[j][i] = energyMaxtrix[i][j];
		}
	}
	// When large number of basis sets are used, it is likely that the overlpa matrix is singular
	// Therefore this needs to be check and sorted
	Vector<double> overlapValues;
	Matrix<double> overlapVector;
	LinearAlgebra::EigenSolver(overlapMatrix, overlapVector, overlapValues);
	// eigen values that are smaller than 1e-12 are removed giving a new rectangle overlap matrix.
	double small = 1e-7 * overlapValues.End();
	Vector<double> reducedOVal;
	int dif(0);
	for (int i = 0; i <overlapValues.Length(); i++)
	{
		if(overlapValues[i] > small)
		{
			reducedOVal.Append(overlapValues[i]);
		} else
		{
			dif++;
		}
	}
	std::cout << "discared: " << dif << std::endl;
	Matrix<double> reducedOVec = Matrix<double>(m_basisSet.Length(), reducedOVal.Length());
	for (int i = 0; i < reducedOVec.GetRows(); i++)
	{
		for (int j = 0; j < reducedOVec.GetColumns(); j++)
		{
			reducedOVec[i][j] = overlapVector[i][j+dif] / std::sqrt(reducedOVal[j]);
		}
	}

	Matrix<double> reducedEnergy = LinearAlgebra::Transpose(reducedOVec) * energyMaxtrix * (reducedOVec);
	LinearAlgebra::EigenSolver(reducedEnergy,m_basisSetCoefficients, m_energyLevels);
	m_basisSetCoefficients = reducedOVec * m_basisSetCoefficients;

*/
}

Vector<double> Molecule::MatrixElement(int level1, int level2)
{
	Vector<double> matrixElement(3);
	if (m_energyLevels.Length() == 0)
	{
		std::cerr << "Error: Must first call CalculateEnergy to set energy levels" << std::endl;
		std::exit(-1);
	}
	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		for (int j = 0; j < m_basisSet.Length(); j++)
		{
			matrixElement = matrixElement + m_basisSet[i].MatrixElement(m_basisSet[j])
						  * (m_basisSetCoefficients[i][level1] * m_basisSetCoefficients[j][level2]);
		}
	}
	return matrixElement;
}

double Molecule::OscilatorStrength(int level1, int level2)
{
	Vector<double> matrixElement = MatrixElement(level1, level2);
	double mE_2 = matrixElement[0] * matrixElement[0] + matrixElement[1] * matrixElement[1] 
				+ matrixElement[2] * matrixElement[2];
	return (2.0 / 3.0) * mE_2 / (m_energyLevels[level2] - m_energyLevels[level1]);
}

Array3D<double> Molecule::CalculateWavefunction(int level)
{
	if (m_xAxis.Length() == 0)
	{
		std::cerr << "Error: xAxis must be set using SetXAxis() function" << std::endl;
	} else if (m_yAxis.Length() == 0)
	{
		std::cerr << "Error: yAxis must be set using SetXAxis() function" << std::endl;
	} else if (m_zAxis.Length() == 0)
	{
		std::cerr << "Error: zAxis must be set using SetXAxis() function" << std::endl;
	}
	Array3D<double> wavefunction = Array3D<double>(m_xAxis.Length(), m_yAxis.Length(), m_zAxis.Length());

	for (int i = 0; i < m_basisSet.Length(); i++)
	{
		Array3D<double> orbitalData =  m_basisSet[i].CalculateDataCartesian(m_xAxis, m_yAxis, m_zAxis);
		wavefunction = wavefunction + (m_basisSetCoefficients[i][level] * orbitalData);
	}
	return wavefunction;
}

void Molecule::SetXAxis(const Vector<double> &xAxis)
{
	m_xAxis = xAxis;
}

void Molecule::SetYAxis(const Vector<double> &yAxis)
{
	m_yAxis = yAxis;
}

void Molecule::SetZAxis(const Vector<double> &zAxis)
{
	m_zAxis = zAxis;
}

void Molecule::OutputData(int level, std::string fileName)
{
	const H5std_string FILE_NAME(fileName + ".h5");
	const H5std_string DATASET_NAMEX( "xAxis" );
	const H5std_string DATASET_NAMEY( "yAxis" );
	const H5std_string DATASET_NAMEZ( "zAxis" );
	const H5std_string DATASET_NAMED( "data" );

	const int vector_RANK = 1;
	const int data_RANK = 3;

	H5::Exception::dontPrint();
	H5::H5File* file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
	hsize_t dataDim[3];
	dataDim[0] = m_xAxis.Length();
	dataDim[1] = m_yAxis.Length();
	dataDim[2] = m_zAxis.Length();

	H5::DataSpace xSpace(vector_RANK, &dataDim[0]);
	H5::DataSpace ySpace(vector_RANK, &dataDim[1]);
	H5::DataSpace zSpace(vector_RANK, &dataDim[2]);
	H5::DataSpace dataSpace(data_RANK, dataDim);

	H5::DataSet* xSet = new H5::DataSet(
	    file->createDataSet(DATASET_NAMEX, H5::PredType::NATIVE_DOUBLE, xSpace));
	H5::DataSet* ySet = new H5::DataSet(
	    file->createDataSet(DATASET_NAMEY, H5::PredType::NATIVE_DOUBLE, ySpace));
	H5::DataSet* zSet = new H5::DataSet(
	    file->createDataSet(DATASET_NAMEZ, H5::PredType::NATIVE_DOUBLE, zSpace));
	H5::DataSet* dataset = new H5::DataSet(
	    file->createDataSet(DATASET_NAMED, H5::PredType::NATIVE_DOUBLE, dataSpace));
	
	// Need to convert data into contiguous array. ANNOYING! means having to double the mem :(
	Array3D<double> data = CalculateWavefunction(level);

	double x_buff[m_xAxis.Length()];
	double y_buff[m_yAxis.Length()];
	double z_buff[m_zAxis.Length()];
	double data_buff[m_xAxis.Length()][m_yAxis.Length()][m_zAxis.Length()];
	for (int i = 0; i < m_xAxis.Length(); i++)
	{
		x_buff[i] = m_xAxis[i];
		for (int j = 0; j < m_yAxis.Length(); j++)
		{
			y_buff[j] = m_yAxis[j];
			for (int k = 0; k < m_zAxis.Length(); k++)
			{
				z_buff[k] = m_zAxis[k];
				data_buff[i][j][k] = data[i][j][k];
			}
		}
	}
	xSet->write(x_buff, H5::PredType::NATIVE_DOUBLE);
	ySet->write(y_buff, H5::PredType::NATIVE_DOUBLE);
	zSet->write(z_buff, H5::PredType::NATIVE_DOUBLE);
	dataset->write(data_buff, H5::PredType::NATIVE_DOUBLE);

	delete xSet;
	delete ySet;
	delete zSet;
	delete dataset;
	delete file;
}

void Molecule::SetBasisSet(std::string basisSetDir)
{
	// Loop over all ions involved
	for (int i = 0; i < m_nuclearPositions.Length(); i++)
	{
		// First check to see if the directory exists.
		Vector<std::string> fileNames;
		std::string ionOrbitals = basisSetDir + std::to_string((int)m_nuclearCharges[i]);
		if(boost::filesystem::is_directory(ionOrbitals) == false)
		{
			std::cerr << "Error: " << ionOrbitals << " cannot be found" << std::endl;
			exit(-1); 
		} else
		{
			for(boost::filesystem::directory_iterator file = boost::filesystem::directory_iterator(ionOrbitals);
				file != boost::filesystem::directory_iterator(); ++file)
			{
				fileNames.Append(file->path().string());
			}
		}
		// Loop over k, m and n such that the sum is less than the maximum L
		for (int n = 0; n <= m_maxL; n++)
		{
			for (int m = 0; m <= m_maxL; m++)
			{
				for (int k = 0; k <= m_maxL; k++)
				{
					if (n + m + k <= m_maxL)
					{
						// Now loop through basis sets
						for (int j = 0; j < fileNames.Length(); ++j)
						{
							STOnGOrbit orbital = STOnGOrbit(fileNames[j], k, m, n,
										     m_nuclearPositions[i]);
							m_basisSet.Append(orbital);
						}

					}
				}
			}
		}
	}
}

double Molecule::CalculatePotential(int z1, int z2, Vector<double> ionLocation1,
                                    Vector<double> ionLocation2)
{
	double ionR = std::sqrt(std::pow(ionLocation1[0] - ionLocation2[0], 2.0)
	                        + std::pow(ionLocation1[1] - ionLocation2[1], 2.0)
	                        + std::pow(ionLocation1[2] - ionLocation2[2], 2.0));
	return z1 * z2 / ionR;
}
