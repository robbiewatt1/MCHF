
#include "Orbital.hh"
#include "H5Cpp.h"
#include "Vector.hh"
#include "Array3D.hh"

Orbital::Orbital()
{
}

Orbital::Orbital(const Vector<double> &xAxis, const Vector<double> &yAxis,
                 const Vector<double> &zAxis):
	m_xAxis(xAxis), m_yAxis(yAxis), m_zAxis(zAxis)
{
	SetData(xAxis, yAxis, zAxis);
}

Orbital::Orbital(const Orbital &orbital)
{
	m_xAxis = orbital.m_xAxis;
	m_yAxis = orbital.m_yAxis;
	m_zAxis = orbital.m_zAxis;
	m_data  = orbital.m_data; 
}

Orbital::~Orbital()
{
}

void Orbital::SetXAxis(const Vector<double> &xAxis)
{
	m_xAxis = xAxis;
}

void Orbital::SetYAxis(const Vector<double> &yAxis)
{
	m_yAxis = yAxis;
}

void Orbital::SetZAxis(const Vector<double> &zAxis)
{
	m_zAxis = zAxis;
}

void Orbital::SetData(const Vector<double> &xAxis, const Vector<double> &yAxis,
                      const Vector<double> &zAxis )
{
	m_data = Array3D<double>(xAxis.Length(), yAxis.Length(), zAxis.Length());
}

Vector<double> Orbital::GetXAxis() const
{
	return m_xAxis;
}

Vector<double> Orbital::GetYAxis() const
{
	return m_yAxis;
}

Vector<double> Orbital::GetZAxis() const
{
	return m_zAxis;
}

Array3D<double> Orbital::GetData() const
{
	return m_data;
}

void Orbital::OutputData()
{
	const H5std_string FILE_NAME("OrbitalData.h5");
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
				data_buff[i][j][k] = m_data[i][j][k];
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


/*
Orbital Orbital::operator=(const Orbital &orbital)
{
	// self assigment guard
	if(this == &orbital)
	{
		return *this;
	}

	delete[] m_data;

	m_xAxis = orbital.m_xAxis;
	m_yAxis = orbital.m_yAxis;
	m_zAxis = orbital.m_zAxis;
	if (orbital.m_data)
	{
		m_data = new double **[m_xAxis.Length()];
		for (int i = 0; i < m_xAxis.Length(); i++)
		{
			m_data[i] = new double *[m_yAxis.Length()];
			for (int j = 0; j < m_yAxis.Length(); j++)
			 {
			 	m_data[i][j] = new double [m_zAxis.Length()];
			 }
		}
	} else
	{
		m_data = 0;
	}
	return *this;
}
*/