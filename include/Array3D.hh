#ifndef ARRAY3D_HH
#define ARRAY3D_HH

#include <iostream>
#include <cassert>
#include <iomanip>
#include <limits>

template<typename T>
class Array3D
{
private:
	int m_nElementsX;
	int m_nElementsY;
	int m_nElementsZ;
	T ***m_data;

public:
	// Defult Constructor
	Array3D():
		m_nElementsX(0), m_nElementsY(0), m_nElementsZ(0), m_data(nullptr)
	{
	}

	// Copy constructor
	Array3D(const Array3D &array)
	{
		m_nElementsX = array.m_nElementsX;
		m_nElementsY = array.m_nElementsY;
		m_nElementsZ = array.m_nElementsZ;

		m_data = new T **[array.m_nElementsX];
		for (int i = 0; i < m_nElementsX; i++)
		{
			m_data[i] = new T *[array.m_nElementsY];
			for (int j = 0; j < m_nElementsY; j++)
			{
				m_data[i][j] = new T[array.m_nElementsZ];
			}
		}

		for (int i = 0; i < m_nElementsX; i++)
		{
			for (int j = 0; j < m_nElementsY; j++)
			{
				for (int k = 0; k < m_nElementsZ; k++)
				{
					m_data[i][j][k] = array.m_data[i][j][k];
				}
			}
		}
	}

	// Constructs an empty vector of nElement length
	Array3D(int nElementsX, int nElementsY, int nElementsZ):
		m_nElementsX(nElementsX), m_nElementsY(nElementsY), m_nElementsZ(nElementsZ)
	{
		// Allocate the vecotr memenory
		m_data = new T **[m_nElementsX];
		for (int i = 0; i < m_nElementsX; i++)
		{
			m_data[i] = new T *[m_nElementsY];
			for (int j = 0; j < m_nElementsY; j++)
			{
				m_data[i][j] = new T[m_nElementsZ];
			}
		}

		// Check if T is double / int /complex. If so then set all values to 0;
		for (int i = 0; i < m_nElementsX; i++)
		{
			for (int j = 0; j < m_nElementsY; j++)
			{
				for (int k = 0; k < m_nElementsZ; k++)
				{
					m_data[i][j][k] = 0;
				}
			}
		}
	}

	// Defult Destructor. Deletes the data to avoid memory leaks
	~Array3D()
	{
		delete [] m_data;
	}

	// Return the length of the Vector
	int Size() const
	{
		return (m_nElementsX * m_nElementsY * m_nElementsZ);
	}

	// Prints the vector to the screan
	void Print() const
	{
		for (int i = 0; i < m_nElementsX; i++)
		{
			for (int j = 0; j < m_nElementsY; j++)
			{
				for (int k = 0; k < m_nElementsZ; k++)
				{
					std::cout << std::fixed << std::setprecision(5) << m_data[i][j][k] << ",";
				}
				if (j < m_nElementsY - 1)
				{
					std::cout << "\n ";
				} else
				{
					std::cout << "\n\n ";
				}
			}
		}
	}

public:

	// returns the value of the vector at elementIndex. This method does not allow you to change the value of vector data.
	// This should be used when passing Vector to class as const.

	// returns the value of the vector at elementIndex. This method allows you to edit the vector data.
	T** operator[](const int elementIndex)
	{
		return m_data[elementIndex];
	}

	T** operator[](const int elementIndex) const
	{
		return m_data[elementIndex];
	}

	// Copy vecotor to new vector
	Array3D<T>& operator=(const Array3D<T> &array)
	{

		// Self assigment gaurd
		if (this == &array)
		{
			return *this;
		}

		// Delete data held by new vector
		delete[] m_data;

		// Copy variables over
		m_nElementsX = array.m_nElementsX;
		m_nElementsY = array.m_nElementsY;
		m_nElementsZ = array.m_nElementsZ;

		if (array.m_data)
		{
			// Allcoate new array space
			m_data = new T **[array.m_nElementsX];
			for (int i = 0; i < m_nElementsX; i++)
			{
				m_data[i] = new T *[array.m_nElementsY];
				for (int j = 0; j < m_nElementsY; j++)
				{
					m_data[i][j] = new T[array.m_nElementsZ];
				}
			}

			for (int i = 0; i < m_nElementsX; i++)
			{
				for (int j = 0; j < m_nElementsY; j++)
				{
					for (int k = 0; k < m_nElementsZ; k++)
					{
						m_data[i][j][k] = array.m_data[i][j][k];
					}
				}
			}
		} else
		{
			m_data = 0;
		}

		return *this;
	}

	/*
		//array-array operations:
		friend Array3D<T> operator+(const Array3D<T> &array1, const Array3D<T> &array2)
		{
			assert(array1.m_nElementsX == array2.m_nElementsX);
			assert(array1.m_nElementsY == array2.m_nElementsY);
			assert(array1.m_nElementsZ == array2.m_nElementsZ);

			Array3D<T> newArray = Array3D<T>(array1.m_nElements);

			for (int i = 0; i < vector1.m_nElements; i++)
			{
				newVector(i) = vector1.m_data[i] + vector2.m_data[i];
			}
			return newVector;
		};

		friend Vector<T> operator-(const Vector<T> &vector1, const Vector<T> &vector2)
		{
			assert(vector1.m_nElements == vector2.m_nElements);
			Vector<T> newVector = Vector<T>(vector1.m_nElements);

			for (int i = 0; i < vector1.m_nElements; i++)
			{
				newVector(i) = vector1.m_data[i] - vector2.m_data[i];
			}
			return newVector;
		};
	*/
	friend Array3D<T> operator*(const Array3D<T> &array1, const Array3D<T> &array2)
	{
		Array3D<T> newArray = Array3D<T>(array1.m_nElementsX, array1.m_nElementsY, array1.m_nElementsZ);

		for (int i = 0; i < array1.m_nElementsX; i++)
		{
			for (int j = 0; j < array1.m_nElementsY; j++)
			{
				for (int k = 0; k < array1.m_nElementsZ; k++)
				{
					newArray[i][j][k] = array1[i][j][k] * array2[i][j][k];
				}
			}
		}
		return newArray;
	};
	/*
		friend Vector<T> operator/(const Vector<T> &vector1, const Vector<T> &vector2)
		{
			assert(vector1.m_nElements == vector2.m_nElements);
			Vector<T> newVector = Vector<T>(vector1.m_nElements);

			for (int i = 0; i < vector1.m_nElements; i++)
			{
				newVector(i) = vector1.m_data[i] / vector2.m_data[i];
			}
			return newVector;
		};

		//T-vector operations:
		friend Vector<T> operator+(const T scalar, const Vector<T> &vector)
		{
			Vector<T> newVector = Vector<T>(vector.m_nElements);

			for (int i = 0; i < vector.m_nElements; i++)
			{
				newVector(i) = scalar + vector.m_data[i];
			}
			return newVector;
		};

		friend Vector<T> operator-(const T scalar, const Vector<T> &vector)
		{
			Vector<T> newVector = Vector<T>(vector.m_nElements);

			for (int i = 0; i < vector.m_nElements; i++)
			{
				newVector(i) = scalar - vector.m_data[i];
			}
			return newVector;
		};

		friend Vector<T> operator*(const T scalar, const Vector<T> &vector)
		{
			Vector<T> newVector = Vector<T>(vector.m_nElements);

			for (int i = 0; i < vector.m_nElements; i++)
			{
				newVector(i) = scalar * vector.m_data[i];
			}
			return newVector;
		};

		friend Vector<T> operator/(const T scalar, const Vector<T> &vector)
		{
			Vector<T> newVector = Vector<T>(vector.m_nElements);

			for (int i = 0; i < vector.m_nElements; i++)
			{
				newVector(i) = scalar / vector.m_data[i];
			}
			return newVector;
		};
	*/
};

#endif