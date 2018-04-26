#ifndef MATRIX_HH
#define MATRIX_HH

#include <iostream>
#include <cassert>
#include <iomanip>

template <typename T>
class Matrix
{
	private:
		int m_nRow;
		int m_nColumn;
		T **m_data;
	
	public:

		Matrix():
		m_nRow(0),m_nColumn(0),m_data(nullptr)
		{
		}

		Matrix(int nRow, int nColumn, bool identity=false):
		m_nRow(nRow),m_nColumn(nColumn)
		{
			// allocate the array space
			m_data = new T *[m_nRow];
			for(int i = 0; i < m_nRow; i++)
			{
				m_data[i] = new T [m_nColumn];
			}

			if(identity)
			{
				// Fill the array as idenity
				for(int i = 0; i < m_nRow; i++)
				{
					for(int j = 0; j < m_nColumn; j++)
					{
						if(i == j)
						{
							m_data[i][j] = 1.0;
						}else
						{
							m_data[i][j] = 0.0;
						}
					}
				}
			}else
			{
				// Fill the array with Zeroes
				for(int i = 0; i < m_nRow; i++)
				{
					for(int j = 0; j < m_nColumn; j++)
					{
						m_data[i][j] = 0;
					}
				}
			}
		}

		// Defult destructor. This deletes the heap memory and avoids memory leaks
		~Matrix()
		{
			delete[] m_data;
		}

		Matrix<T> DeepCopy() const
		{
			Matrix<T> copy(m_nRow, m_nColumn);
			for (int i = 0; i < m_nRow; i++)
			{
				for (int j = 0; j < m_nColumn; j++)
				{
					copy.m_data[i][j] = m_data[i][j];
				}
			}
			return copy;
		}
		
		int GetRows() const
		{
			return m_nRow;
		}

		int GetColumns() const
		{
			return m_nColumn;
		}

		// Prints Matrix to the screen
		void Print() const
		{
			std::cout << "[ ";
			for(int i =0; i < m_nRow; i++)
			{
				for(int j = 0; j < m_nColumn; j++)
				{
					std::cout << std::fixed << std::setprecision(4) << m_data[i][j] << ", ";
				}
				if(i < m_nRow -1)
				{
					std::cout << "\n  ";	
				}else
				{
					std::cout << "]" << std::endl;
				}
			}
		}


	public:

		// return the value of the matrix at point (rowIndex,columnIndex). Much prefered over 
		// [rowIndex,columnIndex] as preforms checks that index is in bound.
		T& operator()(const int rowIndex, const int columnIndex)
		{
			assert(rowIndex >= 0 && rowIndex < m_nRow);
			assert(columnIndex >= 0 && columnIndex < m_nColumn);

			return m_data[rowIndex][columnIndex];
		}

		// return the column at the given row index. This method is no recomeneded as the are 
		// no proper checks to determine if the index is out of range of the array
		T* operator[](const int rowIndex)
		{
			return m_data[rowIndex];
		}

		T* operator[](const int rowIndex) const
		{
			return m_data[rowIndex];
		}

		Matrix<T>& operator=(const Matrix<T> &matrix)
		{
			// Self assigmnet gard
			if( this == &matrix)
			{
				return *this;
			}

			// Delete any data in new matrix is holding
			delete[] m_data;

			// Copy variables over
			m_nRow = matrix.m_nRow;
			m_nColumn = matrix.m_nColumn;
			if(matrix.m_data)
			{
				// allocate the array space
				m_data = new T *[m_nRow];
				for(int i = 0; i < m_nRow; i++)
				{
					m_data[i] = new T [m_nColumn];
				}

				// Copy the data over
				for(int i = 0; i < m_nRow; i++)
				{
					for(int j = 0; j < m_nColumn; j++)
					{
						m_data[i][j] = matrix.m_data[i][j];
					}
				}
			}else
			{
				m_data=0;
			}

			return *this;
			// Copys a matirx to a new variable
		}

	friend Matrix<T> operator*(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
	{
		Matrix<T> newMatrix = Matrix<T>(matrix1.m_nRow,matrix2.m_nColumn);

		for(int i = 0; i < matrix1.m_nRow; i++)
		{
			for(int j = 0; j < matrix2.m_nColumn; j++)
			{
			double elementValue = 0.0;
			for(int k = 0; k < matrix1.m_nColumn; k++)
				{
					elementValue = elementValue + matrix1.m_data[i][k] * matrix2.m_data[k][j];
				}
			newMatrix(i,j) = elementValue;
			}
		}
		return newMatrix;
	}

};
#endif