
#ifndef MATRIX_H
#define MATRIX_H


#include <iostream>

using namespace std;

	class Matrix
	{
	protected:

		int row;
		int column;
		double **double_matrix;

	public:

		Matrix();//Empty Constructor
		Matrix(int n, int m);//Constructor
		Matrix(int n, int m, double val);//Constructor
		Matrix(const Matrix &src);//Copy constructor
		~Matrix();//Destructor
		Matrix& operator=(const Matrix &rhs);//Function which overload '=' .
		void setij(int i, int j, double num);//Function which sets values in the matrix.
		double getij(int i, int j)const;
		int  getRow()const;
		int  getColumn()const;
		void getSumOfColumns(ofstream& file);
		void free_matrix();
		double** getMatrix();
		void print(ofstream &file);
		Matrix operator*(const Matrix & _m) const;
		void operator*=(const Matrix & _m);


	};

#endif
