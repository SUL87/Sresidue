#include <iostream>
#include <string>
#include "header\Matrix.h" 
#include "header\GlobalFunction.h"
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;


Matrix::Matrix()//Empty constructor whitch put null in msa_arr.
{
row=0;
column=0;
double_matrix=NULL;
}


Matrix::Matrix(int n , int m){//Implementation of a costructor which recieve parmater n as the row of msa_arr.


	column = n;
	row=m;
	

	double_matrix=new double *[row];
	for (int i=0; i<row; i++)
	{
		double_matrix[i]=new double[column];
		for(int j=0; j<column ; j++)
		{
			double_matrix[i][j]=0.0;
		}
	}	
}
Matrix::Matrix(int n, int m,double val){//Implementation of a costructor which recieve parmater n as the row of msa_arr.


	column = n;
	row = m;


	double_matrix = new double *[row];
	for (int i = 0; i<row; i++)
	{
		double_matrix[i] = new double[column];
		for (int j = 0; j<column; j++)
		{
			double_matrix[i][j] = val;
		}
	}
}



Matrix::Matrix(const Matrix &src)//Implementation of  a  copy Constructor.
{
	row=src.row;
	column=src.column;
	double_matrix=new double*[src.row];
	for (int i=0; i<src.row; i++)
	{
		double_matrix[i]=new double[column];
		memcpy(double_matrix[i],src.double_matrix[i],src.column*sizeof(double));
	}

	
}


void Matrix::free_matrix(){//Function which free dynamic allocated memory space.
	for (int i=0; i<row; i++)
		delete 	[] double_matrix[i];
	delete double_matrix;
}

Matrix::~Matrix()//Destructor which call free_matrix to deallocate dynamic memory space.
{
	free_matrix(); 
}


void Matrix::setij( int i, int j,double num)//Function which sets values in the matrix.
{
	if(double_matrix==NULL)
		cout <<"double_matrix == NULL" <<endl;
	else 
		if(i<0||i>row )
			cout <<"the index "<<i<<" mast be between 0-"<< row<<endl ;
		else if (j<0||j>12)
				cout <<"the index "<<j<<" mast be between 0-"<< row<<endl;
			else
			{					
				double_matrix[i][j]=num;
			}
}

Matrix& Matrix::operator=(const Matrix &rhs)//Function which overload '=' .
{
	
	if(this==&rhs)
		return *this;
	
	if(double_matrix!=NULL)
		free_matrix();
	
	row=rhs.row;
	column=rhs.column;
	double_matrix=new double*[column];
	for (int i=0; i<rhs.column; i++)
	{
		double_matrix[i]=new double[row];
		memcpy(double_matrix[i],rhs.double_matrix[i],rhs.row*sizeof(double));
	}
	return *this;
}	




double Matrix:: getij(int i, int j)const//Get the value of the place ij in the matrix.
{
	if(double_matrix==NULL)
		cout <<"double_matrix = NULL" <<endl;
	else 
		if(i<0||i>row )
			cout <<"the index "<<i<<" must be between 0-"<< row<<endl ;
		else if (j<0||j>12)
				cout <<"the index "<<j<<" must be between 0-"<< column<<endl; 
			else				
				return double_matrix[i][j];
}

int Matrix:: getRow()const//Get the value of the place ij in the matrix.
{
	return row;
}
int Matrix:: getColumn()const//Get the value of the place ij in the matrix.
{
	return column;
}







double** Matrix:: getMatrix()//Get the value of the place ij in the matrix.
{
	if(double_matrix==NULL)
		cout <<"double_matrix == NULL" <<endl;
	else				
		return double_matrix;
	return NULL;
}

void Matrix::print(ofstream &file)
{
	GlobalFunction g;
	
	for (int i = 0; i<column; i++)
	{
		file << '|' << setw(10) << (i + 1) << "	";
	}
	file << endl;

	for (int k = 0; k<row; k++)//Print to file the PSSM
	{
		for (int i = 0; i<column; i++)
		{
			file << "|" << setw(10) << getij(k, i) << "	";
		}
		file << "----> " << g.intToAminoAcid(k) << endl;
	}
}


void Matrix::getSumOfColumns(ofstream& file)
{
	vector<double> sum(12) ;
	GlobalFunction g;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			sum[j] =sum[j] + getij(i, j);
		}
	}
	g.print_array(file, sum);
}

Matrix Matrix::operator*(const Matrix & _m) const {
	if (column != _m.row) {
		cerr << "ERROR 9232: incorrect size for multiply matrices (" << row << ", " << column << ") (" << _m.row << ", " << _m.column << ") in Matrix Matrix::operator*(Matrix & _m) const" << endl;
		exit(9232);
	}
	Matrix out(row, _m.column, 0.0);
	for ( int i = 0; i<row; i++) {
		for ( int j = 0; j<_m.column; j++) {
			for ( int k = 0; k<_m.row; k++) {
				out.setij(i, j, out.getij(i, j) + getij(i, k) * _m.getij(k, j));
			}
		}
	}

	return out;
}

void Matrix::operator*=(const Matrix & _m) {
	Matrix product = *this * _m;
	*this = product;
}