#include<iostream>
#include<string>
#include "Array.h" 


using namespace std;


Array::Array()//Empty constructor whitch put null in msa_arr.
{
row=0;
double_matrix=NULL;
}


Matrix::Matrix(int n , int m){//Implementation of a costructor which recieve parmater n as the row of msa_arr.
	row=n;
	column=m;
	double_matrix=new double *[row];
	for (int i=0; i<row; i++)

		double_matrix[i]=new double[column];
		for(int j=0; j<column ; j++)
		{
			double_matrix.setij(i,j,0.0);
		}
	}	
}



Matrix::Matrix(const Matrix &src)//Implementation of  a  copy Constructor.
{
	row=src.row;
	column=src.column;
	double_matrix=new double*[src.column];
	for (int i=0; i<src.column; i++)
	{
		double_matrix[i]=new double[row];
		memcpy(double_matrix[i],src.double_matrix[i],src.row*sizeof(double));
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
		free_marix();
	
	row=rhs.row;
	double_matrix=new double*[row];
	for (int i=0; i<rhs.row; i++)
	{
		double_matrix[i]=new double[column];
		memcpy(double_matrix[i],rhs.double_matrix[i],rhs.row*sizeof(msa));
	}
	return *this;
}	


ostream &operator <<(ostream &stream, msa_matrix src)//Function which overload '<<'.
{
	src.print();
	return stream;
}

void msa_matrix:: print()//Function which print the matrix values.
{
	msa temp;
	
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			temp=getij(i,j);
			cout<<temp.aa.c_str()<<" "<<temp.coorZ<<" "<<temp.hydro<<"   ,   ";
		}
	cout<<endl;
	cout<<endl;
	cout<<endl;
	}
}

msa msa_matrix:: getij(int i, int j)const//Get the value of the place ij in the matrix.
{
	if(msa_arr==NULL)
		cout <<"msa_arr == NULL" <<endl;
	else 
		if(i<0||i>row )
			cout <<"the index "<<i<<" mast be between 0-"<< row<<endl ;
		else if (j<0||j>12)
				cout <<"the index "<<j<<" mast be between 0-"<< row<<endl;
			else				
				return msa_arr[i][j];

			
}

string msa_matrix:: getijaa(int i, int j)const//Get the value of the place ij in the matrix.
{
	if(msa_arr==NULL)
		cout <<"msa_arr == NULL" <<endl;
	else 
		if(i<0||i>row )
			cout <<"the index "<<i<<" mast be between 0-"<< row<<endl ;
		else if (j<0||j>12)
				cout <<"the index "<<j<<" mast be between 0-"<< row<<endl;
			else				
				return msa_arr[i][j].aa;
		return NULL;
}

msa** msa_matrix:: getmsa_arr()//Get the value of the place ij in the matrix.
{
	if(msa_arr==NULL)
		cout <<"msa_arr == NULL" <<endl;
	else				
		return msa_arr;
	return NULL;
}

