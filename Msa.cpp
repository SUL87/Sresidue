#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "header\Msa.h"
#include "header\GlobalFunction.h"

using namespace std;

msa_matrix::msa_matrix()//Empty constructor whitch put null in msa_arr.
{
row=0;
column=0;
msa_arr=NULL;
}


msa_matrix::msa_matrix(int n){//Implementation of a costructor which recieve parmater n as the row of msa_arr.
	row=n;
	column = 12;
	msa_arr = new msa*[row];
	for (int i=0; i<row; i++)
	{
		msa_arr[i] = new msa[column];
		for(int j=0; j<column ; j++)
		{
			msa_arr[i][j].aa = "";
			msa_arr[i][j].coorX = 0.0;
			msa_arr[i][j].coorY = 0.0;
			msa_arr[i][j].coorZ = 0.0;
			msa_arr[i][j].hydro = 0.0;
		}
	}	
}



msa_matrix::msa_matrix(const msa_matrix &src)//Implementation of  a  copy Constructor.
{
	row = src.row;
	column = src.column;
	msa_arr=new msa*[row];
	for (int i=0; i<row; i++)
	{
		msa_arr[i]=new msa[column];
		for (int j=0; j<column; j++)
		{
			memcpy(&msa_arr[i][j], &src.msa_arr[i][j],sizeof(msa));
		}
	}
}


void msa_matrix::free_msa(){//Function which free dynamic allocated memory space.
	for (int i=0; i<row; i++)
		delete 	[] msa_arr[i];
	delete msa_arr;
}

msa_matrix::~msa_matrix()//Destructor which call free_matrix to deallocate dynamic memory space.
{
	free_msa();
}


void msa_matrix::setij(string aa ,double coorX,double coorY,double coorZ,double hydro, int i, int j)//Function which sets values in the matrix.
{
	if(msa_arr==NULL)
		cout <<"msa_arr == NULL" <<endl;
	else 
		if(i<0||i>row )
			cout <<"the index "<<i<<" mast be between 0-"<< row<<endl ;
		else if (j<0||j>12)
				cout <<"the index "<<j<<" mast be between 0-12"<<endl;
			else
			{					
				msa_arr[i][j].aa=aa;
				msa_arr[i][j].coorX=coorX; //to change
				msa_arr[i][j].coorY=coorY;
				msa_arr[i][j].coorZ=coorZ;
				msa_arr[i][j].hydro=hydro;
			}
}

msa_matrix& msa_matrix::operator=(const msa_matrix &rhs)//Function which overload '=' .
{
	
	if(this==&rhs)
		return *this;
	
	if(msa_arr!=NULL)
		free_msa();
	
	row = rhs.row;
	column = rhs.column;
	msa_arr = new msa*[row];
	for (int i = 0; i<row; i++)
	{
		msa_arr[i] = new msa[column];
		for (int j = 0; j<column; j++)
		{
			memcpy(&msa_arr[i][j], &rhs.msa_arr[i][j], sizeof(msa));
		}
	}
	return *this;
}	


void msa_matrix::print(ofstream &stream)//Function which print the matrix values.
{
	msa temp;
	stream << "------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			temp=getij(i,j);
			stream << '|' << setw(10) << "msa.aa=" << temp.aa.c_str() << '|' << setw(10) << "  msa.coorX=" << temp.coorX << '|' << setw(10) << "  msa.coorY=" << temp.coorY << '|' << setw(10) << "  msa.coorZ=" << temp.coorZ << '|' << setw(10) << "  msa.hydro=" << temp.hydro << '|' << endl;
			 //stream << "msa.aa=" << temp.aa.c_str() << "  msa.coorZ=" << temp.coorZ << "  msa.hydro=" << temp.hydro << "   ,   ";
		}
		stream << "------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

	}
	
}

msa msa_matrix:: getij(int i, int j)const//Get the value of the place ij in the matrix.
{
	if(msa_arr==NULL)
		cout <<"msa_arr == NULL" <<endl;
	else 
		if(i<0||i>row )
			cout << "the index " << i << " mast be between 0-" << row << endl;
		else if (j<0||j>12)
			cout << "the index " << j << " mast be between 0-" << row << endl;
			else				
				return msa_arr[i][j];
		
}

string msa_matrix:: getijaa(int i, int j)const//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout <<"msa_arr == NULL" <<endl;
	else 
	if (i<0 || i>row)
		cout << "the index " << i << " mast be between 0-" << row << endl;
		else if (j<0||j>12)
			cout << "the index " << j << " mast be between 0-" << row << endl;
			else				
				return msa_arr[i][j].aa;
		return NULL;
}

double msa_matrix::getijCoorZ(int i, int j)const//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout << "msa_arr == NULL" << endl;
	else
	if (i<0 || i>row)
		cout << "the index " << i << " mast be between 0-" << row << endl;
	else if (j<0 || j>12)
		cout << "the index " << j << " mast be between 0-" << row << endl;
	else
		return msa_arr[i][j].coorZ;
	return NULL;
}

double msa_matrix::getijCoorY(int i, int j)const//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout << "msa_arr == NULL" << endl;
	else
	if (i<0 || i>row)
		cout << "the index " << i << " mast be between 0-" << row << endl;
	else if (j<0 || j>12)
		cout << "the index " << j << " mast be between 0-" << row << endl;
	else
		return msa_arr[i][j].coorY;
	return NULL;
}

double msa_matrix::getijCoorX(int i, int j)const//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout << "msa_arr == NULL" << endl;
	else
	if (i<0 || i>row)
		cout << "the index " << i << " mast be between 0-" << row << endl;
	else if (j<0 || j>12)
		cout << "the index " << j << " mast be between 0-" << row << endl;
	else
		return msa_arr[i][j].coorX;
	return NULL;
}

double msa_matrix::getijHydro(int i, int j)const//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout << "msa_arr == NULL" << endl;
	else
	if (i<0 || i>row)
		cout << "the index " << i << " mast be between 0-" << row << endl;
	else if (j<0 || j>12)
		cout << "the index " << j << " mast be between 0-" << row << endl;
	else
		return msa_arr[i][j].hydro;
	return NULL;
}

msa** msa_matrix::getmsa_arr()//Get the value of the place ij in the matrix.
{
	if (msa_arr == NULL)
		cout <<"msa_arr == NULL" <<endl;
	else				
		return msa_arr;
	return NULL;
}

void msa_matrix::init(vector<string>& resSplit, vector<string>& tmData, string flag, int i, ofstream& msaErrorFile)
{
	GlobalFunction globalFunc;


	if (resSplit[0] == "ATOM")
	{
		if (flag == "msaArr")
		{
			vector<string> csvAtomNumber = globalFunc.split(tmData[7], "-");
			int csvAtomNumberStart = atoi(csvAtomNumber[0].c_str());
			int csvAtomNumberEnd = atoi(csvAtomNumber[1].c_str());

			if (resSplit[4] == tmData[6])
			{
				if (strcmp(resSplit[2].c_str(), "CA") == 0)
				{
					int pdbAtomNumber = atoi(resSplit[5].c_str());
					if (pdbAtomNumber >= csvAtomNumberStart && pdbAtomNumber <= csvAtomNumberEnd)
					{
						setij(globalFunc.threeToOne(resSplit[3].c_str(), msaErrorFile),
							atof(resSplit[6].c_str()),//atodouble
							atof(resSplit[7].c_str()),
							atof(resSplit[8].c_str()),
							globalFunc.Hydrophobic(globalFunc.threeToOne(resSplit[3].c_str(), msaErrorFile)),
							i,
							(pdbAtomNumber - csvAtomNumberStart));
					}
				}
			}
		}
		else if (flag == "msaArr2")
		{
			vector<string> csvAtomNumber = globalFunc.split(tmData[9], "-");
			int csvAtomNumberStart = atoi(csvAtomNumber[0].c_str());
			int csvAtomNumberEnd = atoi(csvAtomNumber[1].c_str());
			if (resSplit[4] == tmData[8])
			{
				if (strcmp(resSplit[2].c_str(), "CA") == 0)
				{
					int pdbAtomNumber = atoi(resSplit[5].c_str());
					if (pdbAtomNumber >= csvAtomNumberStart && pdbAtomNumber <= csvAtomNumberEnd)
					{
						setij(globalFunc.threeToOne((resSplit[3]).c_str(), msaErrorFile),
							atof(resSplit[6].c_str()),
							atof(resSplit[7].c_str()),
							atof(resSplit[8].c_str()),
							globalFunc.Hydrophobic(globalFunc.threeToOne((resSplit[3]).c_str(), msaErrorFile)),
							i,
							(pdbAtomNumber - csvAtomNumberStart));
					}
				}
			}
		}
	}

	return;
}
void msa_matrix::hydroZcoorAvg(ofstream&  myfile, int count)
{
	vector<double> zCoorArr(12);
	vector<double> hydroArr(12);
	GlobalFunction g;
	for (int k = 0; k<12; k++)
	{
		for (int i = 0; i<count; i++)
		{
			zCoorArr[k] += msa_arr[i][k].coorZ;
			hydroArr[k] += msa_arr[i][k].hydro;
		}
		zCoorArr[k] = zCoorArr[k] / count;
		hydroArr[k] = hydroArr[k] / count;

	}
	myfile << "************************************************************* avg of Z coordinate *************************************************************" << endl;
	g.print_array(myfile, zCoorArr);
	myfile << "********************************************************** Avg of Hydrophobicity *************************************************************" << endl;
	g.print_array(myfile, hydroArr);
	standard_deviation(zCoorArr, hydroArr, myfile);
}


void msa_matrix::standard_deviation(vector<double>& avgZ, vector<double>& avgHydro, ofstream&  myfile)
{
	GlobalFunction g;
	vector<double> zCoorStd(12);
	vector<double> hydroStd(12);
	for (int k = 0; k<column; k++)//calculate STD
	{
		for (int i = 0; i<row; i++)
		{
			zCoorStd[k] = zCoorStd[k] + pow((getijCoorZ(i, k) - avgZ[k]), 2);
			hydroStd[k] = hydroStd[k] + pow((getijHydro(i, k) - avgHydro[k]), 2);
		}
		zCoorStd[k] = sqrt(zCoorStd[k] / (double)row);
		hydroStd[k] = sqrt(hydroStd[k] / (double)row);
	}
	myfile << "************************************************** standard deviation  for Z coordinate  ***************************************************** " << endl;
	g.print_array(myfile, zCoorStd);
	myfile << "************************************************* standard deviation  for Hydrophobicity ***************************************************** "  << endl;
	g.print_array(myfile, hydroStd);

}