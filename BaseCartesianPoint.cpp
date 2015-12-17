#include <iostream>
#include <iomanip>
#include"header\BaseCartesianPoint.h"
#include "header\CartesianGeometry.h"
#include "header\GlobalFunction.h"

using namespace std;
using namespace MSL;


BaseCartesianPoint::BaseCartesianPoint()//Implementation of a costructor which recieve parmater n as the row of msa_arr.
{
	
	//angle = 0.0;
	//distance = 0.0;
	avgZ1 = 0.0;
	avgZ2 = 0.0;

	msaPoints1 = new CartesianPoint[12];
	msaPoints2 = new CartesianPoint[12];
	center1 = new CartesianPoint[9];
	center2 = new CartesianPoint[9];
	for (int i = 0; i<12; i++)
	{
		msaPoints1[i].setX(0.0);
		msaPoints1[i].setY(0.0);
		msaPoints1[i].setZ(0.0);

		msaPoints2[i].setX(0.0);
		msaPoints2[i].setY(0.0);
		msaPoints2[i].setZ(0.0);
		
	}
	for (int i = 0; i<9; i++)
	{
		center1[i].setX(0.0);
		center1[i].setY(0.0);
		center1[i].setZ(0.0);

		center2[i].setX(0.0);
		center2[i].setY(0.0);
		center2[i].setZ(0.0);

	}
}



BaseCartesianPoint::BaseCartesianPoint(const BaseCartesianPoint &src)//Implementation of  a  copy Constructor.
{
	
	//angle = src.angle;
	//distance = src.distance;
	avgZ1 = src.avgZ1;
	avgZ2 = src.avgZ2;
	center1 = new CartesianPoint[9];
	center2 = new CartesianPoint[9];
	msaPoints1 = new CartesianPoint[12];
	msaPoints2 = new CartesianPoint[12];
	for (int i = 0; i<12; i++)
	{

		memcpy(&msaPoints1[i], &src.msaPoints1[i], sizeof(CartesianPoint));
		memcpy(&msaPoints2[i], &src.msaPoints2[i], sizeof(CartesianPoint));
	}

	for (int i = 0; i<9; i++)
	{

		memcpy(&center1[i], &src.center1[i], sizeof(CartesianPoint));
		memcpy(&center2[i], &src.center2[i], sizeof(CartesianPoint));
	}
}


void BaseCartesianPoint::freeBaseCartesianPoint()//Function which free dynamic allocated memory space.
{
	for (int i = 0; i < 12; i++)
	{
		delete[] msaPoints1;
		delete[] msaPoints2;
	}
	for (int i = 0; i < 9; i++)
	{
		delete[] center1;
		delete[] center2;
	}
}

BaseCartesianPoint::~BaseCartesianPoint()//Destructor which call free_matrix to deallocate dynamic memory space.
{
	freeBaseCartesianPoint();
}


BaseCartesianPoint& BaseCartesianPoint::operator=(const BaseCartesianPoint &rhs)//Function which overload '=' .
{

	if (this == &rhs)
		return *this;

	if (msaPoints1 != NULL || msaPoints2 != NULL)
		freeBaseCartesianPoint();

	//
	//angle = rhs.angle;
	//distance = rhs.distance;
	avgZ1= rhs.avgZ1;
	avgZ2 = rhs.avgZ2;
	msaPoints1 = new CartesianPoint[12];
	msaPoints2 = new CartesianPoint[12];
	center1 = new CartesianPoint[9];
	center2 = new CartesianPoint[9];
	for (int i = 0; i<12; i++)
	{

		memcpy(&msaPoints1[i], &rhs.msaPoints1[i], sizeof(CartesianPoint));
		memcpy(&msaPoints2[i], &rhs.msaPoints2[i], sizeof(CartesianPoint));
	}
	for (int i = 0; i<9; i++)
	{

		memcpy(&center1[i], &rhs.center1[i], sizeof(CartesianPoint));
		memcpy(&center2[i], &rhs.center2[i], sizeof(CartesianPoint));
	}
	return *this;
}

CartesianPoint* BaseCartesianPoint::getMsaPoints1(int i)const
{
	return &msaPoints1[i];
}
CartesianPoint* BaseCartesianPoint::getMsaPoints2(int i)const
{
	return &msaPoints2[i];
}
CartesianPoint* BaseCartesianPoint::getCenter1(int i)const
{
	return &center1[i];
}
CartesianPoint* BaseCartesianPoint::getCenter2(int i)const
{
	return &center2[i];
}

double BaseCartesianPoint::getRmsd()const
{
	return this->rmsd;
}

//double BaseCartesianPoint::getDistance()const
//{
//	return this->distance;
//}
double BaseCartesianPoint::getMsaX1(int i)const
{
	return  msaPoints1[i].getX();
}
double BaseCartesianPoint::getMsaX2(int i)const
{
	return  msaPoints2[i].getX();
}
double BaseCartesianPoint::getMsaY1(int i)const
{
	return  msaPoints1[i].getY();
}
double BaseCartesianPoint::getMsaY2(int i)const
{
	return  msaPoints2[i].getY();
}
double BaseCartesianPoint::getMsaZ1(int i)const
{
	return  msaPoints1[i].getZ();
}
double BaseCartesianPoint::getMsaZ2(int i)const
{
	return  msaPoints2[i].getZ();
}
double BaseCartesianPoint::getAvgZ1()const
{
	return this->avgZ1;
}



double BaseCartesianPoint::getAvgZ2()const
{
	return this->avgZ2;
}



void BaseCartesianPoint::setRmsd(double rmsd)
{
	this->rmsd = rmsd;
}

//void BaseCartesianPoint::setDistance(double distance)
//{
//	this->distance = distance;
//}

void BaseCartesianPoint::setAvgZ1(double avg)
{
	this->avgZ1 = avg;
}



void BaseCartesianPoint::setAvgZ2(double avg)
{
	this->avgZ2 = avg;
}



void BaseCartesianPoint::setPoint1(int i, double coorX, double coorY, double coorZ)
{
	msaPoints1[i].setX(coorX);
	msaPoints1[i].setY(coorY);
	msaPoints1[i].setZ(coorZ);
}

void BaseCartesianPoint::setPoint2(int i, double coorX, double coorY, double coorZ)
{
	msaPoints2[i].setX(coorX);
	msaPoints2[i].setY(coorY);
	msaPoints2[i].setZ(coorZ);
}

void BaseCartesianPoint::setMsaPoints(msa_matrix msa1, msa_matrix msa2, int i)
{
	
	for (int j = 0; j < 12; j++)
	{
		setPoint1(j,msa1.getijCoorX(i, j), msa1.getijCoorY(i, j), msa1.getijCoorZ(i, j));
		setPoint2(j,msa2.getijCoorX(i, j), msa2.getijCoorY(i, j), msa2.getijCoorZ(i, j));
	}
}

void BaseCartesianPoint::setCenter()
{
	Helanal temp;
	for (int j = 0; j < 9; j++)
	{
		Helanal temp1(msaPoints1[j], msaPoints1[j + 1], msaPoints1[j + 2], msaPoints1[j + 3]);
		Helanal temp2(msaPoints2[j], msaPoints2[j + 1], msaPoints2[j + 2], msaPoints2[j + 3]);
		center1[j] = temp1.getCenter();
		center2[j] = temp2.getCenter();
	}
}


void BaseCartesianPoint::calculatAvgZ()
{
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (int i = 0; i < 12; i++)
	{
		sum1 = sum1 + msaPoints1[i].getZ();
		sum2 = sum2 + msaPoints2[i].getZ();
	}
	sum1 = sum1 / 12.0;
	sum2 = sum2 / 12.0;
	setAvgZ1(sum1);
	setAvgZ2(sum2);
}


//void BaseCartesianPoint::calculatAngle()
//{
//	double alpha = CartesianGeometry::dihedral(*getMsaPoints1(0), *getMsaPoints1(1), *getMsaPoints2(0), *getMsaPoints2(1));
//	setAngle(alpha);
//}

void BaseCartesianPoint::print(int clN, BaseCartesianPoint *cluster_arr, int lenght)
{
	GlobalFunction globalFunc;
	ofstream comp;
	string fileName = "Helanal " + to_string(clN) + ".txt";
	comp = globalFunc.creat_file_format_txt(fileName);
	comp << "--------------------------------------- Helanal for cluster " << clN << "------------------------------------------------------" << endl;
	for (int i = 0; i < lenght; i++)
	{
		comp << "------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		
		comp << "******************helix 1 cartesian points******************" << endl;
		for (int j = 0; j < 12; j++)
		{
			comp << "| x= " << cluster_arr[i].getMsaPoints1(j)->getX() << "	y= " << cluster_arr[i].getMsaPoints1(j)->getY() << " z= "
				<< cluster_arr[i].getMsaPoints1(j)->getZ() << " ";
		}
		comp << endl;
		comp << "******************helix 2 cartesian points******************" << endl;
		for (int j = 0; j < 12; j++)
		{
			comp << "| x= " << cluster_arr[i].getMsaPoints2(j)->getX() << "	y= " << cluster_arr[i].getMsaPoints2(j)->getY() << " z= "
				<< cluster_arr[i].getMsaPoints2(j)->getZ() << " ";
		}
		comp << endl;
		comp << "******************helix 1 CENTER cartesian points******************" << endl;
		for (int j = 0; j < 9; j++)
		{
			comp << "| x= " << cluster_arr[i].getCenter1(j)->getX() << "	y= " << cluster_arr[i].getCenter1(j)->getY() << " z= "
				<< cluster_arr[i].getCenter1(j)->getZ() << " ";
		}
		comp << endl;
		comp << "******************helix 2 CENTER cartesian points******************" << endl;
		for (int j = 0; j < 9; j++)
		{
			comp << "| x= " << cluster_arr[i].getCenter2(j)->getX() << " y= " << cluster_arr[i].getCenter2(j)->getY() << " z= "
				<< cluster_arr[i].getCenter2(j)->getZ() << " ";
		}
		comp << endl;
		//comp << "distace:" << endl;
		//comp << cluster_arr[i].getDistance() << endl;
		comp << endl;
		comp << "rmsd:  " << cluster_arr[i].getRmsd() << endl;
		comp <<  endl;
		comp << "average Z1:  " << cluster_arr[i].getAvgZ1() << endl;
		comp <<  endl;
		comp << "average Z2:  " << cluster_arr[i].getAvgZ2() << endl;
	}
}

void BaseCartesianPoint::calculatRmsd()
{
	double total = 0.0;
	
	for (int i= 0; i < 12; i++)
	{
		total = total + (pow(getMsaX1(i) - getMsaX2(i), 2)) + (pow(getMsaY1(i) - getMsaY2(i), 2)) + (pow(getMsaZ1(i) - getMsaZ2(i), 2));
	}
	setRmsd(sqrt(total / 12));
	
}

void BaseCartesianPoint::printPDB(ofstream& pdbFile, vector<string> tmData, ofstream& pdbfullCluster)
{
	GlobalFunction globalFunc;
	vector<string> csvAtomNumber1 = globalFunc.split(tmData[7], "-");
	int csvAtomNumberStart1 = atoi(csvAtomNumber1[0].c_str());




	vector<string> csvAtomNumber2 = globalFunc.split(tmData[9], "-");
	int csvAtomNumberStart2 = atoi(csvAtomNumber2[0].c_str());


	for (int i = 0; i < 9; i++)
	{
		pdbFile << left << setw(6) << "ATOM" << setw(5) << 1 << setw(1) << " " << left << setw(5) << "CA" << setw(3) << "xxx " << setw(1) << tmData[6] << setw(4) << csvAtomNumberStart1 << setw(4) << " " << setw(8) << setprecision(3) << fixed << getCenter1(i)->getX() << setw(8) << setprecision(3) << fixed << getCenter1(i)->getY() << setw(8) << setprecision(3) << fixed << getCenter1(i)->getZ() << endl;
		pdbfullCluster << left << setw(6) << "ATOM" << setw(5) << 1 << setw(1) << " " << left << setw(5) << "CA" << setw(3) << "xxx " << setw(1) << tmData[6] << setw(4) << csvAtomNumberStart1 << setw(4) << " " << setw(8) << setprecision(3) << fixed << getCenter1(i)->getX() << setw(8) << setprecision(3) << fixed << getCenter1(i)->getY() << setw(8) << setprecision(3) << fixed << getCenter1(i)->getZ() << endl;

		//pdbfullCluster << "ATOM  " << 1 << "     CA   xxx " << tmData[6] << "   " << csvAtomNumberStart1  << "  " << setprecision(3) << fixed << getCenter1(i)->getX() << "  " << setprecision(3) << fixed << getCenter1(i)->getY() << "  " << setprecision(3) << fixed << getCenter1(i)->getZ() << endl;
		csvAtomNumberStart1++;
	}

	for (int i = 0; i < 9; i++)
	{
		pdbFile << left << setw(6) << "ATOM" << setw(5) << 2 << setw(1) << " " << left << setw(5) << "CA" << setw(3) << "xxx " << setw(1) << tmData[8] << setw(4) << csvAtomNumberStart2 << setw(4) << " " << setw(8) << setprecision(3) << fixed << getCenter2(i)->getX() << setw(8) << setprecision(3) << fixed << getCenter2(i)->getY() << setw(8) << setprecision(3) << fixed << getCenter2(i)->getZ() << endl;
		pdbfullCluster << left << setw(6) << "ATOM" << setw(5) << 2 << setw(1) << " " << left << setw(5) << "CA" << setw(3) << "xxx " << setw(1) << tmData[8] << setw(4) << csvAtomNumberStart2 << setw(4) << " " << setw(8) << setprecision(3) << fixed << getCenter2(i)->getX() << setw(8) << setprecision(3) << fixed << getCenter2(i)->getY() << setw(8) << setprecision(3) << fixed << getCenter2(i)->getZ() << endl;

		//pdbFile        << "ATOM      " << 2 << "  CA  xxx " << tmData[8] << "  " << csvAtomNumberStart2 << "-" << csvAtomNumberStart2 + 3 << "  " << setprecision(3) << fixed << getCenter2(i)->getX() << "  " << setprecision(3) << fixed << getCenter2(i)->getY() << "  " << setprecision(3) << fixed << getCenter2(i)->getZ() << endl;
		//pdbfullCluster << "ATOM      " << 2 << "  CA  xxx " << tmData[8] << "  " << csvAtomNumberStart2 << "-" << csvAtomNumberStart2 + 3 << "  " << setprecision(3) << fixed << getCenter1(i)->getX() << "  " << setprecision(3) << fixed << getCenter1(i)->getY() << "  " << setprecision(3) << fixed << getCenter1(i)->getZ() << endl;
		csvAtomNumberStart2++;
	}
	
}