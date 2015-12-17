
#ifndef BaseCartesianPoint_H
#define BaseCartesianPoint_H


#include <iostream>
#include "Msa.h"
#include "CartesianPoint.h"

#include "Helanal.h"
using namespace std;




	class BaseCartesianPoint
	{
	private:
		CartesianPoint* msaPoints1;
		CartesianPoint* msaPoints2;
		CartesianPoint* center1;
		CartesianPoint* center2;
		//double angle;
		double rmsd;
		double avgZ1;
		double avgZ2;

	public:

	
		BaseCartesianPoint();//Constructor
		BaseCartesianPoint(const BaseCartesianPoint &src);//Copy constructor
		~BaseCartesianPoint();//Destructor
		void freeBaseCartesianPoint();
		BaseCartesianPoint& operator=(const BaseCartesianPoint &rhs);//Function which overload '=' .

		CartesianPoint* getMsaPoints1(int i)const;
		CartesianPoint* getMsaPoints2(int i)const;
		CartesianPoint* getCenter1(int i)const;
		CartesianPoint* getCenter2(int i)const;

		double getRmsd()const;
		//double getDistance()const;
		double getAvgZ1()const;
		double getAvgZ2()const;
		void setRmsd(double rmsd);
		//void setDistance(double distance);
		void setAvgZ1(double avg);
		void setAvgZ2(double avg);
		void setPoint1(int i, double coorX, double coorY, double coorZ);
		void setPoint2(int i, double coorX, double coorY, double coorZ);
		void setMsaPoints(msa_matrix msa1, msa_matrix msa2, int i);//Function which sets values in the matrix.
		void setCenter();
		void calculatAvgZ();
		double getMsaX1(int i)const;
		double getMsaX2(int i)const;
		double getMsaY1(int i)const;
		double getMsaY2(int i)const;
		double getMsaZ1(int i)const;
		double getMsaZ2(int i)const;

		//*****   to do   *****//
		//void calculatAngle();
		void calculatRmsd();
		void printPDB(ofstream& pdbFile, vector<string> tmData, ofstream& pdbfullCluster);
		//void calculatDistance();
		static void print(int i, BaseCartesianPoint *cluster_arr, int lenght);



	};




#endif
