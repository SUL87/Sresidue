
#ifndef HELIX_H
#define HELIX_H


#include <iostream>

using namespace std;

	typedef struct helixStruct
	{
		string atom;
		string aa;// add all func 
		int sequence;//add all func   the number of the aa in the subunit sequence
		string subUnit;//add all func 
		double coorX;
		double coorY;
		double coorZ;
		double bFactor;
	}*helixP;

	class Helix
	{
	private:
		helixP helix_arr;
		int size;
	public:

		Helix();//Empty Constructor
		Helix(int n);//Constructor
		Helix(const Helix &src);//Copy constructor
		~Helix();//Destructor
		Helix& operator=(const Helix &rhs);//Function which overload '=' .
		void setHelix(string atom, double coorX, double coorY, double coorZ, int i);//Function which sets values in the matrix.
		helixStruct getHleixStruct(int i)const;
		string getHelixAtom(int i)const;
		double getHelixCoorX(int i)const;
		double getHelixCoorY(int i)const;
		double getHelixCoorZ(int i)const;
		void printHelix();
		void freeHelix();
		helixP getHelix_arr()const;
	};

	ostream& operator <<(ostream &stream, Helix s);//Function which overload '<<'.



#endif
