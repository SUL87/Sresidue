#ifndef MSA_H
#define MSA_H


#include <iostream>
#include<vector>
using namespace std;

	struct msa
	{
		string aa;
		double coorX;
		double coorY;
		double coorZ;
		double hydro;
	};
	class msa_matrix
	{
	private:

		int row;
		int column;
		msa **msa_arr;

	public:

		msa_matrix();//Empty Constructor
		msa_matrix(int n);//Constructor
		msa_matrix(const msa_matrix &src);//Copy constructor
		~msa_matrix();//Destructor
		msa_matrix& operator=(const msa_matrix &rhs);//Function which overload '=' .
		void setij(string aa, double coorX, double coorY, double coorZ, double hydro, int i, int j);//Function which sets values in the matrix.
		msa getij(int i, int j)const;
		string getijaa(int i, int j)const;
		double getijCoorZ(int i, int j)const;
		double getijCoorY(int i, int j)const;
		double getijCoorX(int i, int j)const;
		double getijHydro(int i, int j)const;
		void print(ofstream &stream);
		void free_msa();
		msa** getmsa_arr();
		void hydroZcoorAvg(ofstream&  myfile, int count);
		void init(vector<string>& resSplit, vector<string>& tmData, string flag, int i, ofstream& msaErrorFile);
		void standard_deviation(vector<double>& avgZ, vector<double>& avgHydro, ofstream&  myfile);

	};



#endif
