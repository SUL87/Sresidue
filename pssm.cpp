#include "header\Matrix.h"
#include "header\GlobalFunction.h"
#include "header\Msa.h"
#include <fstream>
#include <iomanip>
#include <array>
#include <math.h> 
#include <string>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include "header\pssm.h"
using namespace std;

Pssm::Pssm(int n,int type) : Matrix(12, 20)
{
	flag = type;
	cluster_row_size=n;
}



vector<double> Pssm::average()
{
	vector<double>avg(12);
	for (int k = 0; k<12; k++)//calculate the Mean
	{
		for (int i = 0; i<20; i++)
		{
			avg[k] += getij(i, k);
		}
		avg[k] = avg[k] / 20.0; 
	}
	return avg;
}


vector<double> Pssm::standard_deviation(vector<double>& avg)
{
	vector<double>std(12);
	for (int k = 0; k<12; k++)//calculate STD
	{
		for (int i = 0; i<20; i++)
		{
			std[k] = std[k] + pow((getij(i, k) - avg[k]), 2);

		}
		std[k] = sqrt(std[k] / (double)20);
	}
	return std;

}


vector<int> Pssm::count_the_amino_asid(msa_matrix msa,int i)
{
	GlobalFunction globalFunc;
	vector<int>tempArr(20);

	for (int j = 0; j<cluster_row_size; j++)
	{
		int index = globalFunc.AminoAcidToInt(msa.getijaa(j, i));
			tempArr[index]++;
	}
	return tempArr;
}


void Pssm::set_global_pssm_count(int n)
{
	cluster_row_size = cluster_row_size + 2 * n;
}


int Pssm::get_global_pssm_count()const
{
	return cluster_row_size;
}


void Pssm::initPssm(msa_matrix msa, Pssm global_pssm)
{
	cout << "pssm.init" << endl;
	for (int i = 0; i<12; i++)//2 loops to create the PSSM matrix with data on Amino Acid
	{
		vector<int> row_of_AA = count_the_amino_asid(msa,i);
	
		for (int k = 0; k<20; k++)
		{
			setij(k, i, (double)row_of_AA[k] / cluster_row_size);
			double temp = (global_pssm.getij(k, i) + (double)row_of_AA[k]) / global_pssm.get_global_pssm_count();
			global_pssm.setij(k, i, temp);

		}
	}//Finish Create the PSSM matrix
}


void Pssm::print(ofstream &stream)
{
	GlobalFunction g;
	for (int i = 0; i<12; i++)
	{
		stream << '|' << setw(10) << (i + 1) << "	";
	}
	stream << endl;
	
	
	for (int k = 0; k<20; k++)//Print to file the PSSM
	{
		
		for (int i = 0; i<12; i++)
		{
			stream << '|' << setw(10) << getij(k, i) << "	";
		}
		stream <<"----> "<<g.intToAminoAcid(k)<< endl;
	}
}


double Pssm::sum_of_row(int k)
{
	int i;
	double sum = 0;
	for (i = 0; i < 12; i++)
		sum +=getij(k, i);
	return sum;
}

void Pssm::propensity(ofstream& myfile)
{
	vector<double> percent(20);
	Matrix prop(12,20);
	double sum=0.0;
	GlobalFunction g;

	for(int k=0;k<20;k++)
	{
		sum = sum_of_row(k);
		for(int i=0 ; i<12 ; i++)
		{
			prop.setij(k,i, getij(k,i)/(sum/12));
		}
		sum=0.0;
	}
	myfile<<"************************** propensity for each Amino Acid  ***************************"<<endl;
	prop.print(myfile);

}
void Pssm::gas_propensity(ofstream& myfile)
{
	double totalGAS=0.0;
	vector<double>gas(12);
	GlobalFunction g;

	totalGAS += sum_of_row(1);
	totalGAS += sum_of_row(8);
	totalGAS += sum_of_row(16);

	for (int i = 0; i<12; i++)
	{

		gas[i] = getij(1, i);
		gas[i] += getij(8, i);
		gas[i] += getij(16, i);
		gas[i] = gas[i] / (totalGAS / 12);
	}
	g.print_array(myfile,gas);

}
