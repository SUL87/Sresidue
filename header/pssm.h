
#ifndef PSSM_H
#define PSSM_H

#include <string>
#include <iostream>
#include "Msa.h"
#include "GlobalFunction.h"
#include "Matrix.h"
using namespace std;

	class Pssm : public Matrix
	{
		private:
			int cluster_row_size;
			int flag;
		public:

			Pssm(int n, int type);
			int get_global_pssm_count()const;
			void set_global_pssm_count(int n);
			void initPssm(msa_matrix msa, Pssm global_pssm);
			vector<int> count_the_amino_asid(msa_matrix msa, int i);
			void print(ofstream &stream);
			vector<double> average();
			vector<double> standard_deviation(vector<double>& avg);
			void propensity(ofstream& myfile);
			void gas_propensity(ofstream& myfile);
			double sum_of_row(int k);
	};



#endif
