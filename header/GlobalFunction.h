// Last edit:
// 14/12   14:00

#ifndef GLOBAL_FUNCTION_H
#define GLOBAL_FUNCTION_H
#include <vector>
#include <array>

class Pssm;
using namespace std;

	class GlobalFunction
	{
	public:
		GlobalFunction();//Empty Constructor

		//get path to dir and return vector with all file name in the dir 
		void GetFilesNamesFromFolder(char* path, vector<string>& filename);

		//get path to dir and return the number of file in the dir
		int GetNumberOfFilesInFolder(char* path);

		//get the csv file number and  return the number of row in the file 
		int get_amount_lines_in_file(int i);

		//create txt file open it and return it open 
		ofstream creat_file_format_txt(string fileName);
		//vector<string> split_pdb_row_to_tokens(string readout);
		//ifstream& open_pdb_file(ifstream& file, string pdbName, vector<string>& filename);
		void OpenPdbFile(ifstream& readFile, string filename);
		void OpenFastaFile(ifstream& readFile);
		void ReadProteinsFromFasta();
		void ProteinsCounter();
		void SmallAminoAcidLists();

		vector<string> split(const string& s, const string& delim, bool keep_empty = true);
		void openClusterCSV(ifstream& Cluster, int i);
		string threeToOne(const string &s, ofstream& msaErrorFile);
		string oneToThree(const string &s);
		int AminoAcidToInt(const string &s);
		string intToAminoAcid(const int s);
		double Hydrophobic(const string &s);
		void print_array(ofstream& myFile, vector<double>& vec);
	};

#endif