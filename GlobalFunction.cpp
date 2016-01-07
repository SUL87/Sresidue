// 6.1.16

#include "header\dirent.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "header\Msa.h"
#include "header\GlobalFunction.h" 
#include "header\pssm.h"

using namespace std;

GlobalFunction::GlobalFunction()
{
}

// Get names of files in the folder of char *path
void GlobalFunction::GetFilesNamesFromFolder(char* path, vector<string>& filename)
{
    struct dirent *dp;
    DIR *fd;
    ofstream pdbList("PDB List.txt");
    int i = 0;

    if ((fd = opendir(path)) == NULL)
    {
	   fprintf(stderr, " can't open %s\n", path);
	   cout << "************** Error in function: GetFilesNamesFromFolder ************************" << endl;
	   getchar();
	   return;
    }

    if (!pdbList.is_open())
    {
	   cout << "*********************** Unable to open 'pdb list.txt' file ***********************" << endl;
	   getchar();
    }

    while ((dp = readdir(fd)) != NULL)
    {
	   if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
		  continue;
	   else
	   {
		  pdbList << dp->d_name;
		  pdbList.seekp(-4, ios::cur);	   // go -4 to avoid ".pdb"
		  pdbList << "\n";
		  filename[i] = dp->d_name;
		  i++;
	   }
    }
    closedir(fd);
}

int GlobalFunction::GetNumberOfFilesInFolder(char* path)
{
    struct dirent *dp;
    DIR *fd;
    int count = 0;
    if ((fd = opendir(path)) == NULL)
    {
	   fprintf(stderr, " can't open %s\n", path);
	   cout << "**************error in function: GetNumberOfFilesInFolder ************************" << endl;
	   getchar();
	   return 0;
    }
    while ((dp = readdir(fd)) != NULL)
    {
	   if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
		  continue;    /* skip self and parent */
	   else
		  count++;
    }
    closedir(fd);
    return count;
}

// Just open the pdb file. If cannot open, print error message.
void GlobalFunction::OpenPdbFile(ifstream& readFile, string filename)
{
    ostringstream name;

    name << "C:\\Sresidue\\polytopic\\" + filename;
    string newname = name.str().c_str();
    readFile.open(newname, ios::in);
    if (!readFile.is_open())
    {
	   cout << "Unable to open file: " << "C:\\Sresidue\\polytopic\\" << filename << endl;
	   getchar();
    }
    return;
    getchar();
}

// Just open the fasta file. If cannot open, print error message.
void GlobalFunction::OpenFastaFile(ifstream &readFile)
{
    ostringstream name;

    name << "C:\\fasta\\fasta.fas";
    string newname = name.str().c_str();
    readFile.open(newname, ios::in);
    if (!readFile.is_open())
    {
	   cout << "Unable to open file: " << "C:\\fasta\\fasta.fas" << endl;
	   getchar();
    }
    return;
    getchar();
}

// Create a text file of list of all the proteins in a fasta file
void GlobalFunction::ReadProteinsFromFasta()
{
    ifstream readFile;				// readFile is a pointer to the source fasta file
    ofstream fout;	    				// fout is a pointer to the Proteins Names file (output file)
    char ch, string[5];
    int i = 1;
    bool flag = true;

    printf("Put the fasta file in C:\\fasta\\ and call it fasta.fas\nHit 'Enter' to continue please.");
    getchar();
    OpenFastaFile(readFile);			// open the fasta file with the pointer readFile

    fout.open("Proteins Names.txt", ios::out); // open the output file with the pointer fout

    while (readFile.get(ch))			// read char till the end of file
    {
	   if (ch == '>')				// the name of the protein comes after '>'
		  if (flag == true)			// there is another '>' in the paragraph, so we use a flag to avoid it
		  {
			 //readFile.get(ch);		// avoid the digit before the protein name
			 //readFile.get(line, 4);	// copy the protein name: 3 chars for the name + 1 char for '\0'
			 readFile.get(string, 5);	// copy the protein name: 4 chars for the name + 1 char for '\0'
			 fout << string << endl;	// print the protein name in the output file ("Proteins Names")
			 flag = false;
		  }
		  else
			 flag = true;
    }
    readFile.close();
}

// Count the number of proteins in the file 'Proteins Names.txt' which has the names from the fasta file
void GlobalFunction::ProteinsCounter()
{
    ifstream fin("Proteins Names.txt");
    vector<string> vector;
    string str;

    while (true)
    {
	   getline(fin, str);
	   if (fin.eof())
		  break;
	   vector.push_back(str);
    }
    fin.close();

    sort(vector.begin(), vector.end());

    auto it = unique(begin(vector), end(vector));
    vector.erase(it, vector.end());

    ofstream fout("Proteins Names.txt");
    copy(vector.begin(), vector.end(), ostream_iterator<string>(fout, "\n"));

    printf("\nAmount of proteins in the fasta file: %d\n", vector.size());
    fout.close();
}

// Make list of all the small amino acid (G,A,S,T,C) in the fasta file
void GlobalFunction::SmallAminoAcidLists()
{
    ifstream readFile;					   // readFile is a pointer to the source fasta file
    ofstream fout2, fout3, fout4;	    		   // fout2, fout3, fout4 are pointers to the Small Amino Acid Lists file (output file)
    vector<string> vector2, vector3, vector4;
    string str, proteinName, line;
    char subunit, ch1, ch2, ch3, ch4;
    int cnt = 0, length = 0, position = 0;
    OpenFastaFile(readFile);				    // open the fasta file with the pointer readFile
    fout2.open("2 Small AA List.txt", ios::out); // open the output file with the pointer fout2
    fout3.open("3 Small AA List.txt", ios::out); // open the output file with the pointer fout3
    fout4.open("4 Small AA List.txt", ios::out); // open the output file with the pointer fout4

    while (readFile.get(ch1))				   // read the file till the end of file
    {
	   str.clear();
	   if (ch1 == '>')					   // skip the line of data which start with '>'
	   {
		  position = 0;
		  proteinName.clear();
		  for (int i = 0; i < 4; i++)
		  {
			 readFile.get(ch1);
			 proteinName += ch1;
		  }
		  readFile.get(subunit);
		  getline(readFile, line);		   // skip the rest of the line
		  continue;					   // skip the rest of the loop
	   }
	   if (ch1 != '\n')
		  position++;
	   if (ch1 == 'G' || ch1 == 'A' || ch1 == 'S' || ch1 == 'T' || ch1 == 'C')
		  if (readFile.get(ch2))
		  {
			 if (ch2 != '\n')
				position++;
			 if (ch2 == 'G' || ch2 == 'A' || ch2 == 'S' || ch2 == 'T' || ch2 == 'C')
			 {
				length = 2;
				str += ch1; str += ch2;
				if (readFile.get(ch3))
				{
				    if (ch3 != '\n')
					   position++;
				    if (ch3 == 'G' || ch3 == 'A' || ch3 == 'S' || ch3 == 'T' || ch3 == 'C')
				    {
					   length = 3;
					   str += ch3;
					   if (readFile.get(ch4))
					   {
						  if (ch4 != '\n')
							 position++;
						  if (ch4 == 'G' || ch4 == 'A' || ch4 == 'S' || ch4 == 'T' || ch4 == 'C')
						  {
							 length = 4;
							 str += ch4;
						  }
					   }
					   else break;	    // end of file
				    }
				}
				else break;	    // end of file
			 }
		  }
	   else break;	    // end of file

	   if (length == 2)
	   {
		  fout2 << str << '\t';
		  fout2 << proteinName << '\t';
		  fout2 << subunit << '\t';
		  fout2 << position-2 << '\t';
		  fout2 << '\n';
	   }
	   else if (length == 3)
	   {
		  fout3 << str << '\t';
		  fout3 << proteinName << '\t';
		  fout3 << subunit << '\t';
		  fout3 << position-3 << '\t';
		  fout3 << '\n';
	   }
	   else if (length == 4)
	   {
		  fout4 << str << '\t';
		  fout4 << proteinName << '\t';
		  fout4 << subunit << '\t';
		  fout4 << position-3 << '\t';    // sequences of length of 4 aa don't increase the position, that's why the -3
		  fout4 << '\n';
	   }
	   length = 0;				 // initialize length
    } // end while

    readFile.close();
}

void GlobalFunction::FindWodakInFasta()
{
    ifstream readFile;
    ofstream outFile("Fasta to Wodak.txt");
    string string, num;
    vector<int> result;

    OpenFastaFile(readFile);
    if (readFile.is_open())
    {
	   while (getline(readFile, string))
	   {
		  if (string[0] == '>')
			 {
				outFile << "\n\n" << string << "\n";
				continue;
			 }
		  if (string[0] == 10)	 // if the next char is Line Feed (Enter) - skip it
			 outFile << "\n\n";
		  for (int i = 0; i < string.length(); i++)
		  {
			 num = WodakTableValues(string[i]);
			 if (num != "0")
			 {
				outFile << num << ",";
				//result.push_back(N);
			 }
		  }
		  //copy(result.begin(), result.end(), ostream_iterator<int>(outFile, "\n"));
	   }
	   readFile.close();
	   outFile.close();
    }
    else cout << "Unable to open fasta file.";
}
string GlobalFunction::WodakTableValues(char ch)
{
    switch (ch)
    {
    case 'A':	return "316"; //Ala
	   break;
    case 'R':	return "9"; //Arg
	   break;
    case 'N':	return "27"; //Asn
	   break;
    case 'D':	return "30"; //Aspr
	   break;
    case 'C':	return "34";  //Cysh
	   break;
    case 'Q':	return "12"; //Gln
	   break;
    case 'E':	return "7"; //Glu
	   break;
    case 'G':	return "239"; //Gly
	   break;
    case 'H':	return "22"; //His
	   break;
    case 'I':	return "212"; //Ile
	   break;
    case 'L':	return "226"; //Leu
	   break;
    case 'K':	return "4"; //Lys
	   break;
    case 'M':	return "56"; //Met
	   break;
    case 'F':	return "84"; //Phe
	   break;
    case 'P':	return "25"; //Pro
	   break;
    case 'S':	return "109"; //Ser
	   break;
    case 'T':	return "66";  //Thr
	   break;
    case 'W':	return "27"; //Trp
	   break;
    case 'Y':	return "27";  //Tyr
	   break;
    case 'V':	return "308";  //Val
	   break;
    default:	return "0";
	   break;
    }
}

void GlobalFunction::print_array(ofstream& myFile, vector<double>& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
	   myFile << vec[i] << " , ";
    }
    myFile << endl;
    myFile << endl;
    myFile << endl;
}

ofstream GlobalFunction::creat_file_format_txt(string fileName)
{
    ofstream myfile(fileName);
    if (!myfile.is_open())
    {
	   cout << "*********************** Unable to open file ***********************" << fileName << endl;
	   getchar();
    }

    return myfile;
}

vector<string> GlobalFunction::split(const string& s, const string& delim, bool keep_empty) {
    keep_empty = true;
    vector<string> result;
    if (delim.empty())
    {
	   result.push_back(s);
	   return result;
    }
    string::const_iterator substart = s.begin(), subend;
    while (true)
    {
	   subend = search(substart, s.end(), delim.begin(), delim.end());
	   string temp(substart, subend);
	   if (keep_empty || !temp.empty())
	   {
		  result.push_back(temp);
	   }
	   if (subend == s.end())
	   {
		  break;
	   }
	   substart = subend + delim.size();
    }
    return result;
}

string GlobalFunction::threeToOne(const string &s, ofstream& msaErrorFile)
{
    string result = s;

    if (result.size() != 3)
    {
	   result = result.substr(result.size() - 3);
    }
    if (result == "ALA")
	   return "A";
    if (result == "ARG")
	   return "R";
    if (result == "ASN")
	   return "N";
    if (result == "ASP")
	   return "D";
    if (result == "CYS")
	   return "C";
    if (result == "GLN")
	   return "Q";
    if (result == "GLU")
	   return "E";
    if (result == "GLY")
	   return "G";
    if (result == "HIS")
	   return "H";
    if (result == "ILE")
	   return "I";
    if (result == "LEU")
	   return "L";
    if (result == "LYS")
	   return "K";
    if (result == "MET")
	   return "M";
    if (result == "PHE")
	   return "F";
    if (result == "PRO")
	   return "P";
    if (result == "SER")
	   return "S";
    if (result == "THR")
	   return "T";
    if (result == "TRP")
	   return "W";
    if (result == "TYR")
	   return "Y";
    if (result == "VAL")
	   return "V";
    cout << " ******************error in function : threeToOne  ****************" << endl;
    cout << "Error: there is no amino asid with the name " << s << endl;
    getchar();
    return "error";
}

string GlobalFunction::oneToThree(const string &s)
{

    if (s == "A")
	   return "ALA";
    if (s == "R")
	   return "ARG";
    if (s == "N")
	   return "ASN";
    if (s == "D")
	   return "ASP";
    if (s == "C")
	   return "CYS";
    if (s == "Q")
	   return "GLN";
    if (s == "E")
	   return "GLU";
    if (s == "G")
	   return "GLY";
    if (s == "H")
	   return "HIS";
    if (s == "I")
	   return "ILE";
    if (s == "L")
	   return "LEU";
    if (s == "K")
	   return "LYS";
    if (s == "M")
	   return "MET";
    if (s == "F")
	   return "PHE";
    if (s == "P")
	   return "PRO";
    if (s == "S")
	   return "SER";
    if (s == "T")
	   return "THR";
    if (s == "W")
	   return "TRP";
    if (s == "Y")
	   return "TYR";
    if (s == "V")
	   return "VAL";

    cout << " ******************error in function : oneToThree  ****************" << endl;
    cout << "Error: there is no amino asid with the name " << s << endl;
    getchar();
    return "Error";
}

int GlobalFunction::AminoAcidToInt(const string &s)
{

    if (s == "A")
	   return  0;
    if (s == "R")
	   return  1;
    if (s == "N")
	   return 2;
    if (s == "D")
	   return 3;
    if (s == "C")
	   return 4;
    if (s == "Q")
	   return 5;
    if (s == "E")
	   return 6;
    if (s == "G")
	   return 7;
    if (s == "H")
	   return 8;
    if (s == "I")
	   return 9;
    if (s == "L")
	   return 10;
    if (s == "K")
	   return 11;
    if (s == "M")
	   return 12;
    if (s == "F")
	   return 13;
    if (s == "P")
	   return 14;
    if (s == "S")
	   return 15;
    if (s == "T")
	   return 16;
    if (s == "W")
	   return 17;
    if (s == "Y")
	   return 18;
    if (s == "V")
	   return 19;
    cout << " ******************error in function : AminoAcidToInt  ****************" << endl;
    cout << "Error: there is no amino asid with the name " << s << endl;
    getchar();
    return -1;
}

string GlobalFunction::intToAminoAcid(const int s)
{

    if (s == 0)
	   return "ALA";
    if (s == 1)
	   return "ARG";
    if (s == 2)
	   return "ASN";
    if (s == 3)
	   return "ASP";
    if (s == 4)
	   return "CYS";
    if (s == 5)
	   return "GLN";
    if (s == 6)
	   return "GLU";
    if (s == 7)
	   return "GLY";
    if (s == 8)
	   return "HIS";
    if (s == 9)
	   return "ILE";
    if (s == 10)
	   return "LEU";
    if (s == 11)
	   return "LYS";
    if (s == 12)
	   return "MET";
    if (s == 13)
	   return "PHE";
    if (s == 14)
	   return "PRO";
    if (s == 15)
	   return "SER";
    if (s == 16)
	   return "THR";
    if (s == 17)
	   return "TRP";
    if (s == 18)
	   return "TYR";
    if (s == 19)
	   return "VAL";

    cout << " ******************error in function : intToAminoAcid  ****************" << endl;
    cout << "Error: there is no amino asid with the name " << s << endl;
    getchar();
    return "Error";
}

double GlobalFunction::Hydrophobic(const string &s)
{

    if (s == "A")
	   return 0.62;
    if (s == "R")
	   return -2.53;
    if (s == "N")
	   return -0.78;
    if (s == "D")
	   return -0.90;
    if (s == "C")
	   return 0.29;
    if (s == "Q")
	   return -0.85;
    if (s == "E")
	   return -0.74;
    if (s == "G")
	   return 0.48;
    if (s == "H")
	   return -0.40;
    if (s == "I")
	   return 1.38;
    if (s == "L")
	   return 1.06;
    if (s == "K")
	   return -1.5;
    if (s == "M")
	   return 0.64;
    if (s == "F")
	   return 1.19;
    if (s == "P")
	   return 0.12;
    if (s == "S")
	   return -0.18;
    if (s == "T")
	   return -0.05;
    if (s == "W")
	   return 0.81;
    if (s == "Y")
	   return 0.26;
    if (s == "V")
	   return 1.08;
    cout << " ******************error in function : Hydrophobic  ****************" << endl;
    cout << "Error: there is no amino asid with the name " << s << endl;
    getchar();
    return -1;
}

/*
void GlobalFunction::openClusterCSV(ifstream& clusters, int i)
{
ostringstream name;
name<<"C:\\Sresidue\\cluster\\cluster"<<i<<".csv";
string newname = name.str().c_str();
clusters.open(newname, ios::in);
if (!clusters.is_open())
{
cout << "Unable to open file: " << "C:\\Sresidue\\cluster\\cluster" << i << ".csv" << endl;
getchar();
}

}
*/
/*
// Just open the pdb file. If cannot open, print error message.
ifstream& GlobalFunction::open_pdb_file(ifstream& readFile, string pdbName, vector<string>& filename)
{
for (size_t j = 0; j < filename.size(); j++)	  //looping the files
{
if (filename[j] == pdbName + ".pdb")	  //check if file exist - if true: Create files only with Atoms to parse it
{
//cout << filename[j] << " == " << pdbName << ".pdb" << endl;
ostringstream name;
name << "C:\\Sresidue\\polytopic\\" + pdbName + ".pdb";
string newname = name.str().c_str();
readFile.open(newname, ios::in);
if (!readFile.is_open())
{
cout << "Unable to open file: " << "C:\\Sresidue\\polytopic\\" << pdbName << ".pdb" << endl;
getchar();
}
return readFile;
}
}
cout << "cannot find pdb file : "<<pdbName << " in the folder" << endl;
getchar();
}
*/
/*
vector<string> GlobalFunction::split_pdb_row_to_tokens(string readout)
{
istringstream iss(readout);
vector<string>resSplit{ istream_iterator<string>(iss), istream_iterator<string>() };
if (resSplit.empty())
{
cout << "**** Failed to split in to tokens. resSplit is empty. ****" << endl;
getchar();
}

if (resSplit[0] == "ATOM")
{
if (resSplit[4].size() != 1)
{
string subUnit, aaNum;
subUnit = resSplit[4].substr(0, 1);
aaNum = resSplit[4].substr(1);
resSplit.insert(resSplit.begin()+4, subUnit);
resSplit[5]= aaNum;
}
}
return resSplit;
}
*/
/*
int GlobalFunction::get_amount_lines_in_file(int i)
{
ifstream clusters;
string str;
int count=0;
openClusterCSV(clusters,i);

while(getline(clusters,str))
count++;

clusters.clear();
clusters.seekg(0, ios::beg);
clusters.close();

return count;

}
*/
