// Last edit:
// 14/12   14:00

/*
"Small-residue Dynamic Functional Motifs in Membrane Proteins"
The program will generate an analysis reports on a specific set of clusters and polytopic proteins.

Members of code programmer:
Sagi Sulimani
Marina Katchko
Yogev Levi
Haim
*/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cmath> 
#include "header\GlobalFunction.h" 

#include "header\Msa.h"
#include "header\pssm.h"
#include "header\CartesianGeometry.h"
#include "header\Helanal.h"
#include "header\BaseCartesianPoint.h"

using namespace MSL;
using namespace std;

int main(void)
{
    GlobalFunction globalFunc;

    cout << "******************************************************************" << endl;
    cout << "**                                                              **" << endl;
    cout << "**   Welcome to the scientific and geometric analysis program   **" << endl;
    cout << "**   usege. The program will generate an analysis reports on    **" << endl;
    cout << "**   a specific set of clusters and polytopic proteins          **" << endl;
    cout << "**                                                              **" << endl;
    cout << "**  please make sure that all pdb file for the Protein          **" << endl;
    cout << "**  definition are in folder C:\\Sresidue\\polytopic              **" << endl;
    cout << "**                                                              **" << endl;
    cout << "******************************************************************" << endl;
    cout << endl;
    cout << "Hit 'Enter' to start please ..." << endl;
    getchar();

    // create a vector with the amount of files in the folder
    vector<string> filename(globalFunc.GetNumberOfFilesInFolder("C:\\Sresidue\\polytopic"));

    // initial the vector with the protein pdb file name 
    globalFunc.GetFilesNamesFromFolder("C:\\Sresidue\\polytopic", filename);

    // finds the proteins names for faste file and write it to 'Proteins Names' file  
    globalFunc.ReadProteinsFromFasta();

    // counts how many proteins there is in the 'Proteins Names' file 
    globalFunc.ProteinsCounter();

    // Make list of all the small amino acid (G,A,S,T,C) in the fasta file
    globalFunc.SmallAminoAcidLists();

    /*
    ofstream errorLog;
    errorLog = globalFunc.creat_file_format_txt("errorLog.txt");
    errorLog << "this is the error in the System"<<endl<<endl<<endl;

    ofstream msaErrorLog;
    msaErrorLog = globalFunc.creat_file_format_txt("msaErrorLog.txt");
    msaErrorLog << "this is the error in the MSA" << endl << endl << endl;

    bool msaFlag = true;
    bool pssmFlag = true;
    bool rmsdFlag = true;

    //loop for 16 clusters
    for (int i = 1; i <= 16; i++)
    {

    cout << "main loop Iteration " << i << " of 16" << endl;
    int  clusterleangth;
    //clusterleangth hold the number of row in csv file cluster(i)
    clusterleangth = globalFunc.get_amount_lines_in_file(i);
    errorLog << "clusterleangth = " << clusterleangth << " Iteration " << i << " of 16  " << endl;
    //create msa matrix msaArr[n][12] msaArr2[n][12]
    msa_matrix msaArr(clusterleangth);
    msa_matrix msaArr2(clusterleangth);
    //create the global pssm matrix
    Pssm global_pssm(0,0);
    vector<string> resSplit;
    //create the msa matrix
    if (msaFlag == true)
    {
    ofstream msafile;
    string fileName = "cluster" + to_string(i) + ".txt";
    //create the txt file to print the result
    msafile = globalFunc.creat_file_format_txt(fileName);

    //cluster_row_size = cluster_row_size + 2 * clusterleangth
    global_pssm.set_global_pssm_count(clusterleangth);


    ifstream cluster;
    //open the cluster csv file to read the rows
    globalFunc.openClusterCSV(cluster, i);
    string str, pdbName;
    //loop to analyze the csv file
    int row = 0;

    while (getline(cluster, str))
    {
    vector<string> tmData;
    //split the row in to a tokens
    tmData = globalFunc.split(str, ",");
    //pull the pdb name out
    pdbName = tmData[5].substr(0, 4);

    ifstream readFile;
    //Checks whether the file name is in the folder, then open it
    globalFunc.open_pdb_file(readFile, pdbName, filename);
    if (!readFile.is_open())
    {
    errorLog <<" Fails to open file : " << pdbName  << endl<<endl;
    continue;
    }

    string readout;
    int Column = 0;
    cout << "init msaArr and msaArr2 from file : " << pdbName << " start" << endl;
    //loop to read the pdb file and criate the msa matrix
    while (getline(readFile, readout))
    {


    resSplit=globalFunc.split_pdb_row_to_tokens(readout);

    msaArr.init(resSplit, tmData, "msaArr", row, errorLog);
    msaArr2.init(resSplit, tmData, "msaArr2", row, errorLog);

    }
    for (int t = 0; t < 12; t++)
    {
    if (msaArr.getijaa(row, t) == "" )
    {
    msaErrorLog << "Error in msaArr ----> pdb file name: " << pdbName << " subUnit: " << tmData[6] << " Amino acid number: " << tmData[7] <<" cluster number " << i << endl;
    }
    if (msaArr2.getijaa(row, t) == "")
    {
    msaErrorLog << "Error in msaArr2 ----> pdb file name: " << pdbName << " subUnit: " << tmData[8] << " Amino acid number: " << tmData[9] << " cluster number " << i << endl;
    }
    }
    cout << "init msaArr and msaArr2 from file : " << pdbName << " end  " << "row=" << row << endl;
    row++;
    readFile.close();
    }
    cluster.close();
    msafile << "--------------------------------------- msaArr------------------------------------------------------" << endl;
    msaArr.hydroZcoorAvg(msafile, clusterleangth);

    msafile << "--------------------------------------- msaArr2------------------------------------------------------" << endl;
    msaArr2.hydroZcoorAvg(msafile, clusterleangth);
    cout << "print msa" << endl;
    msafile << "                                        msaArr" << endl;
    msaArr.print(msafile);
    msafile << endl << endl << endl << endl << endl;
    msafile << "                                        msaArr2" << endl;
    msaArr2.print(msafile);
    msafile.close();
    }




    if (pssmFlag == true)
    {
    //create the txt file to print the result
    ofstream pssmFile;
    string fileName = "pssm" + to_string(i) + ".txt";
    cout << "pssm file name is : " << fileName << endl;
    pssmFile = globalFunc.creat_file_format_txt(fileName);

    cout << "enter into pssm Iteration  " << i << " of 16 " << endl;
    Pssm pssm1(clusterleangth,1);
    cout << " pssm1 matrix create Iteration " << i << " of 16  " << endl;
    Pssm pssm2(clusterleangth,2);
    cout << " pssm2 matrix create Iteration  " << i << " of 16  " << endl;
    pssm1.initPssm(msaArr, global_pssm);
    cout << " pssm1.init work  " << i << " of 16  " << endl;
    pssm2.initPssm(msaArr2, global_pssm);
    pssmFile << endl << endl;
    pssmFile <<  "**************************************************** PSSM matrix for Helix 1 *****************************************************" << endl;
    pssmFile << endl << endl;
    pssm1.print(pssmFile);
    pssmFile << endl << endl << endl << endl;
    pssmFile << "**************************************************** PSSM matrix for Helix 2 *****************************************************" << endl;
    pssmFile << endl << endl;
    pssm2.print(pssmFile);

    errorLog << "pssm1 Iteration " << i << " of 16  " << endl;
    pssm1.getSumOfColumns(errorLog);
    errorLog << "pssm2 Iteration " << i << " of 16  " << endl;
    pssm2.getSumOfColumns(errorLog);
    vector<double> avg1;
    vector<double> avg2;
    vector<double> std1;
    vector<double> std2;
    avg1 = pssm1.average();
    std1 = pssm1.standard_deviation(avg1);
    avg2 = pssm2.average();
    std2 = pssm2.standard_deviation(avg2);
    pssmFile << endl << endl << endl << endl;
    pssmFile << "****************************************************** average  for Helix 1 ********************************************************" << endl;
    globalFunc.print_array(pssmFile, avg1);
    pssmFile << endl << endl;
    pssmFile << "****************************************************** average  for Helix 2 ********************************************************" << endl;
    globalFunc.print_array(pssmFile, avg2);
    pssmFile << endl << endl;
    pssmFile << "************************************************* standard deviation  for Helix 1  *************************************************" << endl;
    globalFunc.print_array(pssmFile, std1);
    pssmFile << endl << endl;
    pssmFile << "************************************************* standard deviation  for Helix 2  *************************************************" << endl;
    globalFunc.print_array(pssmFile, std2);
    pssmFile << endl << endl;
    pssmFile << "*********************************************** Helix 1 propensity for GAS Amino Acid ********************************************** " << endl;
    pssm1.gas_propensity(pssmFile);
    pssmFile << endl << endl;
    pssmFile << "*********************************************** Helix 2 propensity for GAS Amino Acid  ********************************************* " << endl;
    pssm2.gas_propensity(pssmFile);
    pssmFile << endl << endl;
    pssmFile << "                      Helix 1 propensity  " << endl;
    pssm1.propensity(pssmFile);
    pssmFile << endl << endl;
    pssmFile << "                      Helix 2 propensity " << endl;
    pssm2.propensity(pssmFile);
    pssmFile << endl << endl;
    pssmFile.close();
    }
    if (rmsdFlag == true)
    {
    ifstream cluster;
    ofstream pdbFile;
    ofstream pdbfullCluster;
    //open the cluster csv file to read the rows
    globalFunc.openClusterCSV(cluster, i);
    string str, pdbName;
    BaseCartesianPoint  *cluster_arr = new BaseCartesianPoint[clusterleangth];
    string pdbCluster = "pdbCluster" + to_string(i) + "\\pdbCluster" + to_string(i) + ".txt";
    pdbfullCluster = globalFunc.creat_file_format_txt(pdbCluster);
    if (!pdbfullCluster.is_open())
    {
    errorLog << " Fails to open file : " << pdbCluster << endl << endl;
    }
    for (int x = 0; x < clusterleangth; x++)
    {

    cluster_arr[x].setMsaPoints(msaArr, msaArr2, x);
    cluster_arr[x].setCenter();
    cluster_arr[x].calculatAvgZ();
    cluster_arr[x].calculatRmsd();
    getline(cluster, str);
    vector<string> tmData;
    //split the row in to a tokens
    tmData = globalFunc.split(str, ",");
    //pull the pdb name out
    string pdbFileName = "pdbCluster" + to_string(i) + "\\" + tmData[5] + ".txt";

    pdbFile = globalFunc.creat_file_format_txt(pdbFileName);
    if (!pdbFile.is_open())
    {
    errorLog << " Fails to open file : " << pdbFileName << endl << endl;
    }
    cluster_arr[x].printPDB(pdbFile, tmData, pdbfullCluster);
    tmData.clear();
    pdbFile.close();
    }
    BaseCartesianPoint::print(i, cluster_arr, clusterleangth); //static function

    pdbfullCluster.close();
    }

    }
    // End of the clusters loop
    msaErrorLog.close();
    errorLog.close();
    */

    cout << endl;
    cout << "******************************************************************" << endl;
    cout << "**   The analysis completed.                                    **" << endl;
    cout << "**   please check the errorLog file and the msaErrorLog file    **" << endl;
    cout << "**   to make sure that all analysis complete successfully.      **" << endl;
    cout << "**                                                              **" << endl;
    cout << "**   All the output files are in the program folder.            **" << endl;
    cout << "******************************************************************" << endl;
    getchar();
    return 0;
}
