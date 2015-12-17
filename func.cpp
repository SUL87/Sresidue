#include "dirent.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <sstream>
#include <iomanip>
#include "Msa.h"
#include <cstring>
#include <cmath> 
#include "func.h"
#include "GlobalFunction.h" 






using namespace std;


int functions::outputClusters()
{
	//GlobalFunction globalFunc;
	//vector<string> filename(globalFunc.get_number_of_files_in_folder("c:\\project\\atoms"));/*Get the number of files in the folder in order to set the Vector array to the Filename*/
	//
	//globalFunc.get_file_name_from_dir("c:\\project\\atoms",filename);/*Get the file name of the files an store them in the Vector array called "filename"*/
	int sum=0;
	int tempArr[21]={0};
	double pssm[21][12]={0.0};
	double pssm2[21][12]={0.0};
	double global_pssm[21][12]={0.0};
	int global_pssm_count=0;
	double stdVal[12]={0.0};
	double avgVal[12]={0.0};
	double stdVal2[12]={0.0};
	double avgVal2[12]={0.0};
	int tempArr2[21]={0};
	struct dirent *ent;
	int i=0,count=0,countForMatrix=0,Container_Count=0,countForMatrix2=0,Container_Count2=0;
	string file_contens;	
	ifstream clusters;
	vector<string> tmData,resSplit,pdbName;
	string str;
	int leangth[17] = {0};
	ostringstream name;

	//for (int i = 1; i<=16 ;i++)
	//{
	//	name<<"C:\\project\\cluster\\cluster"<<i<<".csv";
	//	string newname = name.str().c_str();
	//	clusters.open(newname, ios::in);
	//	while(getline(clusters,str))
	//		count++;
	//	leangth[i] = count;
	//	clusters.clear();
	//	clusters.seekg(0, ios::beg);
	//}
	//for(int z=1;z<=16;z++)
	//{
	//msa_matrix msaArr(leangth[z]);
	//msa_matrix msaArr2(leangth[z]);
	//ostringstream oss;
	//oss<<"cluster"<<z<<".txt";
	//ofstream myfile(oss);
	//string clusterName = oss.str().c_str();
	//myfile.open (clusterName);
//	while(getline(clusters,str)) /////////////////// stoped here
	//{
		//tmData = globalFunc.split(str, ",");
//		pdbName = globalFunc.split(tmData[5],"_");
//		for(int j=0;j<filename.size();j++)//looping threw the  files
//		{
			
//			if(filename[j]==pdbName[0]+".pdb")//check if file exist - if true: Create files only with Atoms to parse it 
//			{
//				ifstream readFile(("c:\\project\\atoms\\"  + pdbName[0]+".pdb").c_str());
//				ofstream outFile;
//				string readout;
//				while(getline(readFile,readout))//Parse to create the new File
//				{
//					string sentence = readout.c_str();
//					istringstream iss(sentence);
//					vector<string> atomFromFile;
//					copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(atomFromFile));
//				    if(!atomFromFile.empty())
//					if(atomFromFile[0]=="ATOM")
//					{
//					    	vector<string> str = globalFunc.split(tmData[7],"-");			
//
//						if(atomFromFile[4] == tmData[6]);
//						{ 
//							if(strcmp(atomFromFile[2].c_str(),"CA") == 0)
//							{
//								string st = atomFromFile[5];
//								if(atoi(st.c_str())>=atoi(str[0].c_str()) && atoi(st.c_str())<=atoi(str[1].c_str()) && Container_Count<12 )
//								{	
//									msaArr.setij(globalFunc.threeToOne((atomFromFile[3]).c_str()),atof(atomFromFile[8].c_str()),globalFunc.Hydrophobic(globalFunc.threeToOne((atomFromFile[3]).c_str())),countForMatrix,Container_Count);
//									Container_Count++;
//								}
//								
//						    }
//						}
//					}
//					
//					if(!atomFromFile.empty())
//					if(atomFromFile[0]=="ATOM")
//					{
//					    	vector<string> str = globalFunc.split(tmData[9],"-");			
//
//						if(atomFromFile[4] == tmData[8]);
//						{ 
//							if(strcmp(atomFromFile[2].c_str(),"CA") == 0)
//							{
//								string st = atomFromFile[5];
//								if(atoi(st.c_str())>=atoi(str[0].c_str()) && atoi(st.c_str())<=atoi(str[1].c_str()) && Container_Count2<12 )
//								{
//									msaArr2.setij(globalFunc.threeToOne((atomFromFile[3]).c_str()),atof(atomFromFile[8].c_str()),globalFunc.Hydrophobic(globalFunc.threeToOne((atomFromFile[3]).c_str())),countForMatrix2,Container_Count2);	
//									Container_Count2++;
//								}
//								
//						    }
//						}
//					}		
//				}
//				Container_Count2=0;
//				Container_Count=0;		
//			    readFile.close();
//			}//end if
//		}//end for
//				countForMatrix++;
//				countForMatrix2++;
//	}//end while
//	msaArr.print();
//	myfile<<"This is The PSSM matrix for Helix 1"<<endl;
//	global_pssm_count=global_pssm_count+count*2;
//	for(int i=0;i<12;i++)//2 loops to create the PSSM matrix with data on Amino Acid
//	{
//		for(int j=0;j<count;j++)
//		{
//			int index = globalFunc.AminoAcidToInt(msaArr.getijaa(j,i));
//			tempArr[index]++;
		
//		}
//		for(int k=0;k<21;k++)
//		{
//			pssm[k][i]=(double)tempArr[k] / (double)count;
//			global_pssm[k][i] = global_pssm[k][i]  + (double)tempArr[k] ;
//		}
//		for(int k=0;k<21;k++)
//			tempArr[k]=0;
//	}//Finish Create the PSSM matrix
//	for(int k=1;k<21;k++)//Print to file the PSSM
//		{
//			for(int i=0;i<12;i++)
//			{
//				myfile<<std::fixed << std::setprecision(3) <<pssm[k][i]<<"	";
//			}
//			myfile<<endl;
//		}
//	
//	myfile<<endl;
//	myfile<<endl;
//
//	sum=0;
//	myfile<<"This is The PSSM matrix for Helix 2"<<endl;
	
//	for(int i=0;i<12;i++)//2 loops to create the PSSM matrix with data on Amino Acid
//	{
//		for(int j=0;j<count;j++)
//		{
//			int index = globalFunc.AminoAcidToInt(msaArr2.getijaa(j,i));
//			tempArr2[index]++;
//		}
//				
//		for(int k=0;k<21;k++)
//		{
//			global_pssm[k][i] = global_pssm[k][i]  + (double)tempArr2[k] ;
//			pssm2[k][i]=(double)tempArr2[k] / (double)count;
//		}
//		for(int k=0;k<21;k++)
//		tempArr2[k]=0;
//	}//Finish Create the PSSM matrix
//	myfile<<endl;
//	myfile<<endl;
	//	for(int k=1;k<21;k++)//Print to file the PSSM
	//{
	//	for(int i=0;i<12;i++)
	//	{
	//		myfile<<std::fixed << std::setprecision(3) <<pssm2[k][i]<<"	";
	//	}
	//	myfile<<endl;
	//}
	
	//for(int k=0;k<12;k++)//calculate the Mean
	//{
	//	for(int i=1;i<21;i++)
	//	{
	//		avgVal[k] = avgVal[k] + pssm[i][k];
	//	}
	//	avgVal[k]= avgVal[k]/(double)21 ;
	//}
	//for(int k=0;k<12;k++)//calculate STD
	//{
	//	for(int i=1;i<21;i++)
	//	{
	//		stdVal[k] = stdVal[k] + (pssm[i][k] - avgVal[k])*(pssm[i][k] - avgVal[k]) ;
	//	}
	//	stdVal[k]=sqrt(stdVal[k]/(double)21);
	//	
	//}
		
	//myfile<<endl;
	//myfile<<endl;
	//myfile<<"The Mean is  for Helix 1"<<endl;
	//for(int i=0;i<12;i++)//print to file the Mean
	//{
	//	myfile<<avgVal[i]<<"	";
	//}
	//myfile<<endl;
	//myfile<<endl;
	//myfile<<"The Std is  for Helix1"<<endl;
	//for(int i=0;i<12;i++)//print to file the STD
	//{
	//	myfile<<stdVal[i]<<"	";
	//}
	//propensity(pssm,myfile);
	msaArr.hydroZcoorAvg(msaArr.getmsa_arr(),myfile,count);
	myfile<<endl;
	
	for(int k=0;k<12;k++)//calculate the Mean
	{
		for(int i=1;i<21;i++)
		{
			avgVal2[k] = avgVal2[k] + pssm2[i][k];
		}
		avgVal2[k]= avgVal2[k]/(double)21 ;
	}
	for(int k=0;k<12;k++)//calculate STD
	{
		for(int i=1;i<21;i++)
		{
			stdVal2[k] = stdVal2[k] + (pssm2[i][k] - avgVal2[k])*(pssm2[i][k] - avgVal2[k]) ;
		}
		stdVal2[k]= sqrt(stdVal2[k]/(double)21);
	}
		
	myfile<<endl;
	myfile<<endl;
	myfile<<"The Mean is  for Helix 2"<<endl;
	for(int i=0;i<12;i++)//print to file the Mean
	{
		myfile<<avgVal2[i]<<"	";
	}
	myfile<<endl;
	myfile<<endl;
	myfile<<"The Std is  for Helix 2"<<endl;
	for(int i=0;i<12;i++)//print to file the STD
	{
		myfile<<stdVal2[i]<<"	";
	}
	propensity(pssm2,myfile);
	msaArr2.hydroZcoorAvg(msaArr2.getmsa_arr(),myfile,count);
	myfile.close();
    sum=0;
	countForMatrix=0;
	Container_Count=0;
	countForMatrix2=0;
	Container_Count2=0;
	for(int i=0;i<12;i++)
		for(int j=0;j<21;j++){
			tempArr[j]=0;
			pssm[j][i]=0;
			tempArr2[j]=0;
			pssm2[j][i]=0;
		}
		count=0;
	}
	
	for(int i=0;i<21 ;i++){
		for(int j=0;j<12 ; j++)
			cout<<std::fixed << std::setprecision(3) <<global_pssm[i][j]/(double)global_pssm_count<<"	";
		cout<<""<<endl;
	}
	return 0;
}
