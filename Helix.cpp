#include<iostream>
#include<string>
#include "header\Helix.h" 


using namespace std;

Helix::Helix()//Empty constructor whitch put null in msa_arr.
{
size=12;
helix_arr = new helixStruct[size];
}


Helix::Helix(int n){//Implementation of a costructor which recieve parmater n as the row of msa_arr.
	size=n;
	helix_arr = new helixStruct[size];
	for (int i=0; i<size; i++)
	{
			helix_arr[i].atom="";
			helix_arr[i].coorX=0.0;
			helix_arr[i].coorY=0.0;
			helix_arr[i].coorZ=0.0;
	}	
}

Helix::Helix(const Helix &src)//Implementation of  a  copy Constructor.
{
	size=src.size;
	helix_arr=new helixStruct[src.size];
	for (int i=0; i<src.size; i++)
		memcpy(&helix_arr[i],&src.helix_arr[i],sizeof(helixStruct));
	
}

void Helix::freeHelix(){//Function which free dynamic allocated memory space.
		delete 	[] helix_arr;
}

Helix::~Helix()//Destructor which call free_matrix to deallocate dynamic memory space.
{
	freeHelix();
}


void Helix::setHelix(string atom ,double x,double y, double z, int i)//Function which sets values in the matrix.
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
			else
			{					
				helix_arr[i].atom=atom;
				helix_arr[i].coorX=x;
				helix_arr[i].coorY=y;
				helix_arr[i].coorZ=z;
			}
}

Helix& Helix::operator=(const Helix &rhs)//Function which overload '=' .
{
	
	if(this==&rhs)
		return *this;
	
	if(helix_arr!=NULL)
		freeHelix();
	
	size=rhs.size;
	helix_arr=new helixStruct[rhs.size];
	for (int i=0; i<rhs.size; i++)
		memcpy(&helix_arr[i],&rhs.helix_arr[i],sizeof(helixStruct));
	return *this;
}	


ostream &operator <<(ostream &stream, Helix src)//Function which overload '<<'.
{
	src.printHelix();
	return stream;
}

void Helix:: printHelix()//Function which print the matrix values.
{
	cout<<"****************** Helix ******************************";
	for(int i=0;i<size;i++)
	{
		cout<<getHelixAtom(i)<<" "<<getHelixCoorX(i)<<" "<<getHelixCoorY(i)<<" "<<getHelixCoorZ(i)<<"   ,   ";
		cout<<endl;
	}
}

helixStruct Helix:: getHleixStruct(int i)const//Get the value of the place ij in the matrix.
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
		else
			return helix_arr[i];
	
}

string Helix:: getHelixAtom(int i)const//Get the value of the place ij in the matrix.
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
			else				
				return helix_arr[i].atom;
	return NULL;
}

double Helix::getHelixCoorX(int i)const
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
			else				
				return helix_arr[i].coorX;
	return NULL;
}

double Helix::getHelixCoorY(int i)const
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
			else				
				return helix_arr[i].coorY;
	return NULL;
}

double Helix::getHelixCoorZ(int i)const
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else 
		if(i<0||i>size )
			cout <<"the index "<<i<<" mast be between 0-"<< size<<endl ;
			else				
				return helix_arr[i].coorZ;
	return NULL;
}

helixP Helix:: getHelix_arr()const//Get the value of the place ij in the matrix.
{
	if(helix_arr==NULL)
		cout <<"helix_arr == NULL" <<endl;
	else				
		return helix_arr;
	return NULL;
}

