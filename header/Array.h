#ifndef Array_H
#define Array_H

#include <iostream>

using namespace std;

typedef class Array
{
	private:

	int row;
	double value;
	
public:
	
	Array();//Empty Constructor
	Array(int n);//Constructor
	Array(const Array &src);//Copy constructor
	~Array();//Destructor
	Array& operator=(const Array &rhs);//Function which overload '=' .
   void seti( int i, double num);//Function which sets values in the matrix.
	double geti(int i)const;
	void print();
	void free_Array();
	double* getArray();
	
	    std::ostream& save(std::ostream& out) const {
			return out << double_matrix;
    }
		
};



 ostream& operator <<(ostream &stream,const Array& s);//Function which overload '<<'.





#endif