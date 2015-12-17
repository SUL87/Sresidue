

#include "header\CartesianPoint.h"
#include "header\CartesianGeometry.h"


using namespace std;
using namespace MSL;


CartesianPoint::CartesianPoint()
	: x(0.0), y(0.0), z(0.0)
{
}

CartesianPoint::CartesianPoint(string _string)
	: x(0.0), y(0.0), z(0.0)
{
	if (_string == "X") {
		x = 1.0;
		return;
	}
	if (_string == "Y") {
		y = 1.0;
		return;
	}
	if (_string == "Z") {
		z = 1.0;
		return;
	}
}

CartesianPoint::CartesianPoint(double _x, double _y, double _z)
	: x(_x), y(_y), z(_z)
{
}

CartesianPoint::CartesianPoint(vector<double> _vec)
	: x(_vec[0]), y(_vec[1]), z(_vec[2])
{
	if (_vec.size() != 3) {
		cerr << "ERROR 2318: incorrect size of vector _vec != 3 in CartesianPoint::CartesianPoint(vector<double> _vec)" << endl;
		exit(2318);
	}
}

CartesianPoint::CartesianPoint(const CartesianPoint & _point) {
	x = _point.x;
	y = _point.y;
	z = _point.z;
}

CartesianPoint::~CartesianPoint() {
}

/*
void CartesianPoint::operator=(const CartesianPoint& _point) {
x = _point.x;
y = _point.y;
z = _point.z;
}
*/

/*
bool CartesianPoint::operator!=(const CartesianPoint& _point) const {
return ((x != _point.x) || (y != _point.y) || (z != _point.z));
}

bool CartesianPoint::operator==(const CartesianPoint &_point) const {
return ((x != _point.x) && (y != _point.y) && (z != _point.z));
}

void CartesianPoint::operator+=(const CartesianPoint& _point){
x += _point.x;
y += _point.y;
z += _point.z;
}

void CartesianPoint::operator-=(const CartesianPoint& _point){
x -= _point.x;
y -= _point.y;
z -= _point.z;
}

CartesianPoint CartesianPoint::operator-(const CartesianPoint &_point) const {
return CartesianPoint((x - _point.x), (y - _point.y), (z - _point.z));
}

CartesianPoint CartesianPoint::operator+(const CartesianPoint &_point) const {
return CartesianPoint((x + _point.x), (y + _point.y), (z + _point.z));
}
*/

CartesianPoint CartesianPoint::operator/(double _factor) const {
	if (_factor == 0.0){
		cerr << "ERROR 9141: Cannot divide by 0 in CartesianPoint CartesianPoint::operator/(double _factor) const" << endl;
		exit(9141);
	}
	return CartesianPoint((x / _factor), (y / _factor), (z / _factor));
}

/*
CartesianPoint CartesianPoint::operator*(double _factor) const {
return CartesianPoint((x*_factor), (y*_factor), (z*_factor));
}

double CartesianPoint::operator*(const CartesianPoint& _point) const {
return ((x*_point.x)+(y*_point.y)+(z*_point.z));
}

void CartesianPoint::operator*=(double _factor) {
x = x*_factor;
y = y*_factor;
z = z*_factor;
}
*/

void CartesianPoint::operator/=(double _factor) {
	if (_factor == 0.0){
		cerr << "ERROR 9141: Cannot divide by 0 in CartesianPoint CartesianPoint::operator/=(double _factor)" << endl;
		exit(9141);
	}
	x = x / _factor;
	y = y / _factor;
	z = z / _factor;
}

CartesianPoint CartesianPoint::operator*(const Matrix & _rotation) const {
	if (_rotation.getRow() != 3 || _rotation.getColumn() != 3) {
		cerr << "ERROR 9142: incorrect matrix size (" << _rotation.getRow() << "x" << _rotation.getColumn() << ") in CartesianPoint CartesianPoint::operator*(const Matrix _rotation) const" << endl;
		exit(9142);
	}
	CartesianPoint out = CartesianGeometry::matrixTimesCartesianPoint(*this, _rotation);
	return out;
}

void CartesianPoint::operator*=(const Matrix & _rotation) {
	if (_rotation.getRow() != 3 || _rotation.getColumn() != 3) {
		cerr << "ERROR 9143: incorrect matrix size (" << _rotation.getRow() << "x" << _rotation.getColumn() << ") in void CartesianPoint::operator*=(const Matrix _rotation)" << endl;
		exit(9143);
	}
	CartesianPoint out = CartesianGeometry::matrixTimesCartesianPoint(*this, _rotation);
	x = out.x;
	y = out.y;
	z = out.z;
}
double & CartesianPoint::operator[](size_t _n) {
	if (_n == 0) {
		return x;
	}
	else if (_n == 1) {
		return y;
	}
	else if (_n == 2) {
		return z;
	}
	cerr << "ERROR 2718: index " << _n << " out of range in double CartesianPoint::operator[](size_t _n)";
	exit(2718);
}

/*
string CartesianPoint::toString() const {
char c [100];
sprintf(c, "[%10.3f %10.3f %10.3f]", x, y, z);
return (string)c;
}
*/

/*
double CartesianPoint::getX() const {
return x;
}

double CartesianPoint::getY() const {
return y;
}

double CartesianPoint::getZ() const {
return z;
}

void CartesianPoint::setX(double _x) {
x = _x;
}

void CartesianPoint::setY(double _y) {
y = _y;
}

void CartesianPoint::setZ(double _z) {
z = _z;
}

void CartesianPoint::setCoor(double _x, double _y, double _z) {
x = _x;
y = _y;
z = _z;
}
*/

void CartesianPoint::setCoor(vector<Real> _vec) {
	if (_vec.size() < 3) {
		while (_vec.size() < 3) {
			_vec.push_back(0.0);
		}
	}
	x = _vec[0];
	y = _vec[1];
	z = _vec[2];
}

/*
void CartesianPoint::setCoor(const CartesianPoint & _vec) {
x = _vec.x;
y = _vec.y;
z = _vec.z;
}
*/

vector<Real> CartesianPoint::getCoor() const {
	vector<double> coordinates(3);
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
	return coordinates;
}

/*
CartesianPoint CartesianPoint::getUnit() const {
double dist = length();
return (*this/dist);
}

double CartesianPoint::length() const {
return sqrt(*this * *this);
}
*/

CartesianPoint CartesianPoint::cross(CartesianPoint _second) const
{
	double newx = y*_second.getZ() - z*_second.getY();
	double newy = z*_second.getX() - x*_second.getZ();
	double newz = x*_second.getY() - y*_second.getX();

	return CartesianPoint(newx, newy, newz);
}

double CartesianPoint::distance2(CartesianPoint _p) const {
	return CartesianGeometry::distance2(*this, _p);
}
double CartesianPoint::distance(CartesianPoint _p) const {
	return CartesianGeometry::distance(*this, _p);
}

double CartesianPoint::angle(CartesianPoint _p) const {
	return CartesianGeometry::angle(*this, _p);
}

double CartesianPoint::angle(CartesianPoint _center, CartesianPoint _third) const {
	return CartesianGeometry::angle(*this, _center, _third);
}

double CartesianPoint::angleRadians(CartesianPoint _p) const {
	return CartesianGeometry::angleRadians(*this, _p);
}

double CartesianPoint::angleRadians(CartesianPoint _center, CartesianPoint _third) const {
	return CartesianGeometry::angleRadians(*this, _center, _third);
}

double CartesianPoint::dihedral(CartesianPoint _second, CartesianPoint _third, CartesianPoint _fourth) const {
	return CartesianGeometry::dihedral(*this, _second, _third, _fourth);
}

