#include"complex.h"
#include <cmath>
#include <stdexcept>

constexpr double EPS = 1e-6;

//{ Constructors

number::number(){

}

number::number(double a){

	real = a;
	imaginary = 0.0;
}

number::number(double a, double b){

	real = a;
	imaginary = b;
}

//}

//{ Utility Functions

bool approx_equal(double a, double b){
	
	return std::fabs(a - b) < EPS;
}

bool is_zero(double number)
{
	return approx_equal(number, 0);
}

//}

//{ Basic Member Functions

number number::from_polar(double r, double t){
	
	return number(r * cos(t), r * sin(t));
}

std::pair<double, double> number::to_polar(const number& a){
	
	return {a.radius(), a.angle()};
}

number number::copy()const{
	
	return number(real, imaginary); 
}

number number::conj()const{

	return number(real, -imaginary);
}

number number::power(double x)const{

	if(x == 2){
		return copy();
	}

	std::pair<double,double> euler_form = to_polar(*this);
	euler_form.second = euler_form.second * x;
	euler_form.first = pow(euler_form.first, x);
	return from_polar(euler_form.first, euler_form.second);
}

number number::normalised()const {

    double r = std::hypot(real, imaginary);
    
	if (is_zero(r))
		return number(0, 0);
    
	return number(real / r, imaginary / r);
}

number number::rotate(double theta)const {

    return *this * from_euler(1.0, theta);
}



//}

//{ Operators

number operator*(number a, const number& b){
	
	a *= b;
	return a;
}

number operator*(double b, number a){
	
	a.real *= b;
	a.imaginary *= b;
	return a;
}

number operator*(number a, double b){

	return b * a;
}

number operator/(double b, number a){

	double r = a.radius_squared();
	
	if(is_zero(r))
		throw std::invalid_argument("Division by zero not allowed");

	a = (b * a.conj())/r;

	return a; 
}

number operator/(number a, double b){
	
	if(is_zero(b))
		throw std::invalid_argument("Division by zero not allowed");

	a.real /= b;
	a.imaginary /= b;
	return a;
}

number operator/(number a, const number& b){

	double r = b.radius_squared();
	
	if(is_zero(r))
		throw std::invalid_argument("Division by zero not allowed");
	
	a = (a * b.conj())/r;
 
	return a;
}

number operator+(number a, const number& b){
	
	a.real += b.real;
	a.imaginary += b.imaginary;
	return a;
}

number operator+(number a, double b){

	a.real += b;
	return a;
}

number operator+(double b, number a){
	
	return a + b;
}

number operator-(number a, const number& b){

	a.real -= b.real;
	a.imaginary -= b.imaginary;
	return a; 

}

number operator-(number a, double b){
	
	a.real -= b;
	return a;
}

number operator-(double b, number a){

	a.real = b - a.re();
	a.imaginary = -a.im();	
	return a;
}

number operator-(const number& a){
	
	return number(-a.re(), -a.im());  
}

bool operator==(const number& a, const number& b){

	return (approx_equal(a.re(), b.re()) && approx_equal(a.im(), b.im()));

}

bool operator-=(const number& a, const number& b){

	return !(a == b);
}

number& number::operator+=(const number& rhs) {
	
	real += rhs.real;
	imaginary += rhs.imaginary;

	return *this; 
}

number& number::operator-=(const number& rhs) {
	
	real -= rhs.real;
	imaginary -= rhs.imaginary;
	
	return *this;
}

number& number::operator*=(const number& rhs) {
	
	double new_real = real * rhs.real - imaginary * rhs.imaginary;
	double new_imag = real * rhs.imaginary + imaginary * rhs.real;
	real = new_real;
	imaginary = new_imag;
	
	return *this;
}

number& number::operator/=(const number& rhs) {

	double denom = rhs.radius_squared();
	double new_real = (real * rhs.real + imaginary * rhs.imaginary) / denom;
	double new_imag = (imaginary * rhs.real - real * rhs.imaginary) / denom;
	real = new_real;
	imaginary = new_imag;
	
	return *this;
}

number& number::operator*=(double s) {

    real *= s;
    imaginary *= s;
    return *this;
}

number& number::operator/=(double s) {

    real /= s;
    imaginary /= s;
    return *this;
}

number& number::operator+=(double s){

	real += s;
	return *this;
}

number& number::operator-=(double s){

	real -= s;
	return *this;
}

std::ostream& operator<<(std::ostream& os, const number& number){

	os << number.re() << " + " << number.im() << "i";
	return os;

}
//}


