#include"number.h"
#include <cmath>

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

//{ Basic Functions

double number::angle()const{

	if(is_zero(imaginary)){
		
		if(is_zero(real) || real > 0)
			return 0.0;
		else 
			return M_PI;

	}
	else if(is_zero(real)){

		return M_PI/2;
	}
	else 
		return atan(imaginary/real);
}

number from_polar(double r, double t){
	
	return number(r * cos(t), r * sin(t));
}

std::pair<double, double> to_polar(const number& a){
	
	return {a.radius(), a.angle()};
}

number number::copy()const{
	
	return number(real, imaginary); 
}

number conj(const number& a){

	return number(a.real, -a.imaginary);
}

number number::normalised()const {

    double r = std::hypot(real, imaginary);
    
	if (is_zero(r))
		return number(0, 0);
    
	return number(real / r, imaginary / r);
}

number number::rotate(double theta)const {

    return *this * from_polar(1.0, theta);
}

//}

//{ Advanced Functions 

number exp(const number& a)
{
    double exp_real = std::exp(a.real);

    return number(exp_real * std::cos(a.imaginary), exp_real * std::sin(a.imaginary));
}

number ln(const number& a)
{
	return number(std::log(a.radius()), a.angle());
}

number number::power(number x)const{

	if(x == 1)
		return copy();
	
	if(x == 2)
		return copy()*copy();

	return exp((ln(*this) * x));
}

number cos(const number& a){
 
	return number(std::cos(a.real) * std::cosh(a.imaginary), -std::sin(a.real) * std::sinh(a.imaginary));
}

number sin(const number& a){
	
	return number(std::sin(a.real) * std::cosh(a.imaginary),std::cos(a.real) * std::sinh(a.imaginary));
}

number tan(const number& a){
	
	return sin(a)/cos(a);
}

std::vector<number> number::roots(unsigned int n)const{

	std::vector<number> roots;
	double r = std::pow(radius(), 1.0/n);
	double t = angle();

	for (unsigned int k = 0; k < n; k++)
		roots.emplace_back(from_polar(r, (t + 2*M_PI*k) / n));

	return roots;
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

	a = (b * conj(a))/r;

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
	
	a = (a * conj(b))/r;
 
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

bool operator==(const number& a, double b){

	return (approx_equal(a.real, b) && is_zero(a.imaginary));
}

bool operator==(double b, const number& a){
	
	return (a == b);
}

bool operator!=(const number& a, const number& b){

	return !(a == b);
}

bool operator!=(const number& a, double b){

	return !(a == b);
}

bool operator!=(double b, const number& a){
	
	return !(b == a);
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

	if(number.im() >= 0){
		os << round_for_zero(number.re()) << " + " << round_for_zero(number.im()) << "i";
	}else {
		os << round_for_zero(number.re()) << " - " << round_for_zero(-number.im()) << "i";
	}
	return os;

}


//}


