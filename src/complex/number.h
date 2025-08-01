#include<cmath>
#include<stdexcept>
#include<vector>
#include <ostream>

constexpr double EPS = 1e-6;

constexpr inline bool approx_equal(double a, double b) {return std::fabs(a - b) < EPS;}

constexpr inline bool is_zero(double x) {return approx_equal(x, 0);}

constexpr inline double round_for_zero(double x) {if(is_zero(x)){return 0.0;}else{return x;}}

class number{

	private:

		double imaginary = 0.0;

		double real = 0.0;

	public:

		number();

		number(double);

		number(double, double);

		constexpr inline void re(double r) noexcept {real = r;} // sets real part

		constexpr inline void im(double i) noexcept {imaginary = i;} // sets imaginary part 

		constexpr inline double re()const noexcept {return real;} // returns the real part

		constexpr inline double im()const noexcept {return imaginary;} // returns the imaginary part 

		constexpr inline bool is_real() const { return is_zero(imaginary); }

		constexpr inline bool is_imag() const { return is_zero(real); }

		constexpr inline double radius_squared()const noexcept {return imaginary * imaginary + real * real;} // returns the radius 

		inline double radius()const {return std::hypot(imaginary, real);}

		std::vector<number> roots(unsigned int)const;

		number power(number)const; // returns the number to any complex power 

		double angle()const; // returns the angle 

		number copy()const;
		
		number normalised()const;

		number rotate(double)const;

		number& operator+=(const number&);

		number& operator/=(const number&);

		number& operator-=(const number&);

		number& operator*=(const number&);

		explicit operator bool() const noexcept { return (*this == 0.0); }

		number& operator*=(double);

		number& operator/=(double);

		number& operator-=(double);
	
		number& operator+=(double);

		friend bool operator==(const number&, const number&);

		friend bool operator==(const number&, double);

		friend bool operator==(double, const number&);

		friend bool operator!=(const number&, const number&);

		friend bool operator!=(const number&, double);

		friend bool operator!=(double, const number&);

		friend number operator*(number, const number&);

		friend number operator*(double, number);

		friend number operator*(number, double);

		friend number operator/(number, const number&);

		friend number operator/(number, double);

		friend number operator/(double, number);

		friend number operator+(number, const number&);

		friend number operator+(double, number);

		friend number operator+(number, double);

		friend number operator-(number, const number&);

		friend number operator-(double, number);

		friend number operator-(number, double);

		friend number operator-(const number&);
		
		friend number conj(const number&);

		friend number from_polar(double, double);

		friend std::pair<double, double> to_polar(const number&); // returns the euler form as std::pair 
	
		friend number ln(const number&);

		friend number exp(const number&);

		friend number cos(const number&);
		
		friend number sin(const number&);

		friend number tan(const number&);

};

inline number i() {return number(0.0, 1.0);}

std::ostream& operator<<(std::ostream&, const number&);


