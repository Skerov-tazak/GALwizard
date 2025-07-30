#include<cmath>
#include <ostream>

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

		inline double angle()const {return (real == 0.0 && imaginary == 0 ? 0.0 : atan(imaginary/real));} // returns the angle 

		constexpr inline double radius_squared()const noexcept {return imaginary * imaginary + real * real;} // returns the radius 

		inline double radius()const {return std::hypot(imaginary, real);}

		number power(double)const; // returns the number to any real power 

		number copy()const;
		
		number conj()const;
		
		number normalised()const;

		number rotate(double)const;

		number& operator+=(const number&);

		number& operator/=(const number&);

		number& operator-=(const number&);

		number& operator*=(const number&);

		number& operator*=(double);

		number& operator/=(double);

		number& operator-=(double);
	
		number& operator+=(double);

		friend bool operator==(const number&, const number&);

		friend bool operator!=(const number&, const number&);

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
	
		static number from_polar(double, double);
		
		static std::pair<double, double> to_polar(const number&); // returns the euler form as std::pair 
																  
		static inline number i() {return number(0.0, 1.0);}

};


bool approx_equal(double a, double b);

bool is_zero(double number);

std::ostream& operator<<(std::ostream&, const number&);


