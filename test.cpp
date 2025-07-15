#include "matrix.h"
#include <limits>
void test_basics();


int main()
{
	test_basics();	
}

void test_basics()
{
	Matrix test;
	test.push_back_row({1,0,1,4});
	test.push_back_row({0,1,2,-1});
	test.push_back_row({-1,-1,1,-2});
	test.push_back_row({-1,-2,1,-1});
	print(test);
	print(row_echelon(test));

}



