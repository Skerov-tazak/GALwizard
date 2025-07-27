#include "matrix.h"
#include <limits>
#include <vector>




void test_basics()
{


}

void test_vect_multiply() 
{
	std::vector<float> difference {1 - length({1,0,-1}), 0, -1};
	Matrix test = Matrix::identity(3);
	float len = length(difference);	

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++){

			test[j][i] = test[j][i] - 2* (difference[j] * difference[i]) / (len * len)  ;
		}
	}

	print(test);


}


void test_QR(){
	
	Matrix test;
	test.push_back_row({1,0,1});
	test.push_back_row({0,1,2});
	test.push_back_row({-1,-1,1});

	std::vector<Matrix> QR = QR_decomposeHS(test);

	print(QR[0]);
	print(QR[1]);

	QR = QR_decomposeGS(test);

	print(QR[0]);
	print(QR[1]);

}

int main()
{
	test_QR();
}
