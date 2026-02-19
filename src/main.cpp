#include"../include/functions.h"

using namespace gal;
using namespace std;

int main() {


		
		Matrix A = Matrix::nullmatrix(4, 4);
		A[0][0] = 0; A[0][1] = 2; A[0][2] = 1; A[0][3] = 3;
		A[1][0] = 0; A[1][1] = 0; A[1][2] = 4; A[1][3] = 2;
		A[2][0] = 5; A[2][1] = 1; A[2][2] = 0; A[2][3] = 1;
		A[3][0] = 1; A[3][1] = 3; A[3][2] = 2; A[3][3] = 0;

		cout << A << endl;
}

