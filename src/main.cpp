#include"../include/functions.h"



int main() {


	
		Matrix A = Matrix::nullmatrix(2, 2);
        A[0][0] = 4; A[0][1] = 7;
        A[1][0] = 2; A[1][1] = 6;

		print(A);
		print(inverse(A));
		
		print(multiply(A,inverse(A)));

}
