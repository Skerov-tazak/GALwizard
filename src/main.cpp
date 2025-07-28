#include"../include/functions.h"



int main() {
	
	

        Matrix A = Matrix::identity(3);
        Matrix result = orthonormalise(A).split_left(0); // Returns a vector with one matrix
														 
		print(result);


}
