#pragma once
#include<vector>
#include<iostream>
//{ MATRIX CLASS DECLARATION
class Matrix
{
	private:
		std::vector<std::vector<float>> data;
		
		Matrix(unsigned int); // Identity Matrix constructor
							  
		Matrix(unsigned int,unsigned int); // Nullmatrix of m by n constructor
	
	public:
		Matrix(); // default constructor
	
		Matrix(std::vector<std::vector<float>>); // constructs a full Matrix (ROW BY ROW)
		
		int rows(); // Returns number of rows
		
		int cols(); // Returns number of columns
	
		void push_back_col(std::vector<float>); // adds a column at the back
		
		void push_back_row(std::vector<float>); // adds a row at the back
		
		void swap_rows(unsigned int,unsigned int); // swaps two rows

		std::vector<float>& operator[](unsigned int); // Operator overload for getting data

		void append(Matrix); // Appends another matrix on the right of this one

		Matrix split_right(unsigned int); // returns submatrix to the right of this index (Included)

		Matrix split_left(unsigned int); // returns submatrix to the left of this index (Included);

		static Matrix identity(unsigned int); // returns an n by n identity matrix
									  
		static Matrix nullmatrix(unsigned int, unsigned int); // return an n by m nullmatrix 

		std::vector<float> get_row(unsigned int); // returns a row in vector form

		std::vector<float> get_col(unsigned int); // returns a column in vector form
};
//}

//DECLARATION OF FUNCTIONS

void print(Matrix); // prints the matrix
void print(std::vector<float>);
bool zero(float);
bool is_nullvector(std::vector<float>);
float test_zero(float);
Matrix clear_zero_floats(Matrix);
std::vector<float> operator+(std::vector<float>,std::vector<float>);
std::vector<float> operator-(std::vector<float>,std::vector<float>);
std::vector<float> operator-(std::vector<float>);
Matrix operator-(Matrix,std::vector<float>);
std::vector<float> operator*(float,std::vector<float>);
Matrix operator*(float,Matrix);
Matrix operator+(Matrix,Matrix);
Matrix operator-(Matrix,Matrix);
Matrix operator-(Matrix);
Matrix operator/(Matrix,float);
std::vector<float> operator/(std::vector<float>, float);
float length(std::vector<float>);
float inner_product(std::vector<float>,std::vector<float>); // returns the dot product of two vectors
Matrix transpose(Matrix); // returns the transpose of the matrix
Matrix multiply(Matrix,Matrix); // multiplies two matrices 
Matrix multiply(std::vector<float>, std::vector<float>); // multiplies two vectors to get a matrix
std::vector<float> multiply(Matrix,std::vector<float>); // multiplies a matrix and a vector
Matrix row_echelon(Matrix); // returns the row echelon form
Matrix reduced_row_echelon(Matrix); // returns the reduced row echelon form
Matrix inverse(Matrix); // returns the inverse
Matrix nullspace(Matrix); // returns the basis vectors for the nullspace
Matrix projection(Matrix); // returns the projection matrix 
std::vector<float> householder_transform(std::vector<float>); // transforms a vector
Matrix householder_matrix(std::vector<float>); // Constructs the Householder Matrix
std::vector<Matrix> QR_decomposeGS(Matrix); // returns the QR decomposition of this matrix
std::vector<Matrix> QR_decomposeHS(Matrix); // returns the QR decomposition of this matrix
std::vector<Matrix> PLU_decompose(Matrix); // returns the LU decomposition of this matrix 
std::vector<Matrix> SΛS_decompose(Matrix); // returns the SΛS decomposition of this matrix
Matrix orthonormalise(Matrix); // returns the orthonormal matrix of bases spaning this column space
Matrix solve_system(Matrix,std::vector<float>); // returns the solution to the system Ax = b (last column in the solution matrix is the particular solution) 
Matrix best_solve(Matrix, std::vector<float>); // returns the best projected solution minimising error
int rank(Matrix); // returns the rank of the matrix
float determinant(Matrix);// returns the determinant
std::vector<float> pivots(Matrix); // returns the list of pivots of this matrix
std::vector<float> eigenvalues(Matrix); // returns the eigenvalues
Matrix eigenvectors(Matrix); // returns the eigenvectors of the matrix


