#include"matrix.h"
#include <vector>

//DECLARATION OF FUNCTIONS


void print(Matrix); // prints the matrix
void print(std::vector<number>);
bool is_nullvector(std::vector<number>);
Matrix clear_zero_numbers(Matrix);
std::vector<number> operator+(std::vector<number>,std::vector<number>);
std::vector<number> operator-(std::vector<number>,std::vector<number>);
std::vector<number> operator-(std::vector<number>);
Matrix operator-(Matrix,std::vector<number>);
std::vector<number> operator*(number,std::vector<number>);
std::vector<number> operator*(std::vector<number>, number);
Matrix operator*(number,Matrix);
Matrix operator*(Matrix, number);
Matrix operator+(Matrix,Matrix);
Matrix operator-(Matrix,Matrix);
Matrix operator-(Matrix);
Matrix operator/(Matrix,number);
bool operator==(const Matrix& lhs, const Matrix& rhs);
bool operator!=(const Matrix& lhs, const Matrix& rhs);
bool operator==(const std::vector<number>& lhs, const std::vector<number>& rhs);
bool operator!=(const std::vector<number>& lhs, const std::vector<number>& rhs);
std::vector<number> operator/(std::vector<number>, number);
double length(std::vector<number>);
number inner_product(std::vector<number>,std::vector<number>); // returns the dot product of two vectors
Matrix transpose(Matrix); // returns the transpose of the matrix
Matrix multiply(Matrix,Matrix); // multiplies two matrices 
Matrix multiply(std::vector<number>, std::vector<number>); // multiplies two vectors to get a matrix
std::vector<number> multiply(Matrix,std::vector<number>); // multiplies a matrix and a vector
Matrix row_echelon(Matrix); // returns the row echelon form
Matrix reduced_row_echelon(Matrix); // returns the reduced row echelon form
Matrix inverse(Matrix); // returns the inverse
bool spaces_equal(Matrix, Matrix); // returns true if two vector sets span the same space
bool vectors_colinear(std::vector<number>, std::vector<number>); // return true if two vectors are colinear
Matrix nullspace(Matrix); // returns the basis vectors for the nullspace
Matrix projection(Matrix); // returns the projection matrix 
std::vector<number> householder_transform(std::vector<number>); // transforms a vector
Matrix householder_matrix(std::vector<number>); // Constructs the Householder Matrix
std::vector<Matrix> QR_decomposeGS(Matrix); // returns the QR decomposition of this matrix
std::vector<Matrix> QR_decomposeHS(Matrix); // returns the QR decomposition of this matrix
std::vector<Matrix> PLU_decompose(Matrix); // returns the LU decomposition of this matrix 
std::vector<Matrix> SΛS_decompose(Matrix); // returns the SΛS decomposition of this matrix
Matrix orthonormalise(Matrix); // returns the orthonormal matrix of bases spaning this column space
std::vector<Matrix> solve_system(Matrix,std::vector<number>); // returns the solution to the system Ax = b (last column in the solution matrix is the particular solution) 
std::vector<Matrix> best_solve(Matrix, std::vector<number>); // returns the best projected solution minimising error
int rank_of_rref(Matrix); // return the rank of rref Matrix
int rank(Matrix); // returns the rank of the matrix
number determinant(Matrix);// returns the determinant
std::vector<number> pivots(Matrix); // returns the list of pivots of this matrix
std::vector<number> eigenvalues(Matrix); // returns the eigenvalues
Matrix eigenvectors(Matrix); // returns the eigenvectors of the matrix


