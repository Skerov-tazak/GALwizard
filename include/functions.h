#include"matrix.h"

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
bool operator==(const Matrix& lhs, const Matrix& rhs);
bool operator!=(const Matrix& lhs, const Matrix& rhs);
bool operator==(const std::vector<float>& lhs, const std::vector<float>& rhs);
bool operator!=(const std::vector<float>& lhs, const std::vector<float>& rhs);
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


