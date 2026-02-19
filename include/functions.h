#include"matrix.h"
#include <vector>

//DECLARATION OF FUNCTIONS

namespace gal {

	void print(const Matrix&); // prints the matrix
	void print(const std::vector<number>&);
	bool is_nullvector(const std::vector<number>&);
	Matrix clear_zero_numbers(const Matrix&);
	std::vector<number> operator+(const std::vector<number>&,const std::vector<number>&);
	std::vector<number> operator-(const std::vector<number>&,const std::vector<number>&);
	std::vector<number> operator-(const std::vector<number>&);
	Matrix operator-(const Matrix&,const std::vector<number>&);
	std::vector<number> operator*(number,const std::vector<number>&);
	std::vector<number> operator*(const std::vector<number>&, number);
	Matrix operator*(number,const Matrix&);
	Matrix operator*(const Matrix&, number);
	Matrix operator+(const Matrix&,const Matrix&);
	Matrix operator-(const Matrix&,const Matrix&);
	Matrix operator-(const Matrix&);
	Matrix operator/(const Matrix&,number);
	bool operator==(const Matrix& lhs, const Matrix& rhs);
	bool operator!=(const Matrix& lhs, const Matrix& rhs);
	bool operator==(const std::vector<number>& lhs, const std::vector<number>& rhs);
	bool operator!=(const std::vector<number>& lhs, const std::vector<number>& rhs);
	std::vector<number> operator/(const std::vector<number>&, number);
	double length(const std::vector<number>&);
	number inner_product(const std::vector<number>&, const std::vector<number>&); // returns the dot product of two vectors
	Matrix transpose(const Matrix&); // returns the transpose of the matrix
	Matrix multiply(const Matrix&,const Matrix&); // multiplies two matrices 
	Matrix multiply(const std::vector<number>&, const std::vector<number>&); // multiplies two vectors to get a matrix
	std::vector<number> multiply(const Matrix&,const std::vector<number>&); // multiplies a matrix and a vector
	Matrix row_echelon(const Matrix&); // returns the row echelon form
	Matrix reduced_row_echelon(const Matrix&); // returns the reduced row echelon form
	Matrix inverse(const Matrix&); // returns the inverse
	bool spaces_equal(const Matrix&, const Matrix&); // returns true if two vector sets span the same space
	bool vectors_colinear(const std::vector<number>&, const std::vector<number>&); // return true if two vectors are colinear
	Matrix nullspace(const Matrix&); // returns the basis vectors for the nullspace
	Matrix projection(const Matrix&); // returns the projection matrix 
	std::vector<number> householder_transform(const std::vector<number>&); // transforms a vector
	Matrix householder_matrix(const std::vector<number>&); // Constructs the Householder Matrix
	std::vector<Matrix> QR_decomposeGS(const Matrix&); // returns the QR decomposition of this matrix
	std::vector<Matrix> QR_decomposeHS(const Matrix&); // returns the QR decomposition of this matrix
	std::vector<Matrix> PLU_decompose(const Matrix&); // returns the LU decomposition of this matrix 
	std::vector<Matrix> SΛS_decompose(const Matrix&); // returns the SΛS decomposition of this matrix
	Matrix orthonormalise(const Matrix&); // returns the orthonormal matrix of bases spaning this column space
	std::vector<Matrix> solve_system(const Matrix&,const std::vector<number>&); // returns the solution to the system Ax = b (last column in the solution matrix is the particular solution) 
	std::vector<Matrix> best_solve(const Matrix&, const std::vector<number>&); // returns the best projected solution minimising error
	int rank_of_rref(const Matrix&); // return the rank of rref Matrix
	int rank(const Matrix&); // returns the rank of the matrix
	number determinant(const Matrix&);// returns the determinant
	std::vector<number> pivots(const Matrix&); // returns the list of pivots of this matrix
	std::vector<number> eigenvalues(const Matrix&); // returns the eigenvalues
	void transform_into_upper_hessenberg(Matrix&);
	void rotate_givens_similarity(unsigned int, unsigned int, unsigned int, Matrix&);
	void rotate_givens(unsigned int, unsigned int, unsigned int, Matrix&);
	Matrix eigenvectors(const Matrix&); // returns the eigenvectors of the matrix

}
