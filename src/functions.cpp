#include "../include/functions.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iomanip> // for std::fixed and std::setprecision

//{ Print functions 
void print(Matrix to_print)
{
	std::cout << std::fixed << std::setprecision(6);
    const int fieldWidth = 13; 
    
    for (size_t i = 0; i < to_print.rows(); ++i) {
        for (size_t j = 0; j < to_print.cols(); ++j) {
          
			std::cout << std::setw(fieldWidth) << to_print[i][j];   
			if (j + 1 < to_print.cols()) std::cout << " ";

        }

        std::cout << '\n';
    }
	std::cout << '\n';
}

void print(std::vector<number> to_print)
{
	std::cout << std::fixed << std::setprecision(6);
    
    const int fieldWidth = 13; 
    
    for (const number& val : to_print) 
        std::cout << std::setw(fieldWidth) << val << '\n';
}
//}

//{ Null validation functions

bool is_nullvector(std::vector<number> vect)
{
	for(number x : vect )
	{
		if(x != 0.0)
			return false;
	}
	return true;
}

number test_zero(number number)
{
	if(number == 0)
		return 0;
	else 
		return number;

}

Matrix clear_zero_numbers(Matrix to_clear)
{
	for(int i = 0; i < to_clear.rows(); i++){
		for(int j = 0; j < to_clear.cols(); j++){
			to_clear[i][j] = test_zero(to_clear[i][j]);
		}
	}
	return to_clear;
}


//}

//{ Operators
std::vector<number> operator+(std::vector<number> one, std::vector<number> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("operator+(): undefined addition of different length vectors");
	std::vector<number> answer(one.size(),0);

	for(int i = 0; i < one.size(); i++)
		answer[i] = one[i] + two[i];
	
	return answer;
}


std::vector<number> operator-(std::vector<number> one, std::vector<number> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("operator+(): undefined addition of different length vectors");
	std::vector<number> answer(one.size(),0);

	for(int i = 0; i < one.size(); i++)
		answer[i] = one[i] - two[i];
	
	return answer;
}

Matrix operator+(Matrix one, Matrix two)
{
	if(one.rows() != two.rows() || one.cols() != two.cols())
		throw std::invalid_argument("operator+(): undefined addition of two different size matrices");
	
	Matrix answer = Matrix::nullmatrix(one.rows(),one.cols());

	for(int i = 0; i < one.rows(); i++)
	{
		for(int j = 0; j < one.cols(); j++)
			answer[i][j] = one[i][j] + two[i][j];
	}

	return answer;
}

Matrix operator-(Matrix one, Matrix two)
{
	if(one.rows() != two.rows() || one.cols() != two.cols())
		throw std::invalid_argument("operator+(): undefined addition of two different size matrices");
	
	Matrix answer = Matrix::nullmatrix(one.rows(),one.cols());

	for(int i = 0; i < one.rows(); i++)
	{
		for(int j = 0; j < one.cols(); j++)
			answer[i][j] = one[i][j] - two[i][j];
	}

	return answer;
}

Matrix operator-(Matrix one)
{
	for(int i = 0; i < one.rows(); i++)
	{
		for(int j = 0; j < one.cols(); j++)
		{
			one[i][j] = -one[i][j];
		}
	}
	return one;
}

std::vector<number> operator-(std::vector<number> one)
{
	for(int i = 0; i < one.size(); i++)
	{
		one[i] = -one[i];
	}
	return one;
}

Matrix operator-(Matrix left, std::vector<number> right)
{
	if(left.cols() > 1)
		throw(std::invalid_argument("operator-(): operation undefined"));

	for(int i = 0; i < left.rows(); i++)
		left[i][0] = left[i][0] - right[i];
	return left;
}

std::vector<number> operator*(number scalar, std::vector<number> vect)
{
	for(auto x : vect)
		x = x * scalar;
	return vect;
}

std::vector<number> operator*(std::vector<number> vect, number scalar)
{
	return scalar*vect;
}

std::vector<number> operator/(std::vector<number> vect, number scalar) 
{
	for(auto x : vect) 
		x = x / scalar;
	return vect;
}

Matrix operator/(Matrix matrix, number scalar)
{
	number real = 1 / scalar;
	return real * matrix;
}

Matrix operator*(number scalar,Matrix matrix)
{
	for(int i = 0; i < matrix.cols();i++)
	{
		for(int j = 0; j < matrix.rows(); j++)
		{
			matrix[j][i] = matrix[j][i] * scalar;
		}
	}

	return matrix;

}

Matrix operator*(Matrix matrix, number scalar){
	return scalar * matrix;
}

bool operator==(const Matrix& lhs, const Matrix& rhs){
	
	if(lhs.rows() != rhs.rows() || rhs.cols() != lhs.cols())
		return false;

	for(int i = 0; i < lhs.cols(); i++){
		for(int j = 0; j < lhs.rows(); j++){
			if(lhs.at(j,i) != rhs.at(j,i))
				return false;
		}
	}

	return true;

}

bool operator!=(const Matrix& lhs, const Matrix& rhs){

	return !(lhs == rhs);
}

bool operator==(const std::vector<number>& lhs, const std::vector<number>& rhs){

	if(lhs.size() != rhs.size())
		return false;
	
	for(int i = 0; i < lhs.size(); i++){
		if(lhs.at(i) != rhs.at(i))
			return false;		
	}
	
	return true;
}

bool operator!=(const std::vector<number>& lhs, const std::vector<number>& rhs){
	
	return !(lhs == rhs);
}



//}

//{ Length

double length(std::vector<number> vect)
{
	double length = 0;

	for(number x : vect)
	{
		length += x.radius_squared();
	}

	length = sqrt(length);
	return length;
}

void normalise(std::vector<number> vect) 
{
	double len = length(vect);

	for(unsigned int i = 0; i < vect.size(); i++) 
	{
		vect[i] = vect[i] / len;
	}		
}

//}

//{ Inner products

number inner_product(std::vector<number> one, std::vector<number> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("inner_product(): undefined product of different length vectors");
	
	number answer = 0;
	for(int i = 0; i < one.size(); i++)
	{
		answer += one[i] * two[i]; 
	}
	return answer;
}
//}

//{ Transpose
Matrix transpose(Matrix original)
{
	Matrix transposed = Matrix::nullmatrix(original.cols(),original.rows());

	for(int i = 0; i < original.rows(); i++){
		for(int j = 0; j < original.cols(); j++){
			transposed[j][i] = original[i][j];
		}
	}

	return transposed;
}
//}

//{ Multiply
Matrix multiply(Matrix left,Matrix right)
{
	if(left.cols() != right.rows())
		throw std::invalid_argument("multiply(): this multiplication is not defined");

	Matrix answer = Matrix::nullmatrix(left.rows(),right.cols());

	for(int i = 0; i < answer.rows(); i++){
		for(int j = 0; j < answer.cols(); j++){
			for(int p = 0; p < left.cols(); p++){
				answer[i][j] = test_zero(answer[i][j] + left[i][p]*right[p][j]);
			}
		}
	}
	return answer;
}

std::vector<number> multiply(Matrix left,std::vector<number> right)
{
	if(left.cols() != right.size())
		throw std::invalid_argument("multiply(): invalid vector length");

	std::vector<number> answer(left.rows(), 0);

	for(int i = 0; i < left.rows(); i++){
		for(int j = 0; j < left.cols(); j++){
			answer[i] = test_zero(answer[i] + left[i][j] * right[j]);
		}
	}

	return answer;
}

Matrix multiply(std::vector<number> left, std::vector<number> right) // treats left as vertical, right as horizontal
{
	Matrix answer = Matrix::nullmatrix(left.size(), right.size());
	
	for(int i = 0; i < answer.rows(); i++){
		for(int j = 0; j < answer.cols(); j++){

			answer[i][j] = left[i] * right[j];
		}
	}

	return answer;
}

//}

//{ Row Echelon form

void clear_pivot_column(Matrix& matrix, unsigned int row, unsigned int col){

	number pivot = matrix[row][col];
	for(int i = row + 1; i < matrix.rows(); i++) // for all rows
	{
		number factor = matrix[i][col]/pivot; // finds a factor

		for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
		{
			matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
		}
	}
}

int find_nonzero_pivot(const Matrix& matrix, unsigned int row, unsigned int col){

		int temp_row = row;
		while(temp_row < matrix.rows() && (0.0 == matrix.at(temp_row, col)))
			temp_row++;

		return temp_row;
}


Matrix row_echelon(Matrix matrix)
{
	int col = 0;
	int row = 0;

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		
		unsigned int pivot_row = find_nonzero_pivot(matrix, row, col);

		if(pivot_row < matrix.rows()){
			
			matrix.swap_rows(pivot_row, row); // if there is one, swaps rows

			clear_pivot_column(matrix, row, col);

		row++;
		}

	col++;
	}

	return matrix;
}
//}

//{ Reduced row echelon form

void clear_over_pivots(Matrix& matrix, std::vector<std::pair<unsigned int, unsigned int>>& pivot_indeces){

	for(int i = pivot_indeces.size() - 1; i >= 0; i--) 	// for all pivots
	{ 

		number pivot = matrix.at(pivot_indeces.at(i).first,pivot_indeces[i].second);

		for(int j = pivot_indeces.at(i).first - 1; j >= 0; j--) // for all rows above this pivot
		{
			number factor = matrix.at(j,pivot_indeces.at(i).second)/pivot;	

			for(int p = matrix.cols() - 1; p >= 0 && (unsigned int)p >= pivot_indeces[i].second; p--) // for all columns in that row
			{
				matrix[j][p] = test_zero(matrix.at(j,p) - factor * matrix.at(pivot_indeces.at(i).first, p));
			}
		}
		for(int g = matrix.cols() - 1; g >= 0 && (unsigned int)g >= pivot_indeces.at(i).second; g--)
			matrix[pivot_indeces.at(i).first][g] = matrix.at(pivot_indeces.at(i).first,g) / pivot;
	}
}

void ref_remember_pivots(Matrix& matrix, std::vector<std::pair<unsigned int, unsigned int>>& pivot_indeces){

	int col = 0;
	int row = 0;

	while(col < matrix.cols() && row < matrix.rows())
	{
		unsigned int pivot_row = find_nonzero_pivot(matrix, row, col);

		if(pivot_row != matrix.rows()) // if none is found, skips this column
		{
			matrix.swap_rows(pivot_row, row); // if there is one, swaps rows
			number pivot = matrix[row][col];
			pivot_indeces.push_back({row,col});

			clear_pivot_column(matrix, row, col);

			row++;
		}	
		col++;
	}
}

Matrix reduced_row_echelon(Matrix matrix)
{	
	std::vector<std::pair<unsigned int, unsigned int>> pivot_indeces; // keep track of pivots to clear 

	// finds the next non-zero pivot

	ref_remember_pivots(matrix, pivot_indeces);

	// clears the upper rows too

	clear_over_pivots(matrix, pivot_indeces);

	return matrix;
}

//}

//{ Inverse
Matrix inverse(Matrix matrix)
{
	if(matrix.cols() != matrix.rows()) // check if the matrix is square
		throw std::invalid_argument("inverse(): the matrix isn't square");

	Matrix identity = Matrix::identity(matrix.cols());
	matrix.append(identity);
	matrix = reduced_row_echelon(matrix);

	for(int i = 0; i < matrix.rows(); i++) // check if this matrix is singular
	{
		if(matrix[i][i] == 0)
			throw std::invalid_argument("inverse(): this matrix is singular and has no inverse");
	}

	Matrix answer = Matrix::empty();
	for(int i = matrix.cols()/2; i < matrix.cols(); i++) // create the answer;
		answer.push_back_col(matrix.get_col(i));

	return answer;
}
//}

//{ Compare Spans 

bool spaces_equal(Matrix A, Matrix B)
{
	
	std::vector<std::pair<unsigned int, unsigned int>> pivots;

	ref_remember_pivots(A, pivots);
	unsigned int rank_A = pivots.size();
	clear_over_pivots(A, pivots);

	pivots.resize(0);
	ref_remember_pivots(B, pivots);
	unsigned int rank_B = pivots.size();
	clear_over_pivots(B, pivots);

	if(rank_A != rank_B)
		return false;

	for(int i = 0; i < A.cols(); i++){
		if(solve_system(B, A.get_col(i)).empty()){
			return false;	
		}
	}
	for(int j = 0; j < B.cols(); j++){
		if(solve_system(A, B.get_col(j)).empty())
			return false;
	}

	return true;
}

bool vectors_colinear(std::vector<number> v, std::vector<number> w){
	
	if(v.size() != w.size())
		return false;

	Matrix together = Matrix::empty();

	together.push_back_col(v);
	together.push_back_col(w);
	
	if(rank(together) > 1)
		return false;

	return true;
}
//}

//{ Nullspace
Matrix nullspace(Matrix matrix)
{
	// First get to the reduced row echelon form and remember pivots

	std::vector<std::pair<unsigned int, unsigned int>> pivot_indeces; // keep track of pivots to clear 

	// clears under the pivots
	
	ref_remember_pivots(matrix, pivot_indeces);

	// clears the upper rows too
	
	clear_over_pivots(matrix, pivot_indeces);

	// After getting to rref, finally the nullspace is generated
	Matrix answer = Matrix::empty();

	std::vector<number> null_vector(matrix.cols(),0);
	answer.push_back_col(null_vector);


	int next_pivot = 0;
	for(int i = 0; i < matrix.cols(); i++) // iterate through, skipping pivot columns
	{
		std::vector<number> partial_answer(matrix.cols(),0);
		while(next_pivot < pivot_indeces.size() && pivot_indeces[next_pivot].second == i)
		{
			i++;
			next_pivot++;
			if(i == matrix.cols()) // in case the last few cols are pivot columns
				return answer;
		}
		partial_answer[i] = 1;

		for(int j = 0; j < pivot_indeces.size(); j++)
		{
			partial_answer[pivot_indeces[j].second] = - matrix[j][i];
		}
		answer.push_back_col(partial_answer);
	}

	return answer;
}

//}

//{ Projections
Matrix projection(Matrix matrix)
{
	// uses the formula A*((At*A)^-1)*At (where At is the transpose of A, and ^-1 is the inverse)
	Matrix transposed = transpose(matrix); // 
	Matrix middle_term = inverse(multiply(transposed,matrix));
	return multiply(multiply(matrix,middle_term),transposed);
}

//}

//{ Givens Rotations, Hessenberg forms 

void transform_into_upper_hessenberg(Matrix& matrix){

	for(int i = 0; i < matrix.cols() - 2; i++){
		for(int j = i + 2; j < matrix.rows(); j++){
			
			rotate_givens_similarity(j, i + 1, i, matrix);
		}
	}
} 

void rotate_givens_similarity(unsigned int row_zero, unsigned int row_len, unsigned int col, Matrix& matrix) // Perform givens rotation
{		   																								 // While right multiplying by transpose

	if(row_zero == row_len)
		throw std::invalid_argument("The row that is zeroed must be different from the row that takes the length");

	number a = matrix.at(row_zero, col);
	number b = matrix.at(row_len, col);

	double normalisation = std::sqrt(a.radius_squared() + b.radius_squared());

	number c_factor = b/normalisation; 
	number s_factor = a/normalisation;

		

	for(int i = 0; i < matrix.cols(); i++){


		number temp = matrix.at(row_zero, i);
		matrix.at(row_zero, i) = test_zero(c_factor * matrix.at(row_zero,i) - s_factor * matrix.at(row_len,i));
		matrix.at(row_len, i) = c_factor * matrix.at(row_len,i) + s_factor * temp; 
	}	

	// Here we transpose the givens matrix and right multiply 
	
	c_factor = conj(c_factor);
	s_factor = conj(s_factor);

	for(int i = 0; i < matrix.rows(); i++){

		number temp = matrix.at(i, row_zero);
		matrix.at(i, row_zero) = c_factor * matrix.at(i, row_zero) - s_factor * matrix.at(i, row_len); 
		matrix.at(i, row_len) = c_factor * matrix.at(i, row_len) + s_factor *  temp;
	}


}

void rotate_givens(unsigned int row_zero, unsigned int row_len, unsigned int col, Matrix& matrix){

	if(row_zero == row_len)
		throw std::invalid_argument("The row that is zeroed must be different from the row that takes the length");

	number a = matrix.at(row_zero, col);
	number b = matrix.at(row_len, col);

	double normalisation = std::sqrt(a.radius_squared() + b.radius_squared());

	number c_factor = b/normalisation; 
	number s_factor = a/normalisation;

	for(int i = 0; i < matrix.cols(); i++){

		number temp = matrix.at(row_zero,i);
		matrix.at(row_zero, i) = test_zero(-s_factor * matrix.at(row_len,i) + c_factor * matrix.at(row_zero,i));
		matrix.at(row_len, i) = c_factor * matrix.at(row_len,i) + s_factor * temp; 

	}	
}


//}

//{ Orthonormalise
Matrix orthonormalise(Matrix matrix)
{
	// uses the Gram-Schmidt method for finding an orthonormal basis of this matrix' column space
	Matrix answer = Matrix::empty();
	
	int i = 0;
	while(is_nullvector(matrix.get_col(i)))
	{
		i++;
		if(i == matrix.cols())
			return matrix;
	}
	answer.push_back_col(matrix.get_col(i)); // finds the first non zero column

	double lngth = length(answer.get_col(0));	
	for(int i = 0; i < answer.rows();i++)
		answer[i][0] = answer[i][0]/lngth;
	

	for(int j = i + 1; j < matrix.cols(); j++)
	{
		Matrix column = Matrix::empty();
		column.push_back_col(matrix.get_col(j)); 
		if(is_nullvector(column.get_col(0))) // skips null columns
			continue;

		for(int p = 0; p < answer.cols(); p++) // loops over the previous vectors
		{
			Matrix temp = Matrix::empty(); 
			temp.push_back_col(answer.get_col(p));
			column = column - multiply(projection(temp), column); // subtracts the part in this direction
		}

		lngth = length(column.get_col(0)); // divides by length
		for(int i = 0; i < column.rows();i++)
			column[i][0] = column[i][0]/lngth;
		answer.append(column);
	}



	return answer;	

}

//}

//{ Decompositions 

//{ QR decomposition

std::vector<number> householder_transform(std::vector<number> seed, std::vector<number> transformed)
{
	return transformed - (2 * (inner_product(seed, transformed)) / inner_product(seed, seed) ) * seed;
}

Matrix householder_matrix(std::vector<number> seed) 
{
	return Matrix::identity(seed.size()) - 2 * (multiply(seed, seed) / inner_product(seed, seed));
}

std::vector<Matrix> QR_decomposeHS(Matrix matrix) // Performs QR decomposition using Householder transformations
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("QR_decompose(): only square matrices can be decomposed this way");

	std::vector<Matrix> QR;
		
	Matrix Q = Matrix::identity(matrix.rows());
	Matrix R = matrix;

	for(int i = 0; i < matrix.cols() - 1; i++){
		
		Matrix H = Matrix::identity(matrix.rows()); 
		
		std::vector<number> difference(matrix.rows() - i,0);

		for(int j = i; j < matrix.rows(); j++){
			
			difference[j - i] = R[j][i];
		}

		double len = length(difference);

		if(difference[0].re() > 0)
			difference[0] = difference[0] + len;
		else 
			difference[0] = difference[0] - len;

		len = length(difference);
		
		if(len != 0){

			for(int p = i; p < matrix.cols(); p++){
				for(int q = i; q < matrix.rows(); q++){

					H[q][p] = H[q][p] - ((2 * (difference[p - i]*difference[q - i])) / (len*len));

				}
			}

		}
		Q = multiply(Q, H); 
		R = multiply(H, R);
			
		
	}

	QR.push_back(Q);
	QR.push_back(R);

	return QR;
}

std::vector<Matrix> QR_decomposeGS(Matrix matrix) // Permorms QR decomposition using the Gram-Schmidt process 
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("QR_decompose(): only square matrices can be decomposed this way");

	std::vector<Matrix> QR;
		
	
	Matrix Q = Matrix::empty();
	Matrix R = Matrix::nullmatrix(matrix.rows(),matrix.cols());

	if(is_nullvector(matrix.get_col(0)))
		throw std::invalid_argument("QR_decompose(): only invertible matrices can be decomposed this way");
	
	Q.push_back_col(matrix.get_col(0));


	double lngth = length(Q.get_col(0));	
	R[0][0] = lngth;
	for(int i = 0; i < Q.rows();i++)
		Q[i][0] = Q[i][0]/lngth;


	for(int j = 1; j < matrix.cols(); j++)
	{
		Matrix column = Matrix::empty();
		column.push_back_col(matrix.get_col(j)); 
		if(is_nullvector(column.get_col(0))) // skips null columns
			throw std::invalid_argument("QR_decompose(): only invertible matrices can be decomposed this way");

		for(int p = 0; p < Q.cols(); p++) // loops over the previous vectors
		{
			Matrix temp = Matrix::empty(); 
			temp.push_back_col(Q.get_col(p));
		 	number dot_product = inner_product(temp.get_col(0),column.get_col(0));
			R[p][j] = dot_product;
			column = column - (dot_product * temp);// subtracts the part in this direction
		}

		lngth = length(column.get_col(0));
		
		R[j][j] = lngth; // sets the diagonal entries

		for(int i = 0; i < column.rows();i++) // divides by length
			column[i][0] = column[i][0]/lngth;
		
		Q.append(column);
	}
	
	QR.push_back(Q);
	QR.push_back(R);
	
	
	return QR;
}
//}

//{ PLU decomposition
std::vector<Matrix> PLU_decompose(Matrix matrix)
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("PLU_decompose(): this operation is defined here only for square matrices");
	std::vector<Matrix> PLU;
	
	Matrix identity = Matrix::identity(matrix.rows());
	Matrix L = identity;
	int col = 0;
	int row = 0;

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		unsigned int pivot_row = find_nonzero_pivot(matrix, row, col);

		if(pivot_row < matrix.rows()) // if none is found, skips this column
		{

			matrix.swap_rows(pivot_row, row); // if there is one, swaps rows
			identity.swap_rows(pivot_row, row);
			L.swap_rows(pivot_row, row);
			
			// Quicker "swap columns" which are standard basis vectors

			L[row][row] = 1;
			L[pivot_row][row] = 0;
			L[row][pivot_row] = 0;
			L[pivot_row][pivot_row] = 1;
			
			number pivot = matrix[row][col];
			for(int i = row + 1; i < matrix.rows(); i++) // for all rows
			{
				number factor = matrix[i][col]/pivot; // finds a factor

				for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
				{
					matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
				}
				L[i][col] = factor;
			}
			row++;
		}

		col++;
	}

	Matrix P = transpose(identity); // the inverse of an orthonormal matrix is its transpose

	PLU.push_back(P);
	PLU.push_back(L);
	PLU.push_back(matrix);
	return PLU;

}
//}

//}

//{ Solve system

std::vector<Matrix> solve_system(Matrix LHS, std::vector<number> RHS)
{
	if(LHS.rows() != RHS.size())
		throw std::invalid_argument("solve_system(): invalid system");

	Matrix NullSpace = nullspace(LHS);
	LHS.push_back_col(-RHS);
	Matrix nullOfAppended = nullspace(LHS);
	Matrix particular_solution = Matrix::empty();
	int i = 0;
	while(nullOfAppended[nullOfAppended.rows() - 1][i] != 1)
	{
		i++;
		if(i == nullOfAppended.cols())
			return {}; 
	}

	std::vector<number> solution_base;
	for(int j = 0; j < nullOfAppended.rows() - 1; j++)
		solution_base.push_back(nullOfAppended[j][i]);

	particular_solution.push_back_col(solution_base);
	
	if(NullSpace.cols() == 1)
		return {Matrix::nullmatrix(particular_solution.rows(), 1), particular_solution};
	else 
		return  {NullSpace, particular_solution};
	
}

//}

//{ Best solution

std::vector<Matrix> best_solve(Matrix LHS, std::vector<number> RHS)
{
	return solve_system(LHS,multiply(projection(LHS),RHS));	

}


//}

//{ Rank

int rank(Matrix matrix)
{
	std::vector<std::pair<unsigned int, unsigned int>> pivot_indeces;
	ref_remember_pivots(matrix, pivot_indeces);
	return pivot_indeces.size();

}
//}

//{ Determinant 
number determinant(Matrix matrix)
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("determinant(): this matrix is not square: determinant not defined");


	number determinant = 1;

	/// HERE PERFORM THE SPECIAL ROW SUBTRACTION TO COMPUTE THE DETERMINANT

	int col = 0;
	int row = 0;

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		unsigned int pivot_row = find_nonzero_pivot(matrix, row, col);

		if(pivot_row < matrix.rows()) // if none is found, skips this column
		{
			if(pivot_row != row){
				matrix.swap_rows(pivot_row, row); // if there is one, swaps rows
				determinant = determinant * (-1);	
			}

			clear_pivot_column(matrix, row, col);

			row++;
		}
		col++;
	}


	for(int i = 0; i < matrix.cols(); i++)
	{
		determinant = determinant * matrix[i][i];
		if(determinant == 0.0)
			return 0;
	}
	return determinant;
}


//}

//{ Eigenvalues UNFINISHED


number power_method(Matrix matrix, int iterations) 
{
	std::vector<number> current(matrix.cols(),1);
	current = current / length(current);
	number current_guess = 0;

	for(int i = 0; i < iterations; i++)
	{
		std::vector<number> next = multiply(matrix, current);
		current_guess = inner_product(next, current);
		current = next / length(next);
	}

	return current_guess;
}



std::vector<number> eigenvalues(Matrix matrix) // utilises the QR algorithm to compute the eigenvalues
{
	transform_into_upper_hessenberg(matrix);
		
		
	


	return {};
}


//}
