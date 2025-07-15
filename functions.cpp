#include "matrix.h"
#include <cmath>
#include <numeric>
#include <stdexcept>

//{ Print functions 
void print(Matrix to_print)
{
	for(int i = 0; i < to_print.rows(); i++){
		for(int j = 0; j < to_print.cols(); j++){
			if(zero(to_print[i][j]))
				std::cout << 0 << " ";
			else
				std::cout << to_print[i][j] << " ";
		}
		std::cout << "\n";
	}
	std::cout<< "\n";
}

void print(std::vector<float> to_print)
{
	for(int i = 0; i < to_print.size(); i++)
		std::cout << to_print[i] << " ";
	std::cout << "\n";
}
//}

//{ Null validation functions
bool zero(float number)
{
	if(number < 0.0001 && number > -0.0001)
		return true;
	else
		return false;

}

bool is_nullvector(std::vector<float> vect)
{
	for(float x : vect )
	{
		if(!zero(x))
			return false;
	}
	return true;
}

float test_zero(float number)
{
	if(zero(number))
		return 0;
	else 
		return number;

}

Matrix clear_zero_floats(Matrix to_clear)
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
std::vector<float> operator+(std::vector<float> one, std::vector<float> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("operator+(): undefined addition of different length vectors");
	std::vector<float> answer(one.size(),0);

	for(int i = 0; i < one.size(); i++)
		answer[i] = one[i] + two[i];
	
	return answer;
}


std::vector<float> operator-(std::vector<float> one, std::vector<float> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("operator+(): undefined addition of different length vectors");
	std::vector<float> answer(one.size(),0);

	for(int i = 0; i < one.size(); i++)
		answer[i] = one[i] - two[i];
	
	return answer;
}

Matrix operator+(Matrix one, Matrix two)
{
	if(one.rows() != two.rows() || one.cols() != two.cols())
		throw std::invalid_argument("operator+(): undefined addition of two different size matrices");
	
	Matrix answer(one.rows(),one.cols());

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
	
	Matrix answer(one.rows(),one.cols());

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

std::vector<float> operator-(std::vector<float> one)
{
	for(int i = 0; i < one.size(); i++)
	{
		one[i] = -one[i];
	}
	return one;
}

Matrix operator-(Matrix left, std::vector<float> right)
{
	if(left.cols() > 1)
		throw(std::invalid_argument("operator-(): operation undefined"));

	for(int i = 0; i < left.rows(); i++)
		left[i][0] = left[i][0] - right[i];
	return left;
}

std::vector<float> operator*(float scalar, std::vector<float> vect)
{
	for(auto x : vect)
		x = x * scalar;
	return vect;
}

Matrix operator*(float scalar,Matrix matrix)
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
//}

//{ length

float length(std::vector<float> vect)
{
	float length = 0;
	for(float x : vect)
	{
		length += x * x;
	}
	length = sqrt(length);
	return length;
}


//}

//{ Inner products

float inner_product(std::vector<float> one, std::vector<float> two)
{
	if(one.size() != two.size())
		throw std::invalid_argument("inner_product(): undefined product of different length vectors");
	
	float answer = 0;
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
	Matrix transposed(original.cols(),original.rows());

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

	Matrix answer(left.rows(),right.cols());

	for(int i = 0; i < answer.rows(); i++){
		for(int j = 0; j < answer.cols(); j++){
			for(int p = 0; p < left.cols(); p++){
				answer[i][j] = test_zero(answer[i][j] + left[i][p]*right[p][j]);
			}
		}
	}
	return answer;
}

std::vector<float> multiply(Matrix left,std::vector<float> right)
{
	if(left.cols() != right.size())
		throw std::invalid_argument("multiply(): invalid vector length");

	std::vector<float> answer(left.rows(), 0);

	for(int i = 0; i < left.rows(); i++){
		for(int j = 0; j < left.cols(); j++){
			answer[i] = test_zero(answer[i] + left[i][j] * right[j]);
		}
	}

	return answer;
}
//}

//{ Row Echelon form
Matrix row_echelon(Matrix matrix)
{
	int col = 0;
	int row = 0;

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		int temp_row = row;
		while(zero(matrix[temp_row][col]) && temp_row < matrix.rows())
		{
			temp_row++;

			if(temp_row == matrix.rows()) // if none is found, skips this column
			{
				col++;
				break;
			}
		}
		matrix.swap_rows(temp_row, row); // if there is one, swaps rows
		
		float pivot = matrix[row][col];
		for(int i = row + 1; i < matrix.rows(); i++) // for all rows
		{
			float factor = matrix[i][col]/pivot; // finds a factor

			for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
			{
				matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
			}
		}

	col++;
	row++;
	}

	return matrix;
}
//}

//{ Reduced row echelon form
Matrix reduced_row_echelon(Matrix matrix)
{	
	int col = 0;
	int row = 0;
	std::vector<std::pair<float,float>> pivot_indeces; // keep track of pivots to clear 

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		int temp_row = row;
		while(zero(matrix[temp_row][col]) && temp_row < matrix.rows())
		{
			temp_row++;

			if(temp_row == matrix.rows()) // if none is found, skips this column
			{
				col++;
				break;
			}
		}
		matrix.swap_rows(temp_row, row); // if there is one, swaps rows
		float pivot = matrix[row][col];
		pivot_indeces.push_back({row,col});

		for(int i = row + 1; i < matrix.rows(); i++) // for all rows
		{
			float factor = matrix[i][col]/pivot; // finds a factor

			for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
			{
				matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
			}
		}

	col++;
	row++;
	}

	// clears the upper rows too
	for(int i = pivot_indeces.size() - 1; i >= 0; i--) 	// for all pivots
	{ 

		float pivot = matrix[pivot_indeces[i].first][pivot_indeces[i].second];

		for(int j = pivot_indeces[i].first - 1; j >= 0; j--) // for all rows above this pivot
		{
			float factor = matrix[j][pivot_indeces[i].second]/pivot;	

			for(int p = matrix.cols() - 1; p >= pivot_indeces[i].second; p--) // for all columns in that row
			{
				matrix[j][p] = test_zero(matrix[j][p] - factor * matrix[pivot_indeces[i].first][p]);
			}
		}
		for(int g = matrix.cols() - 1; g >= pivot_indeces[i].second; g--)
			matrix[pivot_indeces[i].first][g] = matrix[pivot_indeces[i].first][g] / pivot;
	}


	return matrix;
}

//}

//{ Inverse
Matrix inverse(Matrix matrix)
{
	if(matrix.cols() != matrix.rows()) // check if the matrix is square
		throw std::invalid_argument("inverse(): the matrix isn't square");

	Matrix identity(matrix.cols());
	matrix.append(identity);
	matrix = reduced_row_echelon(matrix);

	for(int i = 0; i < matrix.rows(); i++) // check if this matrix is singular
	{
		if(matrix[i][i] == 0)
			throw std::invalid_argument("inverse(): this matrix is singular and has no inverse");
	}

	Matrix answer;
	for(int i = matrix.cols()/2; i < matrix.cols(); i++) // create the answer;
		answer.push_back_col(matrix.get_col(i));

	return answer;
}
//}

//{ Nullspace
Matrix nullspace(Matrix matrix)
{
	// First get to the reduced row echelon form and remember pivots

	int col = 0;
	int row = 0;
	std::vector<std::pair<float,float>> pivot_indeces; // keep track of pivots to clear 

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		int temp_row = row;
		while(zero(matrix[temp_row][col]) && temp_row < matrix.rows())
		{
			temp_row++;

			if(temp_row == matrix.rows()) // if none is found, skips this column
			{
				col++;
				break;
			}
		}
		matrix.swap_rows(temp_row, row); // if there is one, swaps rows
		float pivot = matrix[row][col];
		pivot_indeces.push_back({row,col});

		for(int i = row + 1; i < matrix.rows(); i++) // for all rows
		{
			float factor = matrix[i][col]/pivot; // finds a factor

			for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
			{
				matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
			}
		}

	col++;
	row++;
	}

	// clears the upper rows too
	for(int i = pivot_indeces.size() - 1; i >= 0; i--) 	// for all pivots
	{ 

		float pivot = matrix[pivot_indeces[i].first][pivot_indeces[i].second];

		for(int j = pivot_indeces[i].first - 1; j >= 0; j--) // for all rows above this pivot
		{
			float factor = matrix[j][pivot_indeces[i].second]/pivot;	

			for(int p = matrix.cols() - 1; p >= pivot_indeces[i].second; p--) // for all columns in that row
			{
				matrix[j][p] = test_zero(matrix[j][p] - factor * matrix[pivot_indeces[i].first][p]);
			}
		}
		for(int g = matrix.cols() - 1; g >= pivot_indeces[i].second; g--)
			matrix[pivot_indeces[i].first][g] = matrix[pivot_indeces[i].first][g] / pivot;
	}
	// After getting to rref, finally the nullspace is generated
	Matrix answer;
	std::vector<float> null_vector(matrix.cols(),0);
	answer.push_back_col(null_vector);

	int which_pivot = 0;
	for(int i = 0; i < matrix.cols(); i++) // iterate through, skipping pivot columns
	{
		std::vector<float> partial_answer(matrix.cols(),0);
		while(pivot_indeces[which_pivot].second == i)
		{
			i++;
			which_pivot++;
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

//{ Orthonormalise
Matrix orthonormalise(Matrix matrix)
{
	// uses the Gram-Schmidt method for finding an orthonormal basis of this matrix' column space
	Matrix answer;
	
	int i = 0;
	while(is_nullvector(matrix.get_col(i)))
	{
		i++;
		if(i == matrix.cols())
			return matrix;
	}
	answer.push_back_col(matrix.get_col(i)); // finds the first non zero column

	float lngth = length(answer.get_col(0));	
	for(int i = 0; i < answer.rows();i++)
		answer[i][0] = answer[i][0]/lngth;
	

	for(int j = i + 1; j < matrix.cols(); j++)
	{
		Matrix column;
		column.push_back_col(matrix.get_col(j)); 
		if(is_nullvector(column.get_col(0))) // skips null columns
			continue;

		for(int p = 0; p < answer.cols(); p++) // loops over the previous vectors
		{
			Matrix temp; 
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
std::vector<Matrix> QR_decompose(Matrix matrix)
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("QR_decompose(): only square matrices can be decomposed this way");

	std::vector<Matrix> QR;
		
	
	Matrix Q;
	Matrix R(matrix.rows(),matrix.cols());

	int i = 0;
	while(is_nullvector(matrix.get_col(i)))
	{
		i++;
		if(i == matrix.cols())
			throw std::invalid_argument("QR_decompose(): this matrix cannot be decomposed");
	}
	Q.push_back_col(matrix.get_col(i)); // finds the first non zero column

	float lngth = length(Q.get_col(0));	
	R[0][0] = lngth;
	for(int i = 0; i < Q.rows();i++)
		Q[i][0] = Q[i][0]/lngth;


	for(int j = i + 1; j < matrix.cols(); j++)
	{
		Matrix column;
		column.push_back_col(matrix.get_col(j)); 
		if(is_nullvector(column.get_col(0))) // skips null columns
			continue;

		for(int p = 0; p < Q.cols(); p++) // loops over the previous vectors
		{
			Matrix temp; 
			temp.push_back_col(Q.get_col(p));
		 	float dot_product = inner_product(temp.get_col(0),column.get_col(0));
			R[p][j-i] = dot_product;
			column = column - (dot_product * temp);// subtracts the part in this direction
		}

		lngth = length(column.get_col(0));
		
		R[j - i][j - i] = lngth; // sets the diagonal entries

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
	
	Matrix identity(matrix.rows());
	Matrix L = identity;
	int col = 0;
	int row = 0;

	// finds the next non-zero pivot
	while(col < matrix.cols() && row < matrix.rows())
	{
		int temp_row = row;
		while(zero(matrix[temp_row][col]) && temp_row < matrix.rows())
		{
			temp_row++;
	
			if(temp_row == matrix.rows()) // if none is found, skips this column
			{
				col++;
				break;
			}
		}

		matrix.swap_rows(temp_row, row); // if there is one, swaps rows
		identity.swap_rows(temp_row, row);
		float pivot = matrix[row][col];
		for(int i = row + 1; i < matrix.rows(); i++) // for all rows
		{
			float factor = matrix[i][col]/pivot; // finds a factor
			L[i][row] = factor;

			for(int j = col; j < matrix.cols(); j++) // then for all columns in these rows
			{
				matrix[i][j] = test_zero(matrix[i][j] - matrix[row][j] * factor); // subtracts the pivot row times factor
			}
		}

	col++;
	row++;
	}
	Matrix P = transpose(identity); // the inverse of an orthonormal matrix is its transpose
	Matrix U = matrix.split_left(2); // this is what is left after elimination
	PLU.push_back(P);
	PLU.push_back(L);
	PLU.push_back(U);
	return PLU;

}
//}

//}

//{ Solve system

Matrix solve_system(Matrix LHS, std::vector<float> RHS)
{
	if(LHS.rows() != RHS.size())
		throw std::invalid_argument("solve_system(): invalid system");

	Matrix NullSpace = nullspace(LHS);
	LHS.push_back_col(-RHS);
	Matrix nullOfAppended = nullspace(LHS);
	Matrix particular_solution;
	int i = 0;
	while(nullOfAppended[nullOfAppended.rows() - 1][i] != 1)
	{
		Matrix empty;
		i++;
		if(i == nullOfAppended.cols())
			throw std::invalid_argument("solve_system(): this system has no solution"); 
	}
	std::vector<float> solution_base;
	for(int j = 0; j < nullOfAppended.rows() - 1; j++)
		solution_base.push_back(nullOfAppended[j][i]);

	particular_solution.push_back_col(solution_base);
	
	NullSpace.append(particular_solution);
	
	return NullSpace.split_right(1);
}

//}

//{ Best solution

Matrix best_solve(Matrix LHS, std::vector<float> RHS)
{
	return solve_system(LHS,multiply(projection(LHS),RHS));	

}


//}

//{ Rank
int rank(Matrix matrix)
{
	matrix = row_echelon(matrix);

	int max_rank = std::min(matrix.rows(),matrix.cols());
	for(int i = 0; i < max_rank; i++)
	{
		if(matrix[i][i] == 0)
			max_rank--;
	}

	return max_rank;

}
//}

//{ Determinant CORRECT MISTAKE SIGN
float determinant(Matrix matrix)
{
	if(matrix.rows() != matrix.cols())
		throw std::invalid_argument("determinant(): this matrix is not square: determinant not defined");
	
	float determinant = 1;
	for(int i = 0; i < matrix.cols(); i++)
	{
		determinant = determinant * matrix[i][i];
		if(zero(determinant))
			return 0;
	}
	return determinant;
}


//}

//{ Pivots
std::vector<float> pivots(Matrix matrix)
{
	std::vector<float> pivots;
	matrix = row_echelon(matrix);
	for(int i = 0; i < std::min(matrix.rows(),matrix.cols());i++)
	{
		if(!zero(matrix[i][i]))
			pivots.push_back(matrix[i][i]);
	}

	return pivots;
}
//}

//{ Eigenvalues UNFINISHED
std::vector<float> eigenvalues(Matrix matrix) // utilises the QR algorithm to compute the eigenvalues
{
	Matrix A = matrix;
	std::vector<Matrix> QR = QR_decompose(A);
	
	
		



	return {};
}


//}
