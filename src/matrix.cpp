#include "matrix.h"
#include <stdexcept>

namespace gal {


	Matrix::Matrix() 
	{
		data.resize(0);
		col_len = 0;
	}

	Matrix::Matrix(unsigned int rows, unsigned int cols)
	{
		col_len = rows;
		data.resize(col_len * cols, 0);
	}

	Matrix::Matrix(unsigned int size)
	{
		col_len = size;
		data.resize(size * size);
		for(int i = 0; i < size; i++)
			data.at(real_index(i,i)) = 1;
	}


	Matrix::Matrix(const std::vector<std::vector<number>>& input)
	{ 
		col_len = input.size();
		if(input.size() != 0)
			data.resize(input.size() * input[0].size());
		
		for(unsigned int i = 0; i < input.size(); i++){
			if (input[i].size() != input[0].size())
				throw std::invalid_argument("Matrix(vector<vector<number>): ununiform column length");
			
			for(int j = 0; j < input[0].size(); j++)
			data[real_index(i, j)] = input[i][j];

		}
	}

	Matrix Matrix::column(const std::vector<number>& col) {

			Matrix created;
			created.data = col;
			created.col_len = col.size();

			return created;
	
	}

	Matrix Matrix::identity(unsigned int size)
	{
		return Matrix(size);
	}

	Matrix Matrix::nullmatrix(unsigned int m, unsigned int n)
	{
		return Matrix(m,n);
	}

	bool Matrix::is_empty(){
		
		if(data.size() == 0)
			return true;
		else 
			return false;
	}

	int Matrix::rows()const 
	{
		return col_len;
	}

	int Matrix::cols()const
	{
		if(rows() == 0)
			return 0;
		return data.size()/col_len;
	}

	void Matrix::push_back_col(const std::vector<number>& new_col)
	{
		if(rows() == 0)
		{
			data = new_col;
			col_len = new_col.size();
		}
		else
		{ 
			if(new_col.size() != rows())
				throw std::invalid_argument("push_back_col(): invalid column length");

			data.reserve(data.size() + col_len);
			data.insert(data.end(), new_col.begin(), new_col.end());
		}
	}

	void Matrix::append(const Matrix& to_append)
	{
		if(to_append.rows() != rows())
			throw std::invalid_argument("append(): invalid size matrix");

		for(int i = 0; i < to_append.cols(); i++)
			push_back_col(to_append.get_col(i));
	}

	unsigned int Matrix::real_index(unsigned int row, unsigned int col)const{
		return col * col_len + row;
	}

	Matrix::ConstRowProxy::ConstRowProxy(const Matrix& m, unsigned int r){
		matrix = &m;
		row = r;
	}

	const number& Matrix::ConstRowProxy::operator[](unsigned int col){

		return matrix->data[col * matrix->col_len + row];
	}

	Matrix::RowProxy::RowProxy(Matrix& m, unsigned int r){
		matrix = &m;
		row = r;
	}

	number& Matrix::RowProxy::operator[](unsigned int col){

		return matrix->data[col * matrix->col_len + row];
	}

	Matrix::RowProxy Matrix::operator[](unsigned int index)
	{
		return RowProxy(*this, index);
	}

	Matrix::ConstRowProxy Matrix::operator[](unsigned int index)const
	{
		return ConstRowProxy(*this, index);
	}

	number& Matrix::at(unsigned int row, unsigned int col){

		return data.at(real_index(row, col));
	}

	const number& Matrix::at(unsigned int row, unsigned int col)const{

		return data.at(real_index(row, col));
	}

	Matrix Matrix::split_right(unsigned int index)
	{

		Matrix answer;
		if(index > cols())
			throw std::invalid_argument("split_right(): invalid index");
		
		for(int i = index; i < cols(); i++)
		{
			answer.push_back_col(get_col(i));
		}
		return answer;
	}

	Matrix Matrix::split_left(unsigned int index)
	{
		if(index > cols())
			throw std::invalid_argument("split_left(): invalid index");

		Matrix answer;
		for(unsigned int i = 0; i < index + 1; i++)
		{
			answer.push_back_col(get_col(i));
		}
		return answer;

	}


	Matrix Matrix::empty()
	{
		return Matrix();
	}

	std::vector<number> Matrix::get_row(unsigned int index)const
	{
		std::vector<number> answer;
		answer.reserve(cols());

		if(index >= rows())
			throw std::invalid_argument("get_row(): invalid row index");
		
		for(unsigned int i = 0; i < cols(); i++){
			answer.push_back(data[real_index(index, i)]);
		}

		return answer;
	}

	std::vector<number> Matrix::get_col(unsigned int index)const
	{
		if(index >= cols())
			throw std::invalid_argument("get_col(): invalid col index");

		std::vector<number> answer;
		answer.reserve(rows());
		unsigned int start = index * col_len;

		for(int i = start; i < rows() + start; i++)
			answer.push_back(data[i]);

		return answer;
	}

	void Matrix::swap_rows(unsigned int first, unsigned int second)
	{
		if(first == second)
			return;

		if(first >= rows() || second >= rows())
			throw std::invalid_argument("swap_rows(): invalid row indeces");


		for(unsigned int i = 0; i < cols(); i++){
			number temp = data.at(real_index(first, i));
			data.at(real_index(first, i)) = data.at(real_index(second, i));
			data.at(real_index(second, i)) = temp;
		}
	}

}
