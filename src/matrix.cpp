#include "../include/matrix.h"
#include <initializer_list>
#include <stdexcept>

Matrix::Matrix() 
{

}

Matrix::Matrix(unsigned int rows, unsigned int cols)
{
	std::vector<float> dummy(cols,0);
	data.resize(rows,dummy);
}

Matrix::Matrix(unsigned int size)
{
	std::vector<float> dummy(size,0);
	data.resize(size,dummy);
	for(int i = 0; i < size; i++)
		data[i][i] = 1;
}


Matrix::Matrix(std::vector<std::vector<float>> input)
{ 
	data = input;
}

int Matrix::rows()const 
{
	return data.size();
}

int Matrix::cols()const
{
	if(rows() == 0)
		return 0;
	return data[0].size();
}

void Matrix::push_back_col(std::vector<float> new_col)
{
	if(rows() == 0)
	{
		for(int i = 0; i < new_col.size(); i++)
			data.push_back({new_col[i]});
	}
	else
	{ 
		if(new_col.size() != rows())
			throw std::invalid_argument("push_back_col(): invalid column length");

		for(int i = 0; i < rows(); i++)
		{
			data[i].push_back(new_col[i]);
		}
	
	}
}

void Matrix::push_back_row(std::vector<float> new_row)
{
	if(cols() != 0 && new_row.size() != cols())
		throw std::invalid_argument("push_back_row(): invalid row length");

	data.push_back(new_row);
}

void Matrix::append(Matrix to_append)
{
	if(to_append.rows() != rows())
		throw std::invalid_argument("append(): invalid size matrix");

	for(int i = 0; i < to_append.cols(); i++)
		push_back_col(to_append.get_col(i));
}

std::vector<float>& Matrix::operator[](unsigned int index)
{
	return data[index];	
}

float Matrix::at(unsigned int row, unsigned int col)const{

		
	return data.at(row).at(col);
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
	for(int i = 0; i < index + 1; i++)
	{
		answer.push_back_col(get_col(i));
	}
	return answer;

}

Matrix Matrix::identity(unsigned int size)
{
	return Matrix(size);
}

Matrix Matrix::nullmatrix(unsigned int m, unsigned int n)
{
	return Matrix(m,n);
}

std::vector<float> Matrix::get_row(unsigned int index)const
{
	if(index >= rows())
		throw std::invalid_argument("get_row(): invalid row index");
	return data[index];
}

std::vector<float> Matrix::get_col(unsigned int index)const
{
	if(index >= cols())
		throw std::invalid_argument("get_col(): invalid col index");
	std::vector<float> answer(rows(),0);
	for(int i = 0; i < rows(); i++)
		answer[i] = data[i][index];

	return answer;
}

void Matrix::swap_rows(unsigned int first, unsigned int second)
{
	if(first == second)
		return;

	if(first > rows() || second > rows())
		throw std::invalid_argument("swap_rows(): invalid row indeces");

	std::vector<float> temp = data[first];
	data[first] = data[second];
	data[second] = temp;
}

