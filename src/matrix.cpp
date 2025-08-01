#include "../include/matrix.h"
#include <initializer_list>
#include <stdexcept>

Matrix::Matrix() 
{
	data.resize(0);
}

Matrix::Matrix(unsigned int rows, unsigned int cols)
{
	std::vector<number> dummy(cols,0);
	data.resize(rows,dummy);
}

Matrix::Matrix(unsigned int size)
{
	std::vector<number> dummy(size,0);
	data.resize(size,dummy);
	for(int i = 0; i < size; i++)
		data[i][i] = 1;
}


Matrix::Matrix(const std::vector<std::vector<number>>& input)
{ 
	data = input;
}

Matrix Matrix::column(const std::vector<number>& col) {

		Matrix created;
        unsigned int n = col.size();
        
		created.data = std::vector<std::vector<number>>(n, std::vector<number>(1));

        for (int i = 0; i < n; ++i)
            created.data[i][0] = col[i];

		return created;
}

bool Matrix::is_empty(){
	
	if(data.size() == 0)
		return true;
	else 
		if(data[0].size() == 0)
			return true;

	return false;
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

void Matrix::push_back_col(std::vector<number> new_col)
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

void Matrix::push_back_row(std::vector<number> new_row)
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

std::vector<number>& Matrix::operator[](unsigned int index)
{
	return data[index];	
}

number Matrix::at(unsigned int row, unsigned int col)const{

		
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

Matrix Matrix::empty()
{
	return Matrix();
}

std::vector<number> Matrix::get_row(unsigned int index)const
{
	if(index >= rows())
		throw std::invalid_argument("get_row(): invalid row index");
	return data[index];
}

std::vector<number> Matrix::get_col(unsigned int index)const
{
	if(index >= cols())
		throw std::invalid_argument("get_col(): invalid col index");
	std::vector<number> answer(rows(),0);
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

	std::swap(data[first],data[second]);
}

