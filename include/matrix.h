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
	
		Matrix(); // creates an empty matrix
				  
	public:
	
		Matrix(std::vector<std::vector<float>>); // constructs a full Matrix (ROW BY ROW)
		
		bool is_empty(); // returns whether it is an empty matrix

		int rows()const; // Returns number of rows
		
		int cols()const; // Returns number of columns
	
		void push_back_col(std::vector<float>); // adds a column at the back
		
		void push_back_row(std::vector<float>); // adds a row at the back
		
		void swap_rows(unsigned int,unsigned int); // swaps two rows

		std::vector<float>& operator[](unsigned int); // Operator overload for getting data

		float at(unsigned int, unsigned int)const; // returns value at index

		void append(Matrix); // Appends another matrix on the right of this one

		Matrix split_right(unsigned int); // returns submatrix to the right of this index (Included)

		Matrix split_left(unsigned int); // returns submatrix to the left of this index (Included);

		static Matrix identity(unsigned int); // returns an n by n identity matrix
									  
		static Matrix nullmatrix(unsigned int, unsigned int); // return an n by m nullmatrix 

		static Matrix empty();	// returns an empty matrix;

		std::vector<float> get_row(unsigned int)const; // returns a row in vector form

		std::vector<float> get_col(unsigned int)const; // returns a column in vector form
};
//}

