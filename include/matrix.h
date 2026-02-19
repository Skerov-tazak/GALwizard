#pragma once
#include<vector>
#include"../src/complex/number.h"
#include<iostream>

namespace gal {

	class Matrix
	{
		private:
			
			class ConstRowProxy 
			{
				private:

					const Matrix* matrix;

					unsigned int row;

				public:

					ConstRowProxy(const Matrix& m, unsigned int r);

					const number& operator[](unsigned int col);
			};

			class RowProxy 
			{
				private:

					Matrix* matrix;

					unsigned int row;

				public:

					RowProxy(Matrix& m, unsigned int r);

					number& operator[](unsigned int col);
			};

			std::vector<number> data;

			unsigned int col_len;
			
			Matrix(unsigned int); // Identity Matrix constructor
		
			Matrix(unsigned int,unsigned int); // Nullmatrix of m by n constructor
		
			Matrix(); // creates an empty matrix

			unsigned int real_index(unsigned int, unsigned int)const; // Calculates i,j row/col into flat 1D index
			 
		public:

			friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

			Matrix(const std::vector<std::vector<number>>&); // constructs a full Matrix (ROW BY ROW)
			
			bool is_empty(); // returns whether it is an empty matrix

			int rows()const; // Returns number of rows
			
			int cols()const; // Returns number of columns
		
			void push_back_col(const std::vector<number>&); // adds a column at the back
			
			void swap_rows(unsigned int,unsigned int); // swaps two rows

			RowProxy operator[](unsigned int); // Operator overload for getting data
			
			ConstRowProxy operator[](unsigned int)const; // Operator overload for getting data

			number& at(unsigned int, unsigned int); // returns refernce to value at index

			const number& at(unsigned int, unsigned int)const;
			
			void append(const Matrix&); // Appends another matrix on the right of this one

			Matrix split_right(unsigned int); // returns submatrix to the right of this index (Included)

			Matrix split_left(unsigned int); // returns submatrix to the left of this index (Included);

			static Matrix identity(unsigned int); // returns an n by n identity matrix
										  
			static Matrix nullmatrix(unsigned int, unsigned int); // return an n by m nullmatrix 

			static Matrix empty();	// returns an empty matrix;

			static Matrix column(const std::vector<number>&);

			std::vector<number> get_row(unsigned int)const; // returns a row in vector form

			std::vector<number> get_col(unsigned int)const; // returns a column in vector form
	};
}
	
