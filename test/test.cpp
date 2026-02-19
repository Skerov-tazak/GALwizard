#include "functions.h"
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <vector>

using namespace gal;

//{ Helper Functions

auto check_upper_triangular = [](const Matrix& R) {
	for (int i = 1; i < R.rows(); ++i) {
		for (int j = 0; j < i; ++j) {
			REQUIRE(R.at(i,j) == 0);
		}
	}
};

auto check_lower_triangular = [](const Matrix& L) {
	for (int i = 0; i < L.rows(); ++i) {
		for (int j = i+1; j < L.cols(); ++j) {
			REQUIRE(L.at(i,j) == 0);
		}
	}
};

//}

//{ Arithmetic Tests 

TEST_CASE("Overloading Equals"){

	SECTION(" Square Tests ") {

		Matrix test({{1,2,3},{4,3,2},{5,2,1}});

		SECTION("Testing at"){

			REQUIRE(test.at(2,1) == test[2][1]);
			REQUIRE(test.at(0,0) == test[0][0]);
			REQUIRE(test.at(1,1) == test[1][1]);
			REQUIRE(test.at(0,2) == test[0][2]);
		}
		SECTION( "Testing copy"){

			Matrix test2(test);
			REQUIRE(test2 == test);

		}
		SECTION( "Testing equality" ){

			Matrix test2({{1,2,3},{4,3,2},{5,2,1}});
			REQUIRE(test2 == test);
			test2[0][2] = 9;
			REQUIRE(test2 != test);
		}

	}
	
	SECTION("Rectangular Tests"){
	
		Matrix rect({{3,2},{1,1},{4,3}});
		
		SECTION("copy"){
			
			Matrix second(rect);
			REQUIRE(second == rect);
		}
		SECTION("Null"){

			Matrix nullie({{0,0},{0,0},{0,0}});
			REQUIRE(nullie != rect);
		}
	}
		
}

TEST_CASE("Basic Arithmetic", "[Matrix]") {

	SECTION("Matrix Addition") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 2;
		A[1][0] = 3; A[1][1] = 4;

		Matrix B = Matrix::nullmatrix(2, 2);
		B[0][0] = 5; B[0][1] = 6;
		B[1][0] = 7; B[1][1] = 8;

		Matrix expected = Matrix::nullmatrix(2, 2);
		expected[0][0] = 6; expected[0][1] = 8;
		expected[1][0] = 10; expected[1][1] = 12;

		REQUIRE(A + B == expected);
	}

	SECTION("Matrix Subtraction") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 10; A[0][1] = 9;
		A[1][0] = 8;  A[1][1] = 7;

		Matrix B = Matrix::nullmatrix(2, 2);
		B[0][0] = 1; B[0][1] = 2;
		B[1][0] = 3; B[1][1] = 4;

		Matrix expected = Matrix::nullmatrix(2, 2);
		expected[0][0] = 9; expected[0][1] = 7;
		expected[1][0] = 5; expected[1][1] = 3;

		REQUIRE(A - B == expected);
	}

	SECTION("Matrix Multiplication") {

		SECTION("By identity") {
			Matrix A = Matrix::nullmatrix(2, 2);
			A[0][0] = 3; A[0][1] = 4;
			A[1][0] = 5; A[1][1] = 6;

			Matrix I = Matrix::identity(2);
			Matrix result = multiply(A, I);

			REQUIRE(result == A);
		}

		SECTION("By zero vector") {
			Matrix A = Matrix::nullmatrix(2, 3);
			A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
			A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;

			std::vector<number> v = {0, 0, 0};
			std::vector<number> expected = {0, 0};

			std::vector<number> result = multiply(A, v);

			REQUIRE(result == expected);
		}

		SECTION("By random example") {
			Matrix A = Matrix::nullmatrix(2, 3);
			A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
			A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;

			Matrix B = Matrix::nullmatrix(3, 2);
			B[0][0] = 7;  B[0][1] = 8;
			B[1][0] = 9;  B[1][1] = 10;
			B[2][0] = 11; B[2][1] = 12;

			Matrix expected = Matrix::nullmatrix(2, 2);
			expected[0][0] = 58; expected[0][1] = 64;
			expected[1][0] = 139; expected[1][1] = 154;

			Matrix result = multiply(A, B);

			REQUIRE(result == expected);
		}
	}

	SECTION("Operations on Scalars") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = -2;
		A[1][0] = 3; A[1][1] = -4;

		number scalar = 2.0f;

		Matrix expected = Matrix::nullmatrix(2, 2);
		expected[0][0] = 2; expected[0][1] = -4;
		expected[1][0] = 6; expected[1][1] = -8;

		REQUIRE(A * scalar == expected);
		REQUIRE(scalar * A == expected);
		REQUIRE(expected / scalar == A);
	}
}

//}

//{ Equality of Spaces Tests 

TEST_CASE("Equality of spaces"){
	

	SECTION("Colinearity of Vectors"){

		std::vector<number> one = {1,0};
		std::vector<number> two = {-1,0};

		REQUIRE(vectors_colinear(one,two));

		two[1] = 1;

		REQUIRE(!vectors_colinear(one, two));

		one[1] = -1;

		REQUIRE(vectors_colinear(one,two));

	}
	SECTION("Equality of spans"){
		
		Matrix test1({{2,0,0},{2,1,0},{3,0,1}});
		Matrix test2 = Matrix::identity(3); 

		REQUIRE( test1 != test2 );
		REQUIRE( spaces_equal(test1, test2));

		test1[2][2] = 0;

		REQUIRE( !spaces_equal(test1, test2));

		Matrix zero = Matrix::nullmatrix(3, 3);

		REQUIRE( !spaces_equal(test1, zero) );
		test2[1][1] = 0;
		test2[2][2] = 0;

		REQUIRE( !spaces_equal(test2, zero) );

		zero[0][0] = 1;

		REQUIRE( spaces_equal(zero, test2) );	


	}
	SECTION("Matrices and Vectors"){

		Matrix test1({{1,-1,0},{0,0,0},{0,0,0}});
		Matrix test2({{-1,1},{0,0},{0,0}});
		
		REQUIRE(spaces_equal(test1, test2));

			
	}
}



//}

//{ Gaussian Elimination Tests

TEST_CASE("row_echelon", "[Matrix]") {

	SECTION("Basic 2x2 matrix") {
	
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 2;
		A[1][0] = 3; A[1][1] = 4;

		Matrix expected = Matrix::nullmatrix(2, 2);
		expected[0][0] = 1; expected[0][1] = 2;
		expected[1][0] = 0; expected[1][1] = -2;

		REQUIRE(row_echelon(A) == expected);
	}

	SECTION("Upper triangular matrix remains unchanged") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
		A[1][1] = 4; A[1][2] = 5;
		A[2][2] = 6;

		REQUIRE(row_echelon(A) == A);
	}
}

TEST_CASE("reduced_row_echelon", "[Matrix]") {
	
	SECTION("Basic 2x2 matrix reduces to identity") {
	
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 2;
		A[1][0] = 3; A[1][1] = 4;

		Matrix expected = Matrix::identity(2);

		REQUIRE(reduced_row_echelon(A) == expected);
	}

	SECTION("Identity matrix remains unchanged") {
		Matrix I = Matrix::identity(3);
		REQUIRE(reduced_row_echelon(I) == I);
	}
}

TEST_CASE("inverse", "[Matrix]") {
	
	SECTION("Inverse of 2x2 matrix") {
	
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 4; A[0][1] = 7;
		A[1][0] = 2; A[1][1] = 6;

		Matrix expected = Matrix::nullmatrix(2, 2);
		expected[0][0] = 0.6f;	expected[0][1] = -0.7f;
		expected[1][0] = -0.2f; expected[1][1] = 0.4f;

		REQUIRE(inverse(A) == expected);
	}

	SECTION("Inverse of identity is identity") {
		
		Matrix I = Matrix::identity(3);
		REQUIRE(inverse(I) == I);
	}

	SECTION("Multiplying matrix by its inverse yields identity") {
		
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 3; A[0][1] = 5;
		A[1][0] = 1; A[1][1] = 2;

		Matrix A_inv = inverse(A);
		Matrix I_expected = Matrix::identity(2);

		Matrix product = Matrix::nullmatrix(2, 2);
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 2; ++k)
					product[i][j] += A[i][k] * A_inv[k][j];

		REQUIRE(product == I_expected);
	}
}

//}

//{ Nullspace Tests
TEST_CASE("Nullspace and singular linear systems with number type", "[Matrix][nullspace][solve]") {

	SECTION("Simple singular matrix with nullspace (2x2)") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = number(1); A[0][1] = number(2);
		A[1][0] = number(2); A[1][1] = number(4); // rank 1

		Matrix ns = nullspace(A);

		// Nullspace should include a null vector in column 0, and 1 valid vector in column 1
		REQUIRE(ns.cols() == 2); 
		REQUIRE(ns[0][0] == number(0));
		REQUIRE(ns[1][0] == number(0));

		// Check that column 1 is indeed a valid nullspace vector
		Matrix null_vec = Matrix::nullmatrix(2, 1);
		null_vec[0][0] = ns[0][1];
		null_vec[1][0] = ns[1][1];

		REQUIRE(multiply(A, null_vec) == Matrix::nullmatrix(2, 1));

		// Check orthogonality with transpose (future Hermitian)
		Matrix ns_check = multiply(transpose(null_vec), null_vec);
		REQUIRE(ns_check[0][0] != number(0));
	}

	SECTION("3x3 matrix with imaginary nullspace") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = number(1);  A[0][1] = number(2);	A[0][2] = i();
		A[1][0] = number(2);  A[1][1] = number(4);	A[1][2] = 2*i();
		A[2][0] = number(-1); A[2][1] = number(-2); A[2][2] = -i();

		Matrix ns = nullspace(A);

		// First column should be null vector
		REQUIRE(ns.cols() >= 2);
		for (int r = 0; r < ns.rows(); ++r)
			REQUIRE(ns[r][0] == number(0));

		// Validate remaining nullspace vectors
		for (int c = 1; c < ns.cols(); ++c) {
			Matrix vec = Matrix::nullmatrix(3, 1);
			for (int r = 0; r < 3; ++r) vec[r][0] = ns[r][c];

			REQUIRE(multiply(A, vec) == Matrix::nullmatrix(3, 1));
		}

		// Check Gram matrix of nullspace vectors using transpose
		Matrix ns_gram = multiply(transpose(ns), ns);
		REQUIRE(ns_gram.rows() == ns.cols());
		REQUIRE(ns_gram.cols() == ns.cols());
	}
}

TEST_CASE("Nullspace", "[Matrix][nullspace]") {

	SECTION("Full-rank square matrix so only zero vector in nullspace") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 2;
		A[1][0] = 3; A[1][1] = 4;

		Matrix ns = nullspace(A);

		// Should return 2x1 matrix, all zero
		Matrix expected = Matrix::nullmatrix(2, 1);
		REQUIRE(spaces_equal(ns, expected));
	}

	SECTION("Matrix with a 1D nullspace") {
		// A*x = 0 has solution x = (1, -1)
		Matrix A = Matrix::nullmatrix(1, 2);
		A[0][0] = 1; A[0][1] = 1;

		Matrix ns = nullspace(A);

		// Expected: nullspace has basis vectors (0,0) and (1,-1)
		Matrix expected = Matrix::nullmatrix(2, 1);
		expected[0][0] = 1; expected[1][0] = -1; // basis vector

		REQUIRE(spaces_equal(expected, ns));
	}

	SECTION("Matrix with 2D nullspace (wide zero matrix)") {
		Matrix A = Matrix::nullmatrix(1, 3); // 1 row, 3 columns → rank 0
		// A = [0 0 0], so nullspace is entire ℝ³

		Matrix ns = nullspace(A);

		// Nullspace basis: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
		Matrix expected = Matrix::nullmatrix(3, 4);
		expected[0][0] = 0; expected[1][0] = 0; expected[2][0] = 0; // zero vector
		expected[0][1] = 1; expected[1][1] = 0; expected[2][1] = 0;
		expected[0][2] = 0; expected[1][2] = 1; expected[2][2] = 0;
		expected[0][3] = 0; expected[1][3] = 0; expected[2][3] = 1;

		REQUIRE(spaces_equal(ns, expected));
	}
	SECTION("Zero matrix (e.g. 2x3 all zeros) → full ℝ³ nullspace") {
		Matrix A = Matrix::nullmatrix(2, 3); // All zero

		Matrix ns = nullspace(A);

		// Nullspace basis: zero + 3 standard vectors in ℝ³
		Matrix expected = Matrix::nullmatrix(3, 4);
		expected[0][0] = 0; expected[1][0] = 0; expected[2][0] = 0;
		expected[0][1] = 1; expected[1][1] = 0; expected[2][1] = 0;
		expected[0][2] = 0; expected[1][2] = 1; expected[2][2] = 0;
		expected[0][3] = 0; expected[1][3] = 0; expected[2][3] = 1;

		REQUIRE(spaces_equal(expected, ns));
	}
}

//}

//{ Gram-Schmidt Tests 

TEST_CASE("Orthonormalisation via Gram-Schmidt", "[Matrix][orthonormalise]") {
	SECTION("Standard basis remains unchanged") {
		Matrix A = Matrix::identity(3);
		Matrix result = orthonormalise(A);			 // Returns a vector with one matrix
		REQUIRE(result == A);
	}

	SECTION("Linearly independent non-orthonormal columns") {
		Matrix A = Matrix::nullmatrix(3, 2);
		A[0][0] = 1; A[1][0] = 1; A[2][0] = 0;
		A[0][1] = 1; A[1][1] = 0; A[2][1] = 1;

		Matrix Q = orthonormalise(A);

		Matrix Qt = transpose(Q);
		Matrix I_expected = Matrix::identity(2);
		REQUIRE(multiply(Qt, Q) == I_expected);
	}
}

//}

//{ Decomposition Tests

//{ QR

TEST_CASE("QR decomposition via Gram-Schmidt", "[Matrix][QR_decomposeGS]") {
	SECTION("2x2 matrix") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 1;
		A[1][0] = 1; A[1][1] = -1;

		std::vector<Matrix> qr = QR_decomposeGS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(Q, R) == A);
		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(2));
	}
}

TEST_CASE("QR decomposition via Householder ", "[Matrix][QR_decomposeHS]") {
	SECTION("2x2 identity matrix") {	
		Matrix A = Matrix::identity(2);

		std::vector<Matrix> qr = QR_decomposeHS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(2)); // Q orthogonal
		REQUIRE(multiply(Q, R) == A);							  // Q*R = A
	}

	SECTION("2x2 simple matrix") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = 1; A[0][1] = 1;
		A[1][0] = 1; A[1][1] = -1;

		std::vector<Matrix> qr = QR_decomposeHS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(2));
		REQUIRE(multiply(Q, R) == A);

		// Check R is upper triangular
		REQUIRE(R[1][0] == 0);
	}

	SECTION("3x3 example matrix") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 12;  A[0][1] = -51; A[0][2] = 4;
		A[1][0] = 6;   A[1][1] = 167; A[1][2] = -68;
		A[2][0] = -4;  A[2][1] = 24;  A[2][2] = -41;

		std::vector<Matrix> qr = QR_decomposeHS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);

		// R should be upper triangular
		REQUIRE(R[1][0] == 0);
		REQUIRE(R[2][0] == 0);
		REQUIRE(R[2][1] == 0);
	}
/*	SECTION("3x2 matrix") {
		Matrix A = Matrix::nullmatrix(3, 2);
		A[0][0] = 12; A[1][0] = -51; A[2][0] = 4;
		A[0][1] = 6;  A[1][1] = 167; A[2][1] = -68;

		std::vector<Matrix> qr = QR_decomposeHS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);
	}
	*/
}

TEST_CASE("QR decomposition via Householder (square matrices) - Edge Cases", "[Matrix][QR_decomposeHS]") {


	SECTION("Identity matrix") {
		Matrix A = Matrix::identity(4);
		std::vector<Matrix> qr = QR_decomposeHS(A);
		REQUIRE(qr.size() == 2);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(4));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Zero matrix") {
		Matrix A = Matrix::nullmatrix(3, 3);
		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Upper triangular matrix") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 2; A[0][1] = 3; A[0][2] = 1;
		A[1][1] = 4; A[1][2] = 5;
		A[2][2] = 6;

		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Lower triangular matrix") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 2;
		A[1][0] = 3; A[1][1] = 4;
		A[2][0] = 1; A[2][1] = 5; A[2][2] = 6;

		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Singular matrix (rank deficient)") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
		A[1][0] = 2; A[1][1] = 4; A[1][2] = 6; // Row 2 = 2×Row 1
		A[2][0] = 3; A[2][1] = 6; A[2][2] = 9; // Row 3 = 3×Row 1

		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Negative and mixed values") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = -1; A[0][1] = 4;
		A[1][0] = 2;  A[1][1] = -3;

		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(2));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}

	SECTION("Larger 5x5 matrix (stress test)") {
		Matrix A = Matrix::nullmatrix(5, 5);
		int value = 1;
		for (int i = 0; i < 5; ++i) {
			for (int j = 0; j < 5; ++j) {
				A[i][j] = value++;
			}
		}

		std::vector<Matrix> qr = QR_decomposeHS(A);

		Matrix Q = qr[0], R = qr[1];

		REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(5));
		REQUIRE(multiply(Q, R) == A);
		check_upper_triangular(R);
	}
}

//}

//{ PLU

TEST_CASE("PLU decomposition", "[Matrix][PLU_decompose]") {
	SECTION("3x3 example matrix") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 2; A[0][1] = 0; A[0][2] = 2;
		A[1][0] = 1; A[1][1] = 2; A[1][2] = 3;
		A[2][0] = 3; A[2][1] = 1; A[2][2] = 4;

		std::vector<Matrix> plu = PLU_decompose(A);
		REQUIRE(plu.size() == 3);

		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(multiply(transpose(P), A) == multiply(L, U));
	}

	SECTION("Identity matrix stays the same") {
		Matrix I = Matrix::identity(3);
		std::vector<Matrix> plu = PLU_decompose(I);

		Matrix P = plu[0], L = plu[1], U = plu[2];
		REQUIRE(P == I);
		REQUIRE(L == I);
		REQUIRE(U == I);
	}

	SECTION("Identity matrix") {
		Matrix I = Matrix::identity(3);
		std::vector<Matrix> plu = PLU_decompose(I);

		REQUIRE(plu.size() == 3);
		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(P == I);
		REQUIRE(L == I);
		REQUIRE(U == I);
	}

	SECTION("Simple 3x3 matrix without pivoting") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 2; A[0][1] = 1; A[0][2] = 1;
		A[1][0] = 4; A[1][1] = 3; A[1][2] = 3;
		A[2][0] = 8; A[2][1] = 7; A[2][2] = 9;

		std::vector<Matrix> plu = PLU_decompose(A);
		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(multiply(transpose(P), A) == multiply(L, U));
		check_lower_triangular(L);
		check_upper_triangular(U);
	}

	SECTION("Matrix that requires row swaps (pivoting)") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 0; A[0][1] = 2; A[0][2] = 1;
		A[1][0] = 1; A[1][1] = 0; A[1][2] = 2;
		A[2][0] = 2; A[2][1] = 1; A[2][2] = 0;

		std::vector<Matrix> plu = PLU_decompose(A);
		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(multiply(transpose(P), A) == multiply(L, U));
		check_lower_triangular(L);
		check_upper_triangular(U);

		// Ensure permutation matrix P actually permutes rows
		REQUIRE(P != Matrix::identity(3));
	}

	SECTION("Matrix with multiple pivot swaps") {
		Matrix A = Matrix::nullmatrix(4, 4);
		A[0][0] = 0; A[0][1] = 2; A[0][2] = 1; A[0][3] = 3;
		A[1][0] = 0; A[1][1] = 0; A[1][2] = 4; A[1][3] = 2;
		A[2][0] = 5; A[2][1] = 1; A[2][2] = 0; A[2][3] = 1;
		A[3][0] = 1; A[3][1] = 3; A[3][2] = 2; A[3][3] = 0;

		std::vector<Matrix> plu = PLU_decompose(A);
		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(multiply(transpose(P), A) == multiply(L, U));
		check_lower_triangular(L);
		check_upper_triangular(U);

		// P should definitely not be identity for this one
		REQUIRE(P != Matrix::identity(4));
	}

	SECTION("Singular matrix (pivoting still valid)") {
		Matrix A = Matrix::nullmatrix(3, 3);
		A[0][0] = 2; A[0][1] = 4; A[0][2] = 6;
		A[1][0] = 1; A[1][1] = 2; A[1][2] = 3;
		A[2][0] = 0; A[2][1] = 0; A[2][2] = 0; // rank deficient

		std::vector<Matrix> plu = PLU_decompose(A);
		Matrix P = plu[0], L = plu[1], U = plu[2];

		REQUIRE(multiply(transpose(P), A) == multiply(L, U));
		check_lower_triangular(L);
		check_upper_triangular(U);
	}

}




//}

//}

//{ Solve System 

TEST_CASE("solve_system tests", "[Matrix][solve_system]") {

	SECTION("Unique solution (2x2 system)") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = number(2); A[0][1] = number(1);
		A[1][0] = number(5); A[1][1] = number(7);

		std::vector<number> b = { number(11), number(13) };

		auto result = solve_system(A, b);

		REQUIRE(result.size() == 2);

		Matrix nullspace = result[0];
		Matrix particular = result[1];

		// Nullspace should be only the zero vector
		REQUIRE(nullspace.cols() == 1);
		for (int r = 0; r < nullspace.rows(); ++r)
			REQUIRE(nullspace[r][0] == number(0));

		// Verify that Ax = b
		Matrix check = multiply(A, particular);
		REQUIRE(check[0][0] == b[0]);
		REQUIRE(check[1][0] == b[1]);
	}

	SECTION("Infinite solutions (rank-deficient 2x2)") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = number(1); A[0][1] = number(2);
		A[1][0] = number(2); A[1][1] = number(4); // row 2 = 2 × row 1

		std::vector<number> b = { number(3), number(6) };

		auto result = solve_system(A, b);

		REQUIRE(result.size() == 2);

		Matrix nullspace = result[0];
		Matrix particular = result[1];

		// Nullspace should include the zero vector and at least one valid vector
		REQUIRE(nullspace.cols() >= 2);
		for (int r = 0; r < nullspace.rows(); ++r)
			REQUIRE(nullspace[r][0] == number(0));

		// Verify that Ax = b for particular solution
		REQUIRE(multiply(A, particular)[0][0] == b[0]);
		REQUIRE(multiply(A, particular)[1][0] == b[1]);

		// Adding any nullspace vector (excluding the null vector) should still satisfy Ax = b
		Matrix test_solution = particular + Matrix::column(nullspace.get_col(1)) * number(3);
		REQUIRE(multiply(A, test_solution) == multiply(A, particular));
	}

	SECTION("Inconsistent system") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = number(1); A[0][1] = number(1);
		A[1][0] = number(1); A[1][1] = number(1); // same row

		std::vector<number> b = { number(1), number(2) }; // inconsistent

		auto result = solve_system(A, b);

		REQUIRE(result.empty()); // No solution
	}

	SECTION("Complex system with unique solution") {
		Matrix A = Matrix::nullmatrix(2, 2);
		A[0][0] = number(1);  A[0][1] = i();
		A[1][0] = i();		  A[1][1] = number(1);

		std::vector<number> b = { number(1, 1), number(2, -1) };

		auto result = solve_system(A, b);
		REQUIRE(result.size() == 2);

		Matrix nullspace = result[0];
		Matrix particular = result[1];

		// Nullspace should be only zero vector
		REQUIRE(nullspace.cols() == 1);
		REQUIRE(nullspace[0][0] == number(0));
		REQUIRE(nullspace[1][0] == number(0));

		// Verify Ax = b
		Matrix check = multiply(A, particular);
		REQUIRE(check[0][0] == b[0]);
		REQUIRE(check[1][0] == b[1]);
	}
}

//}
