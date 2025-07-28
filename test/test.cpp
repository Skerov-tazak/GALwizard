#include "../include/functions.h"
#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <vector>


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

            std::vector<float> v = {0, 0, 0};
            std::vector<float> expected = {0, 0};

            std::vector<float> result = multiply(A, v);

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

        float scalar = 2.0f;

        Matrix expected = Matrix::nullmatrix(2, 2);
        expected[0][0] = 2; expected[0][1] = -4;
        expected[1][0] = 6; expected[1][1] = -8;

        REQUIRE(A * scalar == expected);
        REQUIRE(scalar * A == expected);
        REQUIRE(expected / scalar == A);
    }
}


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
        expected[0][0] = 0.6f;  expected[0][1] = -0.7f;
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

TEST_CASE("Nullspace", "[Matrix][nullspace]") {

    SECTION("Full-rank square matrix → only zero vector in nullspace") {
        Matrix A = Matrix::nullmatrix(2, 2);
        A[0][0] = 1; A[0][1] = 2;
        A[1][0] = 3; A[1][1] = 4;

        Matrix ns = nullspace(A);

        // Should return 2x1 matrix, all zero
        Matrix expected = Matrix::nullmatrix(2, 1);
        REQUIRE(ns == expected);
    }

    SECTION("Matrix with a 1D nullspace") {
        // A*x = 0 has solution x = (1, -1)
        Matrix A = Matrix::nullmatrix(1, 2);
        A[0][0] = 1; A[0][1] = 1;

        Matrix ns = nullspace(A);

        // Expected: nullspace has basis vectors (0,0) and (1,-1)
        Matrix expected = Matrix::nullmatrix(2, 2);
        expected[0][0] = 0; expected[1][0] = 0;  // zero vector first
        expected[0][1] = 1; expected[1][1] = -1; // basis vector

        REQUIRE(ns == expected);
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

        REQUIRE(ns == expected);
    }

    SECTION("Matrix with no rows (nullspace = ℝ^n)") {
        Matrix A = Matrix::nullmatrix(0, 2); // 0x2 matrix

        Matrix ns = nullspace(A);

        // Nullspace = ℝ² → zero vector + identity basis vectors
        Matrix expected = Matrix::nullmatrix(2, 3);
        expected[0][0] = 0; expected[1][0] = 0;
        expected[0][1] = 1; expected[1][1] = 0;
        expected[0][2] = 0; expected[1][2] = 1;

        REQUIRE(ns == expected);
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

        REQUIRE(ns == expected);
    }
}

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

TEST_CASE("QR decomposition via Householder", "[Matrix][QR_decomposeHS]") {
    SECTION("3x2 matrix") {
        Matrix A = Matrix::nullmatrix(3, 2);
        A[0][0] = 12; A[1][0] = -51; A[2][0] = 4;
        A[0][1] = 6;  A[1][1] = 167; A[2][1] = -68;

        std::vector<Matrix> qr = QR_decomposeHS(A);
        REQUIRE(qr.size() == 2);

        Matrix Q = qr[0], R = qr[1];

        REQUIRE(multiply(transpose(Q), Q) == Matrix::identity(3));
        REQUIRE(multiply(Q, R) == A);
    }
}

TEST_CASE("PLU decomposition", "[Matrix][PLU_decompose]") {
    SECTION("3x3 example matrix") {
        Matrix A = Matrix::nullmatrix(3, 3);
        A[0][0] = 2; A[0][1] = 0; A[0][2] = 2;
        A[1][0] = 1; A[1][1] = 2; A[1][2] = 3;
        A[2][0] = 3; A[2][1] = 1; A[2][2] = 4;

        std::vector<Matrix> plu = PLU_decompose(A);
        REQUIRE(plu.size() == 3);

        Matrix P = plu[0], L = plu[1], U = plu[2];

        REQUIRE(multiply(P, A) == multiply(L, U));
    }

    SECTION("Identity matrix stays the same") {
        Matrix I = Matrix::identity(3);
        std::vector<Matrix> plu = PLU_decompose(I);

        Matrix P = plu[0], L = plu[1], U = plu[2];
        REQUIRE(P == I);
        REQUIRE(L == I);
        REQUIRE(U == I);
    }
}
