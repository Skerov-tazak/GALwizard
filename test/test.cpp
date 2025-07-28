#include "../include/functions.h"
#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <vector>


TEST_CASE("Overloading Equals"){
	
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

TEST_CASE("Basic Arithmetic"){


}
