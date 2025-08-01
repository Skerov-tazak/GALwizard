#include "number.h"
#include <catch2/catch_test_macros.hpp>


number from_polar(double, double);


TEST_CASE("Basic arithmetic operations", "[number]") {
    number a(3, 4);
    number b(1, -2);

    SECTION("Addition") {
        REQUIRE(a + b == number(4, 2));
        REQUIRE(a + 5 == number(8, 4));
        REQUIRE(5 + a == number(8, 4));
    }

    SECTION("Subtraction") {
        REQUIRE(a - b == number(2, 6));
        REQUIRE(a - 2 == number(1, 4));
        REQUIRE(7 - a == number(4, -4));
    }

    SECTION("Multiplication") {
        REQUIRE(a * b == number(11, -2)); // (3+4i)*(1-2i) = 11-2i
        REQUIRE(a * 2 == number(6, 8));
        REQUIRE(2 * a == number(6, 8));
    }

    SECTION("Division") {
        REQUIRE(a / b == number(-1, 2)); // (3+4i)/(1-2i) = -1+2i
        REQUIRE(a / 2 == number(1.5, 2));
    }
}

TEST_CASE("Comparison operators", "[number]") {
    number a(1.000000000001, 2.0);
    number b(1.0, 2.0);
    REQUIRE(a == b);
    REQUIRE_FALSE(a != b);

 	REQUIRE(number(0.0, 0.0) == 0.0);
    REQUIRE(0.0 == number(0.0, 0.0));
}

TEST_CASE("Conjugate and normalization", "[number]") {
    number a(3, 4);
    REQUIRE(conj(a) == number(3, -4));

    number n = a.normalised();
    REQUIRE(n.radius_squared() == 1.0);
    REQUIRE(n == number(0.6,0.8));
}

TEST_CASE("Rotation", "[number]") {
    number a(1, 0);
    number rotated = a.rotate(M_PI / 2);
    REQUIRE(rotated == number(0.0,1.0));
}

TEST_CASE("Polar form conversion", "[number]") {
    number a(3, 4);
    auto polar = to_polar(a);
    REQUIRE(polar.first == 5.0);
    REQUIRE(polar.second == std::atan2(4, 3));

    number b = from_polar(5, std::atan2(4, 3));
    REQUIRE(b == number(3.0,4.0));
}

TEST_CASE("Power and roots", "[number]") {
    number a(3, 4);
    number squared = a.power(number(2, 0));
    REQUIRE(squared == number(-7, 24));

    auto r = a.roots(2);
    REQUIRE(r.size() == 2);
    for (auto& root : r) {
        REQUIRE(root * root == a);
    }
}

TEST_CASE("Trigonometric functions", "[number]") {
    number z(0, M_PI/2);

    number s = sin(z);
    number c = cos(z);
    number t = tan(z);

    REQUIRE(t == s / c);
}

TEST_CASE("Exponential and logarithmic functions", "[number]") {
    number z(2, 1);

    number e = exp(z);
    number l = ln(e);

    REQUIRE(l == z); // ln(exp(z)) == z
}

TEST_CASE("Mathematical identities", "[number]") {
    number z(3, 4);

    // z * conj(z) = |z|^2
    REQUIRE(z * conj(z) == z.radius_squared());

    // rotation preserves radius
    number rotated = z.rotate(M_PI / 3);
    REQUIRE(approx_equal(rotated.radius_squared(), z.radius_squared()));
}

TEST_CASE("Edge cases", "[number]") {
    number zero(0, 0);
    REQUIRE(approx_equal(zero.radius_squared(), 0.0));
    REQUIRE(approx_equal(zero.angle(), 0.0));

    number purely_real(5, 0);
    REQUIRE(approx_equal(purely_real.angle(),0.0));

    number purely_imag(0, 2);
    REQUIRE(approx_equal(purely_imag.angle(), M_PI/2));

    number neg_real(-5, 0);
    REQUIRE(approx_equal(neg_real.angle(), M_PI));
}
TEST_CASE("Zero and near-zero handling", "[number][edge]") {
    number z(0.0, 0.0);
    REQUIRE(z == 0.0);
    REQUIRE(z.radius_squared() == 0.0);
    REQUIRE(z.angle() == 0.0);

    number tiny(1e-14, -1e-14);
    REQUIRE(tiny == 0.0); // should be considered approximately zero
    REQUIRE(tiny.is_real());
    REQUIRE(tiny.is_imag());
}

TEST_CASE("Very large numbers", "[number][edge]") {
    number big(1e150, 1e150);
    REQUIRE(std::isfinite(big.radius()));
    REQUIRE(big * big != 0.0); // should not underflow to zero
}

TEST_CASE("Very small numbers", "[number][edge]") {
    number small(1e-150, 1e-150);
    REQUIRE(small.radius() < 1e-149);
    REQUIRE(small + small == number(2e-150, 2e-150));
}

TEST_CASE("Negative radius handling in polar form", "[number][edge]") {
    double r = 5.0;
    double theta = -M_PI / 4;
    number z = from_polar(r, theta);

    auto polar = to_polar(z);
    REQUIRE(approx_equal(polar.first, r));
    REQUIRE(approx_equal(polar.second, theta));
}

TEST_CASE("Rotation invariants", "[number][edge]") {
    number z(3, 4);
    double r2 = z.radius_squared();

    number rotated = z.rotate(M_PI);
    REQUIRE(rotated.radius_squared() == r2); // rotation preserves radius
    REQUIRE(rotated == number(-3, -4));      // 180° rotation negates vector

    number full_rotation = z.rotate(2 * M_PI);
    REQUIRE(full_rotation == z);
}

TEST_CASE("Conjugate and multiplication identity", "[number][edge]") {
    number z(3, 4);
    REQUIRE(z * conj(z) == z.radius_squared());
    REQUIRE(conj(conj(z)) == z);
}

TEST_CASE("Power and roots identity", "[number][edge]") {
    number z(3, 4);

    // Square root squared = original
    auto roots = z.roots(2);
    for (const auto& r : roots) {
        REQUIRE(r * r == z);
    }

    // Power with integer exponent
    REQUIRE(z.power(number(2, 0)) == z * z);

    // z^0 == 1
    REQUIRE(z.power(number(0, 0)) == 1.0);
}

TEST_CASE("Logarithmic and exponential inverses", "[number][edge]") {
    number z(2, 1);
    REQUIRE(ln(exp(z)) == z); // ln(exp(z)) == z

    // exp(ln(z)) == z
    REQUIRE(exp(ln(z)) == z);
}

TEST_CASE("Trigonometric symmetries", "[number][edge]") {
    number z(M_PI/2, 0); // purely real

    REQUIRE(sin(z) == 1.0);
    REQUIRE(cos(z) == 0.0);

    // sin(-z) == -sin(z), cos(-z) == cos(z)
    REQUIRE(sin(-z) == -sin(z));
    REQUIRE(cos(-z) == cos(z));
}


