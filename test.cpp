#include "elliptic.h"

#include <iostream>
#include <cassert>
#include <cstdlib>

bool are_asserts_enabled()
{
  bool asserts_are_enabled = false;
  assert(asserts_are_enabled = true); // intentionally assignment operator, not ==
  return asserts_are_enabled;
}

template <typename real>
real random_nonzero_in_symmetric_unit_interval()
{
  real r;
  do {
    r = real(-1) + real(2) * rand() / RAND_MAX;
  } while (std::abs(r) < real(0.1));
  return r;
}

template <typename real>
bool is_negligible(real a, real scale)
{
  static const real tolerance = std::pow(std::numeric_limits<real>::epsilon(), real(0.75));
  return std::abs(a) <= tolerance * std::abs(scale);
}

template <typename real>
bool is_approximately_equal(real a, real b, real scale)
{
  return is_negligible(a - b, scale);
}

template <typename real>
bool is_acceptable_root_of_cubic(real u, real a, real b, real c, real d)
{
  real v = std::max(real(1), std::abs(u));
  real v2 = v * v;
  real v3 = v2 * v;
  real scale = std::abs(d);
  scale = std::max(scale, std::abs(c * v));
  scale = std::max(scale, std::abs(b * v2));
  scale = std::max(scale, std::abs(a * v3));
  real u2 = u * u;
  real u3 = u2 * u;
  return is_negligible(a * u3 + b * u2 + c * u + d, scale);
}

template <typename real>
void test_one_solve_cubic_equation()
{
  real a = 3 * random_nonzero_in_symmetric_unit_interval<real>();
  real b = 3 * random_nonzero_in_symmetric_unit_interval<real>();
  real c = 3 * random_nonzero_in_symmetric_unit_interval<real>();
  real d = 3 * random_nonzero_in_symmetric_unit_interval<real>();

  int determinant_sign;
  const real uninitialized = real(123456789);
  real u1 = uninitialized,
        u2 = uninitialized,
        u3 = uninitialized,
        negative_sum_u2_u3 = uninitialized,
        product_u2_u3 = uninitialized;
  solve_cubic_equation(a, b, c, d,
                       determinant_sign,
                       u1, u2, u3, negative_sum_u2_u3, product_u2_u3);
  real scale = std::abs(d);
  scale = std::max(scale, std::abs(c));
  scale = std::max(scale, std::abs(b));
  scale = std::max(scale, std::abs(a));
  switch(determinant_sign) {
    case 0:
      assert(u1 == u2 || u1 == u3 || u2 == u3);
      // fall through
    case 1: {
        assert(u1 != uninitialized);
        assert(u2 != uninitialized);
        assert(u2 != uninitialized);
        assert(negative_sum_u2_u3 == uninitialized);
        assert(product_u2_u3 == uninitialized);
        assert(is_acceptable_root_of_cubic(u1, a, b, c, d));
        assert(is_acceptable_root_of_cubic(u2, a, b, c, d));
        assert(is_acceptable_root_of_cubic(u3, a, b, c, d));
      break;
    }
    case -1: {
        assert(u1 != uninitialized);
        assert(u2 == uninitialized);
        assert(u2 == uninitialized);
        assert(negative_sum_u2_u3 != uninitialized);
        assert(product_u2_u3 != uninitialized);
        assert(is_acceptable_root_of_cubic(u1, a, b, c, d));
        assert(is_approximately_equal(b, a * (negative_sum_u2_u3 - u1), scale));
        assert(is_approximately_equal(c, a * (product_u2_u3 - u1 * negative_sum_u2_u3), scale));
        assert(is_approximately_equal(d, -a * u1 * product_u2_u3, scale));
      break;
    }
    default:
      assert(false);
  }
}

void test_solve_cubic_equation(int repetitions)
{
  std::cerr << "testing: solve_cubic_equation" << std::endl;
  for (int i = 0; i < repetitions; i++) {
    test_one_solve_cubic_equation<float>();
    test_one_solve_cubic_equation<double>();
  }
}

int main()
{
  if (!are_asserts_enabled()) {
    std::cerr << "error: asserts are disabled" << std::endl;
    return 1;
  }

  test_solve_cubic_equation(10000000);
}