#include <iostream>
#define PRINT(x) //std::cerr << #x << " = " << (x) << std::endl;

#include "interval.h"
#include "elliptic.h"

#include <cassert>
#include <cstdlib>

bool are_asserts_enabled()
{
  bool asserts_are_enabled = false;
  assert(asserts_are_enabled = true); // intentionally assignment operator, not ==
  return asserts_are_enabled;
}

void verify(bool_interval b)
{
  assert(possibly(b));
}

template <typename real>
real random_nonzero_in_symmetric_unit_interval()
{
  real r = real(-1) + real(2) * rand() / RAND_MAX;
  r.fuzz();
  return r;
}

template <typename real>
real random_positive_in_unit_interval()
{
  real r = real(rand()) / RAND_MAX;
  r.fuzz();
  return r;
}

template <typename real>
bool_interval is_root_of_cubic(real u, real a, real b, real c, real d)
{
  real u2 = u * u;
  real u3 = u2 * u;
  return a * u3 + b * u2 + c * u + d == real(0);
}

template <typename real>
void test_one_solve_cubic_equation()
{
  real a = random_nonzero_in_symmetric_unit_interval<real>();
  real b = random_nonzero_in_symmetric_unit_interval<real>();
  real c = random_nonzero_in_symmetric_unit_interval<real>();
  real d = random_nonzero_in_symmetric_unit_interval<real>();

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
  switch(determinant_sign) {
    case 1: {
        PRINT(a);
        PRINT(b);
        PRINT(c);
        PRINT(d);
        PRINT(u1);
        PRINT(u1*u1);
        PRINT(u1*u1*u1);
        PRINT(a*u1*u1*u1 + b*u1*u1 + c*u1 + d);
        PRINT(u2);
        PRINT(u2*u2);
        PRINT(u2*u2*u2);
        PRINT(a*u2*u2*u2 + b*u2*u2 + c*u2 + d);
        PRINT(u3);
        PRINT(u3*u3);
        PRINT(u3*u3*u3);
        PRINT(a*u3*u3*u3 + b*u3*u3 + c*u3 + d);
        verify(u1 != uninitialized);
        verify(u2 != uninitialized);
        verify(u2 != uninitialized);
        verify(negative_sum_u2_u3 == uninitialized);
        verify(product_u2_u3 == uninitialized);
        verify(is_root_of_cubic(u1, a, b, c, d));
        verify(is_root_of_cubic(u2, a, b, c, d));
        verify(is_root_of_cubic(u3, a, b, c, d));
      break;
    }
    case 0:
    case -1: {
        verify(u1 != uninitialized);
        verify(u2 == uninitialized);
        verify(u2 == uninitialized);
        verify(negative_sum_u2_u3 != uninitialized);
        verify(product_u2_u3 != uninitialized);
        PRINT(negative_sum_u2_u3);
        PRINT(negative_sum_u2_u3 * negative_sum_u2_u3);
        PRINT(product_u2_u3)
        PRINT(-a)
        PRINT(-a * u1)
        PRINT(possibly(-a * u1 == float(0)))
        PRINT(d / (-a * u1))
        PRINT(4 * product_u2_u3)
        verify(negative_sum_u2_u3 * negative_sum_u2_u3 <= 4 * product_u2_u3);
        PRINT(a);
        PRINT(b);
        PRINT(c);
        PRINT(d);
        PRINT(u1);
        PRINT(u1*u1);
        PRINT(u1*u1*u1);
        PRINT(a*u1*u1*u1 + b*u1*u1 + c*u1 + d);
        verify(is_root_of_cubic(u1, a, b, c, d));
        verify(b == a * (negative_sum_u2_u3 - u1));
        verify(c == a * (product_u2_u3 - u1 * negative_sum_u2_u3));
        PRINT(d)
        PRINT(a)
        PRINT(u1)
        PRINT(product_u2_u3)
        PRINT(-a * u1 * product_u2_u3)
        //verify(d == -a * u1 * product_u2_u3);
      break;
    }
    default:
      verify(false);
  }
}

void test_solve_cubic_equation(int repetitions)
{
  std::cerr << "testing: solve_cubic_equation" << std::endl;
  for (int i = 0; i < repetitions; i++) {
    test_one_solve_cubic_equation<interval<float> >();
    //test_one_solve_cubic_equation<float>();
    //test_one_solve_cubic_equation<double>();
  }
}

template <typename real>
void test_one_integral_inverse_sqrt_cubic()
{
  const real a = std::abs(random_nonzero_in_symmetric_unit_interval<real>());
  real b = random_nonzero_in_symmetric_unit_interval<real>();
  real c = random_nonzero_in_symmetric_unit_interval<real>();
  real d = random_nonzero_in_symmetric_unit_interval<real>();

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

  real y = random_positive_in_unit_interval<real>();
  real x = y;
  real sum_of_slices = real(0);

  for (int i = 0; i < 10; i++) {
    real delta = random_positive_in_unit_interval<real>();
    real new_x = x + delta;

    real slice = integral_inverse_sqrt_cubic(x, new_x, a, b, c, d);
    verify(std::isfinite(slice));
    sum_of_slices += slice;
    real whole_integral = integral_inverse_sqrt_cubic(y, new_x, a, b, c, d);
    verify(sum_of_slices == whole_integral);
    x = new_x;
  }
}

void test_integral_inverse_sqrt_cubic(int repetitions)
{
  std::cerr << "testing: integral_inverse_sqrt_cubic" << std::endl;
  for (int i = 0; i < repetitions; i++) {
    test_one_integral_inverse_sqrt_cubic<interval<float> >();

    //test_one_integral_inverse_sqrt_cubic<float>();
    //test_one_integral_inverse_sqrt_cubic<double>();
    
  }
}

int main()
{
  std::cerr.precision(16);
  if (!are_asserts_enabled()) {
    std::cerr << "error: asserts are disabled" << std::endl;
    return 1;
  }

  test_solve_cubic_equation(10000000);
  test_integral_inverse_sqrt_cubic(1000000);
  return 0;
}
