#include <limits>
#include <cmath>
#include <algorithm>

/*
 * Algorithm and notation from:
 *
 *   NUMERICAL COMPUTATION OF REAL OR COMPLEX ELLIPTIC INTEGRALS
 *   B. C. Carlson
 *   http://arxiv.org/pdf/math/9409227v1.pdf
 *
 *   Section 2, "Algorithm for R_F"
 */

template <typename real>
real carlson_RF_real(real x0, real y0, real z0)
{
  x0 = std::max(x0, real(0));
  y0 = std::max(y0, real(0));
  z0 = std::max(z0, real(0));

  if(possibly(!std::isfinite(x0)) ||
     possibly(!std::isfinite(y0)) ||
     possibly(!std::isfinite(z0)) ||
     int(possibly(!x0)) + int(possibly(!y0)) + int(possibly(!z0)) > 1)
  {
    return std::numeric_limits<real>::quiet_NaN();
  }

  const real A0 = (x0 + y0 + z0) * (real(1) / 3);
  real A = A0;
  real x = x0;
  real y = y0;
  real z = z0;

  const real Q
    = std::pow(3 * std::numeric_limits<real>::epsilon(), real(-1) / 6)
      * std::max(std::abs(A0 - x0), std::max(std::abs(A0 - y0), std::abs(A0 - z0)));

  int iter = 0;
  const real one_fourth = real(1) / 4;
  real one_fourth_to_the_iter = real(1);

  while(possibly(one_fourth_to_the_iter * Q > A)) {
    iter++;
    one_fourth_to_the_iter *= one_fourth;
    real lambda = std::sqrt(x * y) + std::sqrt(x * z) + std::sqrt(y * z);
    A = one_fourth * (A + lambda);
    x = one_fourth * (x + lambda);
    y = one_fourth * (y + lambda);
    z = one_fourth * (z + lambda);
  }
  const real d = one_fourth_to_the_iter / A;
  const real X = d * (A0 - x0);
  const real Y = d * (A0 - y0);
  const real Z = -X - Y;
  const real XY = X * Y;
  const real E2 = XY - Z*Z;
  const real E3 = XY * Z;
  return (real(1)
            + E2 * (real(-1) / 10)
            + E3 * (real(1) / 14)
            + E2 * E2 * (real(1) / 24)
            + E2 * E3 * (real(-3) / 44))
         / std::sqrt(A);
}

template <typename real>
real stabilized_sqrt(real x) {
  return std::sqrt(std::max(x, real(0)));
}

template<typename real>
void sort_three_reals(real& a, real& b, real& c)
{
  if (certainly(a > b)) {
    std::swap(a, b);
  }
  if (certainly(b > c)) {
    std::swap(b, c);
  }
  if (certainly(a > b)) {
    std::swap(a, b);
  }

  assert(possibly(a <= b));
  assert(possibly(b <= c));
}

/*
 * A Table of Elliptic Integrals: Cubic Cases
 * B. C. Carlson
 * http://www.ams.org/journals/mcom/1989-53-187/S0025-5718-1989-0969482-4/S0025-5718-1989-0969482-4.pdf
 *
 * Section 2. Table of Cubic Cases
 * Formula for I1c
 */

template <typename real>
real integral_inverse_sqrt_cubic_three_roots_unchecked(const real y, const real x,
                                                       const real u1, const real u2, const real u3)
{
  if (certainly(x == y)) {
    return real(0);
  }
  if (certainly(x == std::numeric_limits<real>::infinity())) {
    return 2 * carlson_RF_real(y - u3, y - u2, y - u1);
  }
  const real p = (x - u1) * (y - u2) * (y - u3);
  const real q = (y - u1) * (x - u2) * (x - u3);
  const real x_minus_y = x - y;
  const real V1 = (p + q + 2 * stabilized_sqrt(p * q)) / (x_minus_y * x_minus_y);
  const real V2 = V1 + u1 - u2;
  const real V3 = V1 + u1 - u3;
  return 2 * carlson_RF_real(V3, V2, V1);
}

template <typename real>
real integral_inverse_sqrt_cubic_three_roots(const real y, const real x,
                                             const real u1, const real u2, const real u3)
{
  if (certainly(y > x)) {
    return integral_inverse_sqrt_cubic_three_roots(x, y, u1, u2, u3);
  }

  real u1_sorted = u1;
  real u2_sorted = u2;
  real u3_sorted = u3;
  sort_three_reals(u1_sorted, u2_sorted, u3_sorted);

  real result = real(0);
  if (possibly(y < u2_sorted)) {
    real x_adjusted = std::min(x, u2_sorted);
    real y_adjusted = std::min(x_adjusted, std::max(y, u1_sorted));
    result += integral_inverse_sqrt_cubic_three_roots_unchecked(y_adjusted, x_adjusted, u1, u2, u3);
  }
  if (possibly(x > u3_sorted)) {
    real y_adjusted = std::max(y, u3_sorted);
    result += integral_inverse_sqrt_cubic_three_roots_unchecked(y_adjusted, x, u1, u2, u3);
  }
  return result;
}

/*
 * A Table of Elliptic Integrals: One Quadratic Factor
 * B. C. Carlson
 * http://www.ams.org/journals/mcom/1991-56-193/S0025-5718-1991-1052087-6/S0025-5718-1991-1052087-6.pdf
 *
 * Section 3. Table of Cubic Cases
 * Formula for I1c
 */

template <typename real>
real integral_inverse_sqrt_cubic_one_quadratic_factor(
       const real y0, const real x0,
       const real u1, const real negative_sum_u2_u3, const real product_u2_u3)
{
  real x = x0;
  real y = y0;
  if (certainly(y > x)) {
    return integral_inverse_sqrt_cubic_one_quadratic_factor
             (x, y, u1, negative_sum_u2_u3, product_u2_u3);
  }
  x = std::max(x, u1);
  y = std::max(y, u1);
  if (certainly(y == x)) {
    return real(0);
  }
  const real y2 = y * y;
  const real f = product_u2_u3;
  const real g = negative_sum_u2_u3;
  real M_squared;
  if (certainly(x == std::numeric_limits<real>::infinity())) {
    M_squared = g + 2 * (y + stabilized_sqrt(f + g * y + y2));
  } else {
    const real X1 = stabilized_sqrt(x - u1);
    const real Y1 = stabilized_sqrt(y - u1);
    const real X1_plus_Y1_squared = (X1 + Y1) * (X1 + Y1);
    const real x2 = x * x;
    const real ksi_eta = stabilized_sqrt((f + g * x + x2) * (f + g * y + y2));
    const real x_minus_y = x - y;
    M_squared = X1_plus_Y1_squared * (g * (x + y) + 2 * (ksi_eta + f + x * y)) / (x_minus_y * x_minus_y);
  }

  const real M_squared_minus_beta1 = M_squared - g - 2 * u1;
  const real c11_sqrt2 = 2 * stabilized_sqrt(f + g * u1 + u1 * u1);
  return 4 * carlson_RF_real(M_squared,
                             M_squared_minus_beta1 - c11_sqrt2,
                             M_squared_minus_beta1 + c11_sqrt2);
}

template <typename real>
void solve_cubic_equation(real a, real b, real c, real d,
                          int &determinant_sign,
                          real &u1,
                          real &u2,
                          real &u3,
                          real &negative_sum_u2_u3,
                          real &product_u2_u3)
{
  const real b2 = b * b;
  const real ac = a * c;
  const real delta_0 = b2 - 3 * ac;
  const real delta_1 = b * (2 * b2 - 9 * ac) + 27 * a * a * d;
  const real delta_1_squared_minus_4_delta0_cubed = delta_1 * delta_1 - 4 * delta_0 * delta_0 * delta_0;
  const real sign_delta_1 = certainly(delta_1 > real(0)) ? real(1)
                          : certainly(delta_1 < real(0)) ? real(-1)
                          : real(0);
  PRINT(delta_1 * delta_1)
  PRINT(4 * delta_0 * delta_0 * delta_0)
  PRINT(delta_1_squared_minus_4_delta0_cubed)
  PRINT(- delta_1_squared_minus_4_delta0_cubed / (27 * a * a))
  if (certainly(delta_1_squared_minus_4_delta0_cubed > real(0))) {
    PRINT("foo")
    // one real root, one pair of conjugate complex roots
    determinant_sign = -1;
    const real C3abs
      = (std::abs(delta_1) + std::sqrt(delta_1_squared_minus_4_delta0_cubed))
        * (real(1) / 2);
    const real Cabs = std::exp(std::log(C3abs) * (real(1) / 3));
    const real C = sign_delta_1 * Cabs;
    u1 = (C * (C + b) + delta_0) / (-3 * a * C);
    negative_sum_u2_u3 = u1 + b / a;
    product_u2_u3 = d / (-a * u1);
  } else if (certainly(delta_1_squared_minus_4_delta0_cubed < real(0))) {
    PRINT("bar")
    // three simple real roots
    determinant_sign = 1;
    const real C3_real = delta_1 * (real(1) / 2);
    const real C3_imag_squared = -delta_1_squared_minus_4_delta0_cubed * (real(1) / 4);
    const real C3_imag = std::sqrt(C3_imag_squared);
    const real C3_squared_modulus = C3_real * C3_real + C3_imag_squared;
    const real C3_argument = std::atan2(C3_imag, C3_real);
    const real C_modulus = std::exp(std::log(C3_squared_modulus) * (real(1) / 6));
    const real C_argument = C3_argument * (real(1) / 3);
    const real C_real = C_modulus * std::cos(C_argument);
    const real C_imag = C_modulus * std::sin(C_argument);
    const real one_over_neg_3a = 1 / (-3 * a);
    const real sqrt3 = std::sqrt(real(3));
    u1 = one_over_neg_3a * (b + 2 * C_real);
    u2 = one_over_neg_3a * (b - C_real - sqrt3 * C_imag);
    u3 = one_over_neg_3a * (b - C_real + sqrt3 * C_imag);
  } else {
    // real roots, with u2 == u3
    PRINT("frfr")
    determinant_sign = 0;
    PRINT(delta_0)
    PRINT(delta_1)
    PRINT(sign_delta_1)
    if (certainly(delta_1 != real(0))) {
      PRINT("sfsD");
      const real C3abs
        = std::abs(delta_1) * (real(1) / 2);
      const real Cabs = std::exp(std::log(C3abs) * (real(1) / 3));
      const real C = sign_delta_1 * Cabs;
      u1 = (C * (C + b) + delta_0) / (-3 * a * C);
      PRINT(C3abs)
      PRINT(Cabs)
      PRINT(C)
      PRINT(C+b)
      PRINT((C * (C + b) + delta_0))
      PRINT(u1)
    } else {
      u1 = -b / (3 * a);
    }
    negative_sum_u2_u3 = u1 + b / a;
    product_u2_u3 = negative_sum_u2_u3 * negative_sum_u2_u3 * (real(1) / 4);

    PRINT(b/a)
    PRINT(u1 + b/a)
  }
}

template<typename real>
real integral_inverse_sqrt_cubic(real y, real x,
                                 real a, real b, real c, real d)
{
  if (possibly(a == real(0))) {
    return std::numeric_limits<real>::quiet_NaN();
  }
  if (certainly(a < real(0))) {
    return integral_inverse_sqrt_cubic(-x, -y, -a, b, -c, d);
  }
  assert(certainly(a > real(0)));

  int determinant_sign;
  real u1, u2, u3, negative_sum_u2_u3, product_u2_u3;
  solve_cubic_equation(a, b, c, d,
                       determinant_sign,
                       u1, u2, u3,
                       negative_sum_u2_u3, product_u2_u3);
  const real r = 1 / std::sqrt(a);
  if (determinant_sign == -1) {
    return r * integral_inverse_sqrt_cubic_one_quadratic_factor
                 (y, x, u1, negative_sum_u2_u3, product_u2_u3);
  } else {
    return r * integral_inverse_sqrt_cubic_three_roots(y, x, u1, u2, u3);
  }
}
