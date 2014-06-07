#ifndef interval_h_
#define interval_h_

#include "interval-ready.h"

#include <limits>
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cassert>

#include <iostream>

template <typename real>
class interval
{
  real m_lo, m_hi;

public:

  interval(real lo, real hi)
    : m_lo(lo)
    , m_hi(hi)
  {
    assert(lo <= hi ||
           std::isnan(lo) ||
           std::isnan(hi));
  }

  interval& operator=(real a)
  {
    m_lo = m_hi = a;
    return *this;
  }

  interval(real a)
  {
    *this = a;
  }

  interval()
  {
    *this = std::numeric_limits<real>::quiet_NaN();
  }

  real lo() const
  {
    return m_lo;
  }

  real hi() const
  {
    return m_hi;
  }

  real mid() const
  {
    return real(0.5) * (lo() + hi());
  }

  interval lower_half() const
  {
    return interval(lo(), mid());
  }

  interval upper_half() const
  {
    return interval(mid(), hi());
  }

  real diameter() const
  {
    return hi() - lo();
  }

  real radius() const
  {
    return real(0.5) * diameter();
  }

  static interval construct_NaN()
  {
    return std::numeric_limits<real>::quiet_NaN();
  }

  static interval construct_quantum(real x)
  {
    const int saved_rounding_mode = std::fegetround();
    std::fesetround(FE_UPWARD);
    real y = x + std::numeric_limits<real>::min();
    std::fesetround(saved_rounding_mode);
    return interval(x, y);
  }

  bool is_quantum() const
  {
    const int saved_rounding_mode = std::fegetround();
    std::fesetround(FE_UPWARD);
    real lo_plus_quantum = lo() + std::numeric_limits<real>::min();
    std::fesetround(saved_rounding_mode);
    return hi() == lo_plus_quantum;
  }
};

template <typename real>
interval<real> interval_union(const interval<real>& i, const interval<real>& j)
{
  return interval<real>(std::min(i.lo(), j.lo()), std::max(i.hi(), j.hi()));
}

class bool_interval
{
  bool m_lo, m_hi;

  void assert_consistent() const
  {
    assert(!m_lo || m_hi);
  }

public:

  bool_interval()
    : m_lo(false)
    , m_hi(true)
  {}

  bool_interval(bool lo, bool hi)
    : m_lo(lo)
    , m_hi(hi)
  {
    assert_consistent();
  }

  bool_interval& operator=(bool value)
  {
    m_lo = value;
    m_hi = true;
    return *this;
  }

  bool_interval(bool value)
  {
    *this = value;
  }

  bool_interval operator!() const
  {
    return bool_interval(!m_hi, !m_lo);
  }

  bool lo() const { return m_lo; }
  bool hi() const { return m_hi; }
};

inline bool certainly(const bool_interval& b)
{
  return b.lo();
}

inline bool possibly(const bool_interval& b)
{
  return b.hi();
}

inline bool_interval operator&&(const bool_interval& a,
                                const bool_interval& b)
{
  return bool_interval(a.lo() && b.lo(), a.hi() && b.hi());
}

inline bool_interval operator||(const bool_interval& a,
                                const bool_interval& b)
{
  return bool_interval(a.lo() || b.lo(), a.hi() || b.hi());
}

template <typename real>
bool_interval operator>=(const interval<real>& a, const interval<real>& b)
{
  if (std::isnan(a) || std::isnan(b)) {
    return bool_interval(false, true);
  }
  return bool_interval(a.lo() >= b.hi(), a.hi() >= b.lo());
}

template <typename real>
bool_interval operator<=(const interval<real>& a, const interval<real>& b)
{
  return b >= a;
}

template <typename real>
bool_interval operator==(const interval<real>& a, const interval<real>& b)
{
  return a >= b && b >= a;
}

template <typename real>
bool_interval operator>(const interval<real>& a, const interval<real>& b)
{
  return !(b >= a);
}

template <typename real>
bool_interval operator<(const interval<real>& a, const interval<real>& b)
{
  return !(a >= b);
}

template <typename real>
bool_interval operator!=(const interval<real>& a, const interval<real>& b)
{
  return !(a == b);
}

template <typename real>
bool_interval operator!(const interval<real>& a)
{
  return a == real(0);
}

template <typename real>
bool_interval operator>=(const interval<real>& a, real b)
{
  if (std::isnan(a) || std::isnan(b)) {
    return bool_interval(false, true);
  }
  return bool_interval(a.lo() >= b, a.hi() >= b);
}

template <typename real>
bool_interval operator>=(real a, const interval<real>& b)
{
  if (std::isnan(a) || std::isnan(b)) {
    return bool_interval(false, true);
  }
  return bool_interval(a >= b.hi(), a >= b.lo());
}

template <typename real>
bool_interval operator<=(const interval<real>& a, real b)
{
  return b >= a;
}

template <typename real>
bool_interval operator<=(real a, const interval<real>& b)
{
  return b >= a;
}

template <typename real>
bool_interval operator==(const interval<real>& a, real b)
{
  return a >= b && b >= a;
}

template <typename real>
bool_interval operator==(real a, const interval<real>& b)
{
  return a >= b && b >= a;
}

template <typename real>
bool_interval operator>(const interval<real>& a, real b)
{
  return !(b >= a);
}

template <typename real>
bool_interval operator>(real a, const interval<real>& b)
{
  return !(b >= a);
}

template <typename real>
bool_interval operator<(const interval<real>& a, real b)
{
  return !(a >= b);
}

template <typename real>
bool_interval operator<(real a, const interval<real>& b)
{
  return !(a >= b);
}

template <typename real>
bool_interval operator!=(const interval<real>& a, real b)
{
  return !(a == b);
}

template <typename real>
bool_interval operator!=(real a, const interval<real>& b)
{
  return !(a == b);
}

template <typename real, typename other_type>
interval<real>& operator+=(interval<real>& i, const other_type& other)
{
  return i = i + other;
}

template <typename real, typename other_type>
interval<real>& operator-=(interval<real>& i, const other_type& other)
{
  return i = i - other;
}

template <typename real, typename other_type>
interval<real>& operator*=(interval<real>& i, const other_type& other)
{
  return i = i * other;
}

template <typename real, typename other_type>
interval<real>& operator/=(interval<real>& i, const other_type& other)
{
  return i = i / other;
}

template <typename real>
interval<real> operator-(const interval<real>& i)
{
  return interval<real>(-i.hi(), -i.lo());
}

template <typename real>
interval<real> operator+(const interval<real>& i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i.lo() + j.lo();
  std::fesetround(FE_UPWARD);
  const real b = i.hi() + j.hi();
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator-(const interval<real>& i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i.lo() - j.hi();
  std::fesetround(FE_UPWARD);
  const real b = i.hi() - j.lo();
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator*(const interval<real>& i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real all = i.lo() * j.lo();
  const real ahl = i.hi() * j.lo();
  const real alh = i.lo() * j.hi();
  const real ahh = i.hi() * j.hi();
  const real a = std::min(std::min(all, ahl), std::min(alh, ahh));
  std::fesetround(FE_UPWARD);
  const real bll = i.lo() * j.lo();
  const real bhl = i.hi() * j.lo();
  const real blh = i.lo() * j.hi();
  const real bhh = i.hi() * j.hi();
  const real b = std::max(std::max(bll, bhl), std::max(blh, bhh));
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator/(const interval<real>& i, const interval<real>& j)
{
  if (possibly(j == real(0))) {
    return interval<real>::construct_NaN();
  }

  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real all = i.lo() / j.lo();
  const real ahl = i.hi() / j.lo();
  const real alh = i.lo() / j.hi();
  const real ahh = i.hi() / j.hi();
  const real a = std::min(std::min(all, ahl), std::min(alh, ahh));
  std::fesetround(FE_UPWARD);
  const real bll = i.lo() / j.lo();
  const real bhl = i.hi() / j.lo();
  const real blh = i.lo() / j.hi();
  const real bhh = i.hi() / j.hi();
  const real b = std::max(std::max(bll, bhl), std::max(blh, bhh));
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator+(const interval<real>& i, real j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i.lo() + j;
  std::fesetround(FE_UPWARD);
  const real b = i.hi() + j;
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator-(const interval<real>& i, real j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i.lo() - j;
  std::fesetround(FE_UPWARD);
  const real b = i.hi() - j;
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator*(const interval<real>& i, real j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::min(i.lo() * j, i.hi() * j);
  std::fesetround(FE_UPWARD);
  const real b = std::max(i.lo() * j, i.hi() * j);
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator/(const interval<real>& i, real j)
{
  if (j == real(0)) {
    return interval<real>::construct_NaN();
  }

  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::min(i.lo() / j, i.hi() / j);
  std::fesetround(FE_UPWARD);
  const real b = std::max(i.lo() / j, i.hi() / j);
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator+(real i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i + j.lo();
  std::fesetround(FE_UPWARD);
  const real b = i + j.hi();
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator-(real i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = i - j.hi();
  std::fesetround(FE_UPWARD);
  const real b = i - j.lo();
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator*(real i, const interval<real>& j)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::min(i * j.lo(), i * j.hi());
  std::fesetround(FE_UPWARD);
  const real b = std::max(i * j.lo(), i * j.hi());
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator/(real i, const interval<real>& j)
{
  if (possibly(j == real(0))) {
    return interval<real>::construct_NaN();
  }

  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::min(i / j.lo(), i / j.hi());
  std::fesetround(FE_UPWARD);
  const real b = std::max(i / j.lo(), i / j.hi());
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> operator+(const interval<real>& i, int j)
{
  return i + real(j);
}

template <typename real>
interval<real> operator-(const interval<real>& i, int j)
{
  return i - real(j);
}

template <typename real>
interval<real> operator*(const interval<real>& i, int j)
{
  return i * real(j);
}

template <typename real>
interval<real> operator/(const interval<real>& i, int j)
{
  return i / real(j);
}

template <typename real>
interval<real> operator+(int i, const interval<real>& j)
{
  return real(i) + j;
}

template <typename real>
interval<real> operator-(int i, const interval<real>& j)
{
  return real(i) - j;
}

template <typename real>
interval<real> operator*(int i, const interval<real>& j)
{
  return real(i) * j;
}

template <typename real>
interval<real> operator/(int i, const interval<real>& j)
{
  return real(i) / j;
}

namespace std {

template <typename real>
interval<real> sqrt(const interval<real>& i)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::sqrt(i.lo());
  std::fesetround(FE_UPWARD);
  const real b = std::sqrt(i.hi());
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> exp(const interval<real>& i)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::exp(i.lo());
  std::fesetround(FE_UPWARD);
  const real b = std::exp(i.hi());
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> log(const interval<real>& i)
{
  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = std::log(i.lo());
  std::fesetround(FE_UPWARD);
  const real b = std::log(i.hi());
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> cos(const interval<real>& i)
{
  const real width = i.hi() - i.lo();

  if (width >= real(2 * M_PI)) {
    return interval<real>(real(-1), real(1));
  }

  if (width < real(M_PI / 2)) {
    const real derivative_lo = sin(i.lo());
    const real derivative_hi = sin(i.hi());

    const int saved_rounding_mode = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    const real a = derivative_lo < real(0) && derivative_hi > real(0)
                   ? real(-1)
                   : std::min(std::cos(i.lo()), std::cos(i.hi()));
    std::fesetround(FE_UPWARD);
    const real b = derivative_lo > real(0) && derivative_hi < real(0)
                   ? real(1)
                   : std::max(std::cos(i.lo()), std::cos(i.hi()));
    std::fesetround(saved_rounding_mode);

    return interval<real>(a, b);
  }

  return interval_union(i.lower_half(), i.upper_half());
}

template <typename real>
interval<real> sin(const interval<real>& i)
{
  const real width = i.hi() - i.lo();

  if (width >= real(2 * M_PI)) {
    return interval<real>(real(-1), real(1));
  }

  if (width < real(M_PI / 2)) {
    const real derivative_lo = -cos(i.lo());
    const real derivative_hi = -cos(i.hi());

    const int saved_rounding_mode = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    const real a = derivative_lo < real(0) && derivative_hi > real(0)
                   ? real(-1)
                   : std::min(std::sin(i.lo()), std::sin(i.hi()));
    std::fesetround(FE_UPWARD);
    const real b = derivative_lo > real(0) && derivative_hi < real(0)
                   ? real(1)
                   : std::max(std::sin(i.lo()), std::sin(i.hi()));
    std::fesetround(saved_rounding_mode);

    return interval<real>(a, b);
  }

  return interval_union(i.lower_half(), i.upper_half());
}

template <typename real>
interval<real> pow(const interval<real>& i, const interval<real>& j)
{
  if (possibly(i == real(0) && j == real(0))) {
    return interval<real>::construct_NaN();
  }

  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = min(min(pow(i.lo(), j.lo()), pow(i.hi(), j.lo())),
                     min(pow(i.lo(), j.hi()), pow(i.hi(), j.hi())));
  std::fesetround(FE_UPWARD);
  const real b = max(max(pow(i.lo(), j.lo()), pow(i.hi(), j.lo())),
                     max(pow(i.lo(), j.hi()), pow(i.hi(), j.hi())));
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
real foo(real x, real y)
{
  real result = atan2(x, y);
  std::cerr << "atan2(" << x << ", " << y << ") = " << result << std::endl;
  return result;
}

template <typename real>
interval<real> atan2(const interval<real>& i, const interval<real>& j)
{
  if (possibly(i == real(0) && j == real(0))) {
    return interval<real>::construct_NaN();
  }

  const int saved_rounding_mode = std::fegetround();
  std::fesetround(FE_DOWNWARD);
  const real a = min(min(foo(i.lo(), j.lo()), foo(i.hi(), j.lo())),
                     min(foo(i.lo(), j.hi()), foo(i.hi(), j.hi())));
  std::fesetround(FE_UPWARD);
  const real b = max(max(foo(i.lo(), j.lo()), foo(i.hi(), j.lo())),
                     max(foo(i.lo(), j.hi()), foo(i.hi(), j.hi())));
  std::fesetround(saved_rounding_mode);
  return interval<real>(a, b);
}

template <typename real>
interval<real> min(const interval<real>& i, const interval<real>& j)
{
  return interval<real>(min(i.lo(), j.lo()),
                        min(i.hi(), j.hi()));
}

template <typename real>
interval<real> min(const interval<real>& i, real j)
{
  return interval<real>(min(i.lo(), j),
                        min(i.hi(), j));
}

template <typename real>
interval<real> min(real i, const interval<real>& j)
{
  return interval<real>(min(i, j.lo()),
                        min(i, j.hi()));
}

template <typename real>
interval<real> max(const interval<real>& i, const interval<real>& j)
{
  return interval<real>(max(i.lo(), j.lo()),
                        max(i.hi(), j.hi()));
}

template <typename real>
interval<real> max(const interval<real>& i, real j)
{
  return interval<real>(max(i.lo(), j),
                        max(i.hi(), j));
}

template <typename real>
interval<real> max(real i, const interval<real>& j)
{
  return interval<real>(max(i, j.lo()),
                        max(i, j.hi()));
}

template <typename real>
interval<real> abs(const interval<real>& i)
{
  const real a = abs(i.lo());
  const real b = abs(i.hi());
  return interval<real>(min(a, b), max(a, b));
}

template <typename real>
bool isfinite(const interval<real>& i)
{
  return isfinite(i.lo()) && isfinite(i.hi());
}

template <typename real>
bool isnan(const interval<real>& i)
{
  return isnan(i.lo()) || isnan(i.hi());
}

template <typename real>
class numeric_limits<interval<real>>
  : public numeric_limits<real>
{
  typedef numeric_limits<real> base;
public:
  static const interval<real> epsilon() { return base::epsilon(); }
  static const interval<real> infinity() { return base::infinity(); }
  static const interval<real> quiet_NaN() { return base::quiet_NaN(); }
  static const interval<real> min() { return base::min(); }
  static const interval<real> max() { return base::max(); }
};

} // namespace std

#endif // interval_h_
