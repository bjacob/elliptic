#include <limits>
#include <algorithm>
#include <cmath>
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
           isnan(lo) ||
           isnan(hi));
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

  real diameter() const
  {
    return hi() - lo();
  }

  real radius() const
  {
    return real(0.5) * diameter();
  }

  void fuzz()
  {
    static const real one_plus_epsilon  = real(1) + std::numeric_limits<real>::epsilon();
    static const real one_minus_epsilon = real(1) - std::numeric_limits<real>::epsilon();

    if (m_hi > 0) {
      m_hi *= one_plus_epsilon;
    } else if (m_hi < 0) {
      m_hi *= one_minus_epsilon;
    }

    if (m_lo < 0) {
      m_lo *= one_plus_epsilon;
    } else if (m_lo > 0) {
      m_lo *= one_minus_epsilon;
    }
  }

  template <typename other_type>
  interval& operator+=(const other_type& other)
  {
    return *this = *this + other;
  }

  template <typename other_type>
  interval& operator-=(const other_type& other)
  {
    return *this = *this - other;
  }

  template <typename other_type>
  interval& operator*=(const other_type& other)
  {
    return *this = *this * other;
  }

  template <typename other_type>
  interval& operator/=(const other_type& other)
  {
    return *this = *this / other;
  }

  interval operator-() const
  {
    return interval<real>(-hi(), -lo());
  }
};

template <typename real>
interval<real> operator+(const interval<real>& i, const interval<real>& j)
{
  return interval<real>(i.lo() + j.lo(), i.hi() + j.hi());
}

template <typename real>
interval<real> operator-(const interval<real>& i, const interval<real>& j)
{
  return interval<real>(i.lo() - j.hi(), i.hi() - j.lo());
}

template <typename real>
interval<real> operator*(const interval<real>& i, const interval<real>& j)
{
  const real a = i.lo() * j.lo();
  const real b = i.hi() * j.lo();
  const real c = i.lo() * j.hi();
  const real d = i.hi() * j.hi();
  return interval<real>(std::min(std::min(a, b), std::min(c, d)),
                        std::max(std::max(a, b), std::max(c, d)));
}

template <typename real>
interval<real> operator/(const interval<real>& i, const interval<real>& j)
{
  const real a = i.lo() / j.lo();
  const real b = i.hi() / j.lo();
  const real c = i.lo() / j.hi();
  const real d = i.hi() / j.hi();
  return interval<real>(std::min(std::min(a, b), std::min(c, d)),
                        std::max(std::max(a, b), std::max(c, d)));
}

template <typename real>
interval<real> operator+(const interval<real>& i, real j)
{
  return interval<real>(i.lo() + j, i.hi() + j);
}

template <typename real>
interval<real> operator-(const interval<real>& i, real j)
{
  return interval<real>(i.lo() - j, i.hi() - j);
}

template <typename real>
interval<real> operator*(const interval<real>& i, real j)
{
  const real a = i.lo() * j;
  const real b = i.hi() * j;
  return interval<real>(std::min(a, b), std::max(a, b));
}

template <typename real>
interval<real> operator/(const interval<real>& i, real j)
{
  const real a = i.lo() / j;
  const real b = i.hi() / j;
  const real c = i.lo() / j;
  const real d = i.hi() / j;
  return interval<real>(std::min(a, b), std::max(a, b));

}

template <typename real>
interval<real> operator+(real i, const interval<real>& j)
{
  return interval<real>(i + j.lo(), i + j.hi());
}

template <typename real>
interval<real> operator-(real i, const interval<real>& j)
{
  return interval<real>(i - j.hi(), i - j.lo());
}

template <typename real>
interval<real> operator*(real i, const interval<real>& j)
{
  const real a = i * j.lo();
  const real b = i * j.hi();
  return interval<real>(std::min(a, b), std::max(a, b));
}

template <typename real>
interval<real> operator/(real i, const interval<real>& j)
{
  const real a = i / j.lo();
  const real b = i / j.hi();
  return interval<real>(std::min(a, b), std::max(a, b));
}



template <typename real>
interval<real> operator+(const interval<real>& i, int j)
{
  return interval<real>(i.lo() + j, i.hi() + j);
}

template <typename real>
interval<real> operator-(const interval<real>& i, int j)
{
  return interval<real>(i.lo() - j, i.hi() - j);
}

template <typename real>
interval<real> operator*(const interval<real>& i, int j)
{
  const real a = i.lo() * j;
  const real b = i.hi() * j;
  return interval<real>(std::min(a, b), std::max(a, b));
}

template <typename real>
interval<real> operator/(const interval<real>& i, int j)
{
  const real a = i.lo() / j;
  const real b = i.hi() / j;
  return interval<real>(std::min(a, b), std::max(a, b));

}

template <typename real>
interval<real> operator+(int i, const interval<real>& j)
{
  return interval<real>(i + j.lo(), i + j.hi());
}

template <typename real>
interval<real> operator-(int i, const interval<real>& j)
{
  return interval<real>(i - j.hi(), i - j.lo());
}

template <typename real>
interval<real> operator*(int i, const interval<real>& j)
{
  const real a = i * j.lo();
  const real b = i * j.hi();
  return interval<real>(std::min(a, b), std::max(a, b));
}

template <typename real>
interval<real> operator/(int i, const interval<real>& j)
{
  const real a = i / j.lo();
  const real b = i / j.hi();
  return interval<real>(std::min(a, b), std::max(a, b));
}

namespace std {

template <typename real>
interval<real> sqrt(const interval<real>& i)
{
  return interval<real>(sqrt(i.lo()), sqrt(i.hi()));
}

template <typename real>
interval<real> exp(const interval<real>& i)
{
  return interval<real>(exp(i.lo()), exp(i.hi()));
}

template <typename real>
interval<real> log(const interval<real>& i)
{
  return interval<real>(log(i.lo()), log(i.hi()));
}

template <typename real>
interval<real> cos(const interval<real>& i)
{
  const real a = cos(i.lo());
  const real b = cos(i.hi());
  return interval<real>(min(a, b), max(a, b));
}

template <typename real>
interval<real> sin(const interval<real>& i)
{
  const real a = sin(i.lo());
  const real b = sin(i.hi());
  return interval<real>(min(a, b), max(a, b));
}

template <typename real>
interval<real> pow(const interval<real>& i, real j)
{
  const real a = pow(i.lo(), j);
  const real b = pow(i.hi(), j);
  return interval<real>(min(a, b), max(a, b));
}

template <typename real>
interval<real> pow(real i, const interval<real>& j)
{
  const real a = pow(i, j.lo());
  const real b = pow(i, j.hi());
  return interval<real>(min(a, b), max(a, b));
}

template <typename real>
interval<real> pow(const interval<real>& i, const interval<real>& j)
{
  const real a = pow(i.lo(), j.lo());
  const real b = pow(i.hi(), j.lo());
  const real c = pow(i.lo(), j.hi());
  const real d = pow(i.hi(), j.hi());
  return interval<real>(min(min(a, b), min(c, d)),
                        max(max(a, b), max(c, d)));
}

template <typename real>
interval<real> atan2(const interval<real>& i, const interval<real>& j)
{
  const real a = atan2(i.lo(), j.lo());
  const real b = atan2(i.hi(), j.lo());
  const real c = atan2(i.lo(), j.hi());
  const real d = atan2(i.hi(), j.hi());
  return interval<real>(min(min(a, b), min(c, d)),
                        max(max(a, b), max(c, d)));
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
class numeric_limits<interval<real> >
  : public numeric_limits<real>
{
  typedef numeric_limits<real> base;
public:
  static const interval<real> infinity() { return base::infinity(); }
  static const interval<real> quiet_NaN() { return base::quiet_NaN(); }
  static const interval<real> min() { return base::min(); }
  static const interval<real> max() { return base::max(); }
};

} // namespace std

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

inline bool certainly(bool b)
{
  return b;
}

inline bool possibly(const bool_interval& b)
{
  return b.hi();
}

inline bool possibly(bool b)
{
  return b;
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
std::ostream& operator<<(std::ostream& s, const interval<real>& i)
{
  return s << "[ " << i.lo() << " .. " << i.hi() << " ]";
}

template <typename real>
bool_interval operator>=(const interval<real>& a, const interval<real>& b)
{
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
  return a == interval<real>(0);
}
