#ifndef interval_ready_h_
#define interval_ready_h_

template <typename real> class interval;
class bool_interval;

namespace std {

template <typename real>
bool isnan(const interval<real>& i);

template <typename real>
interval<real> sqrt(const interval<real>& i);

template <typename real>
interval<real> exp(const interval<real>& i);

template <typename real>
interval<real> log(const interval<real>& i);

template <typename real>
interval<real> cos(const interval<real>& i);

template <typename real>
interval<real> sin(const interval<real>& i);

template <typename real>
interval<real> pow(const interval<real>& i, real j);

template <typename real>
interval<real> pow(real i, const interval<real>& j);

template <typename real>
interval<real> pow(const interval<real>& i, const interval<real>& j);

template <typename real>
interval<real> atan2(const interval<real>& i, const interval<real>& j);

template <typename real>
interval<real> min(const interval<real>& i, const interval<real>& j);

template <typename real>
interval<real> min(const interval<real>& i, real j);

template <typename real>
interval<real> min(real i, const interval<real>& j);

template <typename real>
interval<real> max(const interval<real>& i, const interval<real>& j);

template <typename real>
interval<real> max(const interval<real>& i, real j);

template <typename real>
interval<real> max(real i, const interval<real>& j);

template <typename real>
interval<real> abs(const interval<real>& i);

template <typename real>
bool isfinite(const interval<real>& i);

} // namespace std

inline bool certainly(bool b)
{
  return b;
}

inline bool possibly(bool b)
{
  return b;
}

#endif // interval_ready_h_
