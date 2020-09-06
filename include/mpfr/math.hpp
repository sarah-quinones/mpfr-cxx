#ifndef MATH_HPP_IED2CNIL
#define MATH_HPP_IED2CNIL

#include "mpfr/mp_float.hpp"
#include "mpfr/detail/prologue.hpp"

namespace mpfr {

/// \return `true` if the argument is negative, `false` otherwise.\n
/// If the argument is a NaN, detects whether the sign bit is set.
template <precision_t P> auto signbit(mp_float_t<P> const& arg) noexcept -> bool {
  return _::impl_access::actual_prec_sign_const(arg) < 0;
}

/// \return `true` if the argument is infinite, `false` otherwise.
template <precision_t P> auto isinf(mp_float_t<P> const& arg) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(arg);
  return mpfr_inf_p(&x_.m);
}

/// \return `true` if the argument is finite, `false` otherwise.
template <precision_t P> auto isfinite(mp_float_t<P> const& arg) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(arg);
  return mpfr_number_p(&x_.m);
}

/// \return `true` if the argument is a NaN, `false`, otherwise.
template <precision_t P> auto isnan(mp_float_t<P> const& arg) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(arg);
  return mpfr_nan_p(&x_.m);
}

/// \return The absolute value of the argument.
template <precision_t P> auto fabs(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = arg;
  _::impl_access::actual_prec_sign_mut(out) =
      _::prec_abs(_::impl_access::actual_prec_sign_const(arg));
  return out;
}

/// \return The absolute value of the argument.
template <precision_t P> auto abs(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return mpfr::fabs(arg);
}

/// \return The base to the power of the exponent.
template <precision_t P>
auto pow(mp_float_t<P> const& base, mp_float_t<P> const& exponent) noexcept -> mp_float_t<P> {
  return _::apply_binary_op(base, exponent, mpfr_pow);
}

/// If `sin_dest` and `cos_dest` alias, the behavior is undefined.
/// @param[in]  arg      Argument at which the sine and cosine are to be computed.
/// @param[out] sin_dest Reference to the number that will be set to the sine value.
/// @param[out] cos_dest Reference to the number that will be set to the cosine value.
template <precision_t P>
void sincos(mp_float_t<P> const& arg, mp_float_t<P>& sin_dest, mp_float_t<P>& cos_dest) noexcept {
  if (sin_dest == cos_dest) {
    _::crash_with_message("sin_dest and cos_dest cannot alias.");
  }
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(arg);
  _::mpfr_raii_setter_t sg = _::impl_access::mpfr_setter(sin_dest);
  _::mpfr_raii_setter_t cg = _::impl_access::mpfr_setter(cos_dest);
  mpfr_sin_cos(&sg.m, &cg.m, &x_.m, MPFR_RNDN);
}

/// If `sinh_dest` and `cosh_dest` alias, the behavior is undefined.
/// @param[in]  arg       Argument at which the hyperbolic sine and cosine are to be computed.
/// @param[out] sinh_dest Reference to the number that will be set to the hyperbolic sine value.
/// @param[out] cosh_dest Reference to the number that will be set to the hyperbolic cosine value.
template <precision_t P>
void sinhcosh(
    mp_float_t<P> const& arg, mp_float_t<P>& sinh_dest, mp_float_t<P>& cosh_dest) noexcept {
  if (&sinh_dest == &cosh_dest) {
    _::crash_with_message("sinh_dest and cosh_dest cannot alias.");
  }
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(arg);
  _::mpfr_raii_setter_t sg = _::impl_access::mpfr_setter(sinh_dest);
  _::mpfr_raii_setter_t cg = _::impl_access::mpfr_setter(cosh_dest);
  mpfr_sinh_cosh(&sg.m, &cg.m, &x_.m, MPFR_RNDN);
}

/// \return Arc tangent of y/x in the correct quadrant depending on the signs of the arguments.
template <precision_t P>
auto atan2(mp_float_t<P> const& y, mp_float_t<P> const& x) noexcept -> mp_float_t<P> {
  mp_float_t<P> out;
  {
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
    _::mpfr_cref_t y_ = _::impl_access::mpfr_cref(y);
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);

    mpfr_atan2(&g.m, &y_.m, &x_.m, MPFR_RNDN);
  }
  return out;
}

/// \return Square root of the sum of the squares of the arguments.
template <precision_t P>
auto hypot(mp_float_t<P> const& x, mp_float_t<P> const& y) noexcept -> mp_float_t<P> {
  mp_float_t<P> out;
  {
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
    _::mpfr_cref_t y_ = _::impl_access::mpfr_cref(y);
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);

    mpfr_hypot(&g.m, &x_.m, &y_.m, MPFR_RNDN);
  }
  return out;
}

/// \return The next representable number of `from` in the direction of `to`.
template <precision_t P>
auto nexttoward(mp_float_t<P> const& from, mp_float_t<P> const& to) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = from;
  {
    _::mpfr_raii_inout_t&& g = _::impl_access::mpfr_inout_setter(from);
    _::mpfr_cref_t y_ = _::impl_access::mpfr_cref(to);
    mpfr_nexttoward(&g.m, &y_.m);
  }
  return out;
}

/// \return The next representable number of `from` in the direction of \f$+\infty\f$.
template <precision_t P> auto nextabove(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = arg;
  {
    _::mpfr_raii_inout_t&& g = _::impl_access::mpfr_inout_setter(arg);
    mpfr_nextabove(&g.m);
  }
  return out;
}

/// \return The next representable number of `from` in the direction of \f$-\infty\f$.
template <precision_t P> auto nextbelow(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = arg;
  {
    _::mpfr_raii_inout_t&& g = _::impl_access::mpfr_inout_setter(arg);
    mpfr_nextbelow(&g.m);
  }
  return out;
}

/// \return Square root of the argument.
template <precision_t P> auto sqrt(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_sqrt);
}

/// \return Cubic root of the argument.
template <precision_t P> auto cbrt(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_cbrt);
}

/// \return Sine of the argument.
template <precision_t P> auto sin(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_sin);
}
/// \return Cosine of the argument.
template <precision_t P> auto cos(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_cos);
}
/// \return Tangent of the argument.
template <precision_t P> auto tan(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_tan);
}
/// \return Inverse of the sine of the argument.
template <precision_t P> auto asin(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_asin);
}
/// \return Inverse of the cosine of the argument.
template <precision_t P> auto acos(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_acos);
}
/// \return Inverse of the tangent of the argument.
template <precision_t P> auto atan(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_atan);
}

/// \return Hyperbolic sine of the argument.
template <precision_t P> auto sinh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_sinh);
}
/// \return Hyperbolic cosine of the argument.
template <precision_t P> auto cosh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_cosh);
}
/// \return Hyperbolic tangent of the argument.
template <precision_t P> auto tanh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_tanh);
}
/// \return Inverse of the hyperbolic sine of the argument.
template <precision_t P> auto asinh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_asinh);
}
/// \return Inverse of the hyperbolic cosine of the argument.
template <precision_t P> auto acosh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_acosh);
}
/// \return Inverse of the hyperbolic tangent of the argument.
template <precision_t P> auto atanh(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_atanh);
}

/// \return Exponential of the argument with base \f$e := e^{\text{arg}}\f$.
template <precision_t P> auto exp(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_exp);
}
/// \return Exponential of the argument with base \f$2 := 2^{\text{arg}}\f$.
template <precision_t P> auto exp2(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_exp2);
}
/// \return Exponential of the argument with base \f$10 := 10^{\text{arg}}\f$.
template <precision_t P> auto exp10(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_exp10);
}
/// \return Exponential of the argument with base \f$e\f$ minus 1 \f$:= e^{\text{arg}} - 1\f$.
template <precision_t P> auto expm1(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_expm1);
}
/// \return Logarithm of the argument to base \f$e := \log_e(\text{arg})\f$.
template <precision_t P> auto log(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_log);
}
/// \return Logarithm of the argument to base \f$2 := \log_2(\text{arg})\f$.
template <precision_t P> auto log2(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_log2);
}
/// \return Logarithm of the argument to base \f$10 := \log_10(\text{arg})\f$.
template <precision_t P> auto log10(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_log10);
}
/// \return Logarithm to base \f$e\f$ of \f$1\f$ plus the argument \f$:=\log_e(1 + \text{arg})\f$.
template <precision_t P> auto log1p(mp_float_t<P> const& arg) noexcept -> mp_float_t<P> {
  return _::apply_unary_op(arg, mpfr_log1p);
}

} // namespace mpfr

#include "mpfr/detail/epilogue.hpp"

#endif /* end of include guard MATH_HPP_IED2CNIL */
