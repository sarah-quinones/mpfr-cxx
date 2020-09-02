#ifndef MPFR_HPP_IKTJCQIU
#define MPFR_HPP_IKTJCQIU

#ifndef CXX_MPFR_SINGLE_HEADER
#include "mpfr/detail/mpfr.hpp"
#endif

namespace mpfr {
template <precision_t Precision> struct mp_float_t {
private:
  friend struct _::impl_access;

  mp_limb_t m_mantissa[_::prec_to_nlimb(static_cast<std::uint64_t>(Precision))];
  mpfr_exp_t m_exponent;
  mpfr_exp_t m_actual_prec_sign;

  [[nodiscard]] static auto arithmetic_op(
      mp_float_t const& a,
      mp_float_t const& b,
      void (*op)(_::mpfr_raii_setter_t&, _::mpfr_cref_t, _::mpfr_cref_t)) -> mp_float_t {
    mp_float_t out;
    {
      _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
      _::mpfr_cref_t a_ = _::impl_access::mpfr_cref(a);
      _::mpfr_cref_t b_ = _::impl_access::mpfr_cref(b);
      op(g, a_, b_);
    }
    return out;
  }

  [[nodiscard]] static auto
  comparison_op(mp_float_t const& a, mp_float_t const& b, int (*comp)(mpfr_srcptr, mpfr_srcptr))
      -> bool {
    _::mpfr_cref_t a_ = _::impl_access::mpfr_cref(a);
    _::mpfr_cref_t b_ = _::impl_access::mpfr_cref(b);
    return comp(&a_.m, &b_.m) != 0;
  }

  HEDLEY_ALWAYS_INLINE static auto
  mul_b_is_pow2(mp_float_t& a, mpfr_exp_t b_exponent, mpfr_exp_t b_actual_prec_sign, bool div)
      -> bool {
    _::mpfr_pretend_t ea_{a.m_exponent};
    mpfr_prec_t ea{a.m_exponent};
    mpfr_prec_t eb{div ? (1 - b_exponent) : b_exponent - 1};
    if (HEDLEY_LIKELY(
            mpfr_regular_p(&ea_) and            //
            (ea + eb - 1) < mpfr_get_emax() and //
            (ea + eb - 1) > mpfr_get_emin())) {

      a.m_actual_prec_sign *= (b_actual_prec_sign < 0) ? -1 : 1;
      a.m_exponent += eb;
      return true;
    }
    return false;
  }

public:
  static constexpr mpfr_prec_t precision = static_cast<mpfr_prec_t>(Precision);
  mp_float_t() noexcept = default;
  mp_float_t // NOLINT(hicpp-explicit-conversions,cppcoreguidelines-pro-type-member-init)
             // mantissa is set by assignment
      (double num) noexcept {
    *this = num;
  }
  inline auto operator=(double num) noexcept -> mp_float_t& {
    int exponent{};
    double normalized = ::fabs(::frexp(num, &exponent));
    if (normalized == 0.5) {
      bool signbit = std::signbit(num);

      // a is a power of two
      // a = signbit * 2^(exp-1)

      auto* xp = static_cast<mp_limb_t*>(m_mantissa);
      constexpr uint64_t size = sizeof(m_mantissa) / sizeof(m_mantissa[0]);

      if (size >= 1) {
        std::memset(xp, 0, sizeof(mp_limb_t) * (size - 1));
      }
      m_exponent = exponent;
      xp[size - 1] = _::pow2_mantissa_last;
      m_actual_prec_sign = signbit ? -1 : 1;

      return *this;
    } else {
      _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(*this);
      _::set_d(g, num);
    }
    return *this;
  }
  [[nodiscard]] explicit inline operator double() const noexcept {
    _::mpfr_cref_t m = _::impl_access::mpfr_cref(*this);
    return mpfr_get_d(&m.m, MPFR_RNDN);
  }

  [[nodiscard]] inline auto operator+() const noexcept -> mp_float_t { return *this; }
  [[nodiscard]] inline auto operator-() const noexcept -> mp_float_t {
    mp_float_t out{*this};
    out.m_actual_prec_sign = _::prec_negate_if(out.m_actual_prec_sign, true);
    return out;
  }

  [[nodiscard]] friend auto operator+(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    return arithmetic_op(a, b, _::set_add);
  }
  [[nodiscard]] friend auto operator-(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    return arithmetic_op(a, b, _::set_sub);
  }
  [[nodiscard]] friend auto operator*(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    mp_float_t const* pow2 = nullptr;
    mp_float_t const* other; // NOLINT
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      pow2 = &b;
      other = &a;
    } else if (_::prec_abs(a.m_actual_prec_sign) == 1) {
      pow2 = &a;
      other = &b;
    }

    if (pow2 != nullptr) {
      mp_float_t out = *other;
      if (mul_b_is_pow2(out, pow2->m_exponent, pow2->m_actual_prec_sign, false)) {
        return out;
      }
    }

    return arithmetic_op(a, b, _::set_mul);
  }

  [[nodiscard]] friend auto operator/(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      mp_float_t out = a;
      if (mul_b_is_pow2(out, b.m_exponent, b.m_actual_prec_sign, true)) {
        return out;
      }
    }

    return arithmetic_op(a, b, _::set_div);
  }

  [[nodiscard]] inline auto operator+=(mp_float_t const& other) noexcept -> mp_float_t& {
    *this = *this + other;
    return *this;
  }
  [[nodiscard]] inline auto operator-=(mp_float_t const& other) noexcept -> mp_float_t& {
    *this = *this - other;
    return *this;
  }
  [[nodiscard]] inline auto operator*=(mp_float_t const& other) noexcept -> mp_float_t& {
    if (_::prec_abs(other.m_actual_prec_sign) == 1) {
      if (mul_b_is_pow2(*this, other.m_exponent, other.m_actual_prec_sign, false)) {
        return *this;
      }
    } else {
      *this = *this * other;
    }
    return *this;
  }
  [[nodiscard]] inline auto operator/=(mp_float_t const& other) noexcept -> mp_float_t& {
    if (_::prec_abs(other.m_actual_prec_sign) == 1) {
      mul_b_is_pow2(*this, other, true);
    } else {
      *this = *this / other;
    }
    return *this;
  }

  friend auto operator==(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_equal_p);
  }
  friend auto operator!=(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_lessgreater_p);
  }
  friend auto operator<(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_less_p);
  }
  friend auto operator<=(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_lessequal_p);
  }
  friend auto operator>(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_greater_p);
  }
  friend auto operator>=(mp_float_t const& a, mp_float_t const& b) noexcept -> bool {
    return comparison_op(a, b, mpfr_greaterequal_p);
  }

  template <typename CharT, typename Traits>
  friend auto operator<<(std::basic_ostream<CharT, Traits>& out, mp_float_t const& a)
      -> std::basic_ostream<CharT, Traits>& {
    constexpr std::size_t stack_bufsize =
        _::digits2_to_10(static_cast<std::size_t>(precision) + 64);
    char stack_buffer[stack_bufsize];
    _::write_to_ostream(out, _::impl_access::mpfr_cref(a), stack_buffer, stack_bufsize);
    return out;
  }
  static auto from_string(char const* str, char** end_ptr = nullptr) -> mp_float_t {
    mp_float_t out;
    {
      _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
      mpfr_strtofr(&g.m, str, end_ptr, 0, MPFR_RNDN);
    }
    return out;
  }
};

// math functions

template <precision_t P> auto signbit(mp_float_t<P> const& x) noexcept -> bool {
  return _::impl_access::actual_prec_sign_const(x) < 0;
}

template <precision_t P> auto isinf(mp_float_t<P> const& x) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
  return mpfr_inf_p(&x_.m);
}

template <precision_t P> auto isfinite(mp_float_t<P> const& x) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
  return mpfr_number_p(&x_.m);
}

template <precision_t P> auto isnan(mp_float_t<P> const& x) noexcept -> bool {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
  return mpfr_nan_p(&x_.m);
}

template <precision_t P> auto fabs(mp_float_t<P> const& x) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = x;
  _::impl_access::actual_prec_sign_mut(out) =
      _::prec_abs(_::impl_access::actual_prec_sign_const(x));
  return out;
}

template <precision_t P>
auto pow(mp_float_t<P> const& x, mp_float_t<P> const& y) noexcept -> mp_float_t<P> {
  return _::apply_binary_op(x, y, mpfr_pow);
}

template <precision_t P>
void sincos(mp_float_t<P> const& x, mp_float_t<P>* s, mp_float_t<P>* c) noexcept {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
  _::mpfr_raii_setter_t sg = _::impl_access::mpfr_setter(*s);
  _::mpfr_raii_setter_t cg = _::impl_access::mpfr_setter(*c);
  mpfr_sin_cos(&sg.m, &cg.m, &x_.m, MPFR_RNDN);
}

template <precision_t P>
void sinhcosh(mp_float_t<P> const& x, mp_float_t<P>* s, mp_float_t<P>* c) noexcept {
  _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
  _::mpfr_raii_setter_t sg = _::impl_access::mpfr_setter(*s);
  _::mpfr_raii_setter_t cg = _::impl_access::mpfr_setter(*c);
  mpfr_sinh_cosh(&sg.m, &cg.m, &x_.m, MPFR_RNDN);
}

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

template <precision_t P>
auto nexttoward(mp_float_t<P> const& x, mp_float_t<P> const& y) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = x;
  {
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
    _::mpfr_cref_t y_ = _::impl_access::mpfr_cref(y);

    mpfr_custom_get_exp(&g.m) = mpfr_custom_get_exp(&x_.m);
    MPFR_SIGN(&g.m) = MPFR_SIGN(&x_.m);

    mpfr_nexttoward(&g.m, &y_.m);
  }
  return out;
}

template <precision_t P> auto nextabove(mp_float_t<P> const& x) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = x;
  {
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);

    mpfr_custom_get_exp(&g.m) = mpfr_custom_get_exp(&x_.m);
    MPFR_SIGN(&g.m) = MPFR_SIGN(&x_.m);

    mpfr_nextabove(&g.m);
  }
  return out;
}

template <precision_t P> auto nextbelow(mp_float_t<P> const& x) noexcept -> mp_float_t<P> {
  mp_float_t<P> out = x;
  {
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);

    mpfr_custom_get_exp(&g.m) = mpfr_custom_get_exp(&x_.m);
    MPFR_SIGN(&g.m) = MPFR_SIGN(&x_.m);

    mpfr_nextbelow(&g.m);
  }
  return out;
}

#define CXX_MPFR_UNARY_OP(Name)                                                                    \
  template <precision_t P> auto Name(mp_float_t<P> const& x) noexcept->mp_float_t<P> {             \
    return _::apply_unary_op(x, mpfr_##Name);                                                      \
  }                                                                                                \
  static_assert(true, "")

CXX_MPFR_UNARY_OP(sqrt);
CXX_MPFR_UNARY_OP(cbrt);

CXX_MPFR_UNARY_OP(sin);
CXX_MPFR_UNARY_OP(cos);
CXX_MPFR_UNARY_OP(tan);
CXX_MPFR_UNARY_OP(asin);
CXX_MPFR_UNARY_OP(acos);
CXX_MPFR_UNARY_OP(atan);

CXX_MPFR_UNARY_OP(sinh);
CXX_MPFR_UNARY_OP(cosh);
CXX_MPFR_UNARY_OP(tanh);
CXX_MPFR_UNARY_OP(asinh);
CXX_MPFR_UNARY_OP(acosh);
CXX_MPFR_UNARY_OP(atanh);

CXX_MPFR_UNARY_OP(exp);
CXX_MPFR_UNARY_OP(exp2);
CXX_MPFR_UNARY_OP(exp10);
CXX_MPFR_UNARY_OP(expm1);
CXX_MPFR_UNARY_OP(log);
CXX_MPFR_UNARY_OP(log2);
CXX_MPFR_UNARY_OP(log10);
CXX_MPFR_UNARY_OP(log1p);

} // namespace mpfr

// stl numeric_limits
namespace std {
template <mpfr::precision_t Precision> struct numeric_limits<mpfr::mp_float_t<Precision>> {
  using T = mpfr::mp_float_t<Precision>;
  static constexpr bool is_specialized = true;
  static auto(max)() noexcept -> T {
    T out = 0.5L;
    out._m_exponent = mpfr_get_emax();
  }
  static auto(min)() noexcept -> T {
    T out = 0.5L;
    out._m_exponent = mpfr_get_emin();
  }
  static auto lowest() noexcept -> T { return -(max)(); }

  static constexpr int digits = static_cast<int>(Precision);
  static constexpr int digits10 =
      static_cast<int>(::mpfr::_::digits2_to_10(static_cast<uint64_t>(Precision)) + 1);
  static constexpr int max_digits10 = digits10;
  static constexpr bool is_signed = true;
  static constexpr bool is_integer = false;
  static constexpr bool is_exact = true;
  static constexpr int radix = 2;
  static auto epsilon() noexcept -> T {
    static auto eps = mpfr::nextabove(T{1.0});
    return eps;
  }
  static constexpr auto round_error() noexcept -> T { return T{0.5L}; }

  static constexpr int min_exponent = (MPFR_EMIN_DEFAULT >= INT_MIN) ? MPFR_EMIN_DEFAULT : INT_MIN;
  static constexpr int min_exponent10 =
      -static_cast<int>(::mpfr::_::digits2_to_10(static_cast<uint64_t>(-min_exponent)) + 1);
  static constexpr int max_exponent = (MPFR_EMAX_DEFAULT <= INT_MAX) ? MPFR_EMAX_DEFAULT : INT_MAX;
  static constexpr int max_exponent10 =
      static_cast<int>(::mpfr::_::digits2_to_10(static_cast<uint64_t>(max_exponent)) + 1);

  static constexpr bool has_infinity = false;
  static constexpr bool has_quiet_NaN = true;
  static constexpr bool has_signaling_NaN = true;
  static constexpr float_denorm_style has_denorm = denorm_absent;
  static constexpr bool has_denorm_loss = false;
  static auto infinity() noexcept -> T { return T{1.0L} / T{0.0L}; }
  static auto quiet_NaN() noexcept -> T { return T{0.0L} / T{0.0L}; }
  static constexpr auto signaling_NaN() noexcept -> T { return quiet_NaN(); }
  static constexpr auto denorm_min() noexcept -> T { return T{}; }

  static constexpr bool is_iec559 = true;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = false;

  static constexpr bool traps = false;
  static constexpr bool tinyness_before = false;
  static constexpr float_round_style round_style = round_toward_zero;
};
} // namespace std

#endif /* end of include guard MPFR_HPP_IKTJCQIU */
