#ifndef MP_FLOAT_HPP_KC35IAEF
#define MP_FLOAT_HPP_KC35IAEF

#include "mpfr/detail/handle_as_mpfr.hpp"
#include "mpfr/detail/prologue.hpp"

namespace mpfr {

/// Stack allocated fixed precision floating point.\n
/// Arithmetic and comparison operators follow IEEE 754 rules.
template <precision_t Precision> struct mp_float_t {
  static constexpr precision_t precision = Precision;
  static_assert(static_cast<mpfr_prec_t>(Precision) > 0, "precision must be positive.");

  /// Default/Zero initialization sets the number to positive zero.
  ///
  /// `mp_float_t<_> a;`
  ///
  /// `mp_float_t<_> a{};`
  mp_float_t() noexcept { std::memset(this, 0, sizeof(*this)); };

  /// \n
  template <typename T, _::enable_if_t<_::integral_or_floating_point<T>::value, void*> = nullptr>
  mp_float_t // NOLINT(hicpp-explicit-conversions)
      (T a) noexcept
      : mp_float_t() {
    *this = a;
  }

  /// \n
  template <typename T, _::enable_if_t<_::integral_or_floating_point<T>::value, void*> = nullptr>
  auto operator=(T const& a) noexcept -> mp_float_t& {
    _::integral_or_floating_point<T>::set(
        m_exponent,
        m_actual_prec_sign,
        precision_mpfr,
        m_mantissa,
        sizeof(m_mantissa) / sizeof(m_mantissa[0]),
        a);
    return *this;
  }

  /// \n
  [[MPFR_CXX_NODISCARD]] explicit operator long double() const noexcept {
    _::mpfr_cref_t m = _::impl_access::mpfr_cref(*this);
    return mpfr_get_ld(&m.m, _::get_rnd());
  }
  [[MPFR_CXX_NODISCARD]] explicit operator intmax_t() const noexcept {
    _::mpfr_cref_t m = _::impl_access::mpfr_cref(*this);
    return mpfr_get_sj(&m.m, _::get_rnd());
  }
  [[MPFR_CXX_NODISCARD]] explicit operator uintmax_t() const noexcept {
    _::mpfr_cref_t m = _::impl_access::mpfr_cref(*this);
    return mpfr_get_uj(&m.m, _::get_rnd());
  }
  /** @name Arithmetic operators
   */
  ///@{
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator+(mp_float_t const& a) noexcept -> mp_float_t {
    return a;
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator-(mp_float_t const& a) noexcept -> mp_float_t {
    mp_float_t out{a};
    out.m_actual_prec_sign = _::prec_negate_if(out.m_actual_prec_sign, true);
    return out;
  }

  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator+(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    return arithmetic_op(a, b, _::set_add);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator-(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    return arithmetic_op(a, b, _::set_sub);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator*(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    mp_float_t const* pow2 = nullptr;
    mp_float_t const* other; // NOLINT(cppcoreguidelines-init-variables)
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      pow2 = &b;
      other = &a;
    } else if (_::prec_abs(a.m_actual_prec_sign) == 1) {
      pow2 = &a;
      other = &b;
    }

    if (pow2 != nullptr) {
      mp_float_t out = *other;
      if (_::mul_b_is_pow2(
              &out.m_exponent,
              &out.m_actual_prec_sign,
              pow2->m_exponent,
              pow2->m_actual_prec_sign,
              false)) {
        return out;
      }
    }

    return arithmetic_op(a, b, _::set_mul);
  }

  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator/(mp_float_t const& a, mp_float_t const& b) noexcept
      -> mp_float_t {
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      mp_float_t out;
      out.m_exponent = a.m_exponent;
      out.m_actual_prec_sign = a.m_actual_prec_sign;
      if (_::mul_b_is_pow2(
              &out.m_exponent, &out.m_actual_prec_sign, b.m_exponent, b.m_actual_prec_sign, true)) {
        std::memcpy(out.m_mantissa, a.m_mantissa, sizeof(out.m_mantissa));
        return out;
      }
    }

    return arithmetic_op(a, b, _::set_div);
  }
  ///@}

  /** @name Assignment arithmetic operators
   */
  ///@{
  /// \n
  auto operator+=(mp_float_t const& b) noexcept -> mp_float_t& {
    *this = *this + b;
    return *this;
  }
  /// \n
  auto operator-=(mp_float_t const& b) noexcept -> mp_float_t& {
    *this = *this - b;
    return *this;
  }
  /// \n
  auto operator*=(mp_float_t const& b) noexcept -> mp_float_t& {
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      if (_::mul_b_is_pow2(
              &m_exponent, &m_actual_prec_sign, b.m_exponent, b.m_actual_prec_sign, false)) {
        return *this;
      }
    }
    *this = *this * b;
    return *this;
  }
  /// \n
  auto operator/=(mp_float_t const& b) noexcept -> mp_float_t& {
    if (_::prec_abs(b.m_actual_prec_sign) == 1) {
      if (_::mul_b_is_pow2( //
              &m_exponent,
              &m_actual_prec_sign,
              b.m_exponent,
              b.m_actual_prec_sign,
              true)) {
        return *this;
      }
    }
    *this = *this / b;
    return *this;
  }
  ///@}

  /** @name Comparison operators
   */
  ///@{
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator==(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_equal_p);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator!=(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_lessgreater_p);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator<(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_less_p);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator<=(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_lessequal_p);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator>(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_greater_p);
  }
  /// \n
  [[MPFR_CXX_NODISCARD]] friend auto operator>=(mp_float_t const& a, mp_float_t const& b) noexcept
      -> bool {
    return comparison_op(a, b, mpfr_greaterequal_p);
  }
  ///@}

  /// Write the number to an output stream.
  template <typename CharT, typename Traits>
  friend auto operator<<(std::basic_ostream<CharT, Traits>& out, mp_float_t const& a)
      -> std::basic_ostream<CharT, Traits>& {
    constexpr std::size_t stack_bufsize =
        _::digits2_to_10(static_cast<std::size_t>(precision_mpfr) + 64);
    char stack_buffer[stack_bufsize];
    _::write_to_ostream(out, _::impl_access::mpfr_cref(a), stack_buffer, stack_bufsize);
    return out;
  }

private:
  friend struct _::impl_access;

  [[MPFR_CXX_NODISCARD]] static auto arithmetic_op(
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

  [[MPFR_CXX_NODISCARD]] static auto
  comparison_op(mp_float_t const& a, mp_float_t const& b, int (*comp)(mpfr_srcptr, mpfr_srcptr))
      -> bool {
    _::mpfr_cref_t a_ = _::impl_access::mpfr_cref(a);
    _::mpfr_cref_t b_ = _::impl_access::mpfr_cref(b);
    return comp(&a_.m, &b_.m) != 0;
  }
  static constexpr mpfr_prec_t precision_mpfr = static_cast<mpfr_prec_t>(Precision);

  mp_limb_t m_mantissa[_::prec_to_nlimb(static_cast<std::uint64_t>(Precision))]{};
  mpfr_exp_t m_exponent{};
  mpfr_exp_t m_actual_prec_sign{};
}; // namespace mpfr

/// Allows handling `mp_float_t<_>` objects through `mpfr_ptr`/`mpfr_srcptr` proxy objects.
/// If after the callable is executed, one of the arguments has been modified by mpfr, the
/// corresponding `mp_float_t<_>` object is set to the equivalent value.
///
/// If any of the following occurs, the behavior is undefined:
/// * A parameter that is modified aliases another parameter.
/// * The `mp_float_t<_>` object referenced by a parameter that is modified is accessed
/// while the callable is being run (except through the `mpfr_ptr` proxy).
///
/// \par Side effects
/// Arguments are read before the callable is executed regardless of whether the
/// corresponding `mpfr_srcptr`/`mpfr_ptr` is accessed by the callable
//
/// \return The return value of the callable.
///
/// @param[in] fn   Arguments of type `mp_float_t<_>`.
/// @param[in] fn   A callable that takes arguments of type `mpfr_ptr` for mutable
/// parameters, or `mpfr_srcptr` for immutable parameters.
template <typename Fn, typename... Args>
callable_return_type
handle_as_mpfr_t(Fn&& fn, Args&&... args) // NOLINT(modernize-use-trailing-return-type)
    noexcept(callable_is_noexcept) {
  return _::impl_handle_as_mpfr_t<callable_is_noexcept>(static_cast<Fn&&>(fn), args...);
}

#undef callable_is_noexcept
#undef callable_return_type
} // namespace mpfr

#if defined(min) or defined(max)
HEDLEY_WARNING("min/max macros definitions are undone.")
#undef min
#undef max
#endif
namespace std {
/// Specialization of the standard library numeric limits
template <mpfr::precision_t Precision> struct numeric_limits<mpfr::mp_float_t<Precision>> {
  using T = mpfr::mp_float_t<Precision>;
  static constexpr bool is_specialized = true;

  /// Largest finite number.
  static auto max() noexcept -> T {
    T out = one_m_eps_impl();
    mpfr::_::impl_access::exp_mut(out) = mpfr_get_emax();
    return out;
  }
  /// Smallest strictly positive number.
  static auto min() noexcept -> T {
    T out = 0.5L;
    mpfr::_::impl_access::exp_mut(out) = mpfr_get_emin();
    return out;
  }
  /// Smallest finite number.
  static auto lowest() noexcept -> T { return -(max)(); }

  static constexpr int digits = static_cast<int>(Precision);
  static constexpr int digits10 =
      static_cast<int>(::mpfr::_::digits2_to_10(static_cast<uint64_t>(Precision)) + 1);
  static constexpr int max_digits10 = digits10;
  static constexpr bool is_signed = true;
  static constexpr bool is_integer = false;
  static constexpr bool is_exact = true;
  static constexpr int radix = 2;
  /// Distance between 1 and the smallest number larger than 1.
  static auto epsilon() noexcept -> T {
    static auto eps = epsilon_impl();
    return eps;
  }
  /// Largest possible error in ULP.
  static constexpr auto round_error() noexcept -> T {
    return T{mpfr::_::get_rnd() == MPFR_RNDN ? 0.5 : 1.0};
  }

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
  /// Positive infinity.
  static auto infinity() noexcept -> T { return T{1.0L} / T{0.0L}; }
  /// Not a number.
  static auto quiet_NaN() noexcept -> T { return T{0.0L} / T{0.0L}; }
  /// Not a number.
  static constexpr auto signaling_NaN() noexcept -> T { return quiet_NaN(); }
  static constexpr auto denorm_min() noexcept -> T { return T{}; }

  static constexpr bool is_iec559 = true;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = false;

  static constexpr bool traps = false;
  static constexpr bool tinyness_before = false;
  static constexpr float_round_style round_style = round_toward_zero;

private:
  static auto epsilon_impl() noexcept -> T {
    T x = 1;
    {
      mpfr::_::mpfr_raii_setter_t&& g = mpfr::_::impl_access::mpfr_setter(x);
      mpfr_nextabove(&g.m);
      mpfr_sub_ui(&g.m, &g.m, 1, mpfr::_::get_rnd());
    }
    return x;
  }
  static auto one_m_eps_impl() noexcept -> T {
    T x = 1;
    {
      mpfr::_::mpfr_raii_setter_t&& g = mpfr::_::impl_access::mpfr_setter(x);
      mpfr_nextbelow(&g.m);
    }
    return x;
  }
};
} // namespace std
#include "mpfr/detail/epilogue.hpp"

#endif /* end of include guard MP_FLOAT_HPP_KC35IAEF */
