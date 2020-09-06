#ifndef MPFR_HPP_NZTOL31N
#define MPFR_HPP_NZTOL31N

#ifndef CXX_MPFR_SINGLE_HEADER
#include "mpfr/enums.hpp"
#endif

#if CXX_MPFR_DEBUG == 1
#define CXX_MPFR_STRINGIZE2(...) #__VA_ARGS__
#define CXX_MPFR_STRINGIZE(...) CXX_MPFR_STRINGIZE2(__VA_ARGS__)
#define CXX_MPFR_ASSERT(...)                                                                       \
  (static_cast<bool>(__VA_ARGS__)                                                                  \
       ? (void)0                                                                                   \
       : ::mpfr::_::crash_with_message("assertion_failed at " __FILE__                             \
                                       ":" CXX_MPFR_STRINGIZE(__LINE__) "\n" #__VA_ARGS__))
#else
#define CXX_MPFR_ASSERT(...) ((void)0)
#endif

namespace mpfr {

template <precision_t> struct mp_float_t;

namespace _ {

[[noreturn]] inline void crash_with_message(char const* message) {
  ::fprintf(stderr, "%s\n", message);
  std::terminate();
}

using std::size_t;
using std::uint64_t;

constexpr auto round_up_to_multiple(uint64_t n, uint64_t k) -> uint64_t {
  return (n % k == 0) ? n : (n / k * k + k);
}

constexpr auto digits10_to_2(uint64_t n) -> uint64_t {
  return static_cast<uint64_t>(
             3.321928094887362347870319429489390175864831393024580612054756395L *
             static_cast<long double>(n)) +
         1;
}

constexpr auto midpoint(uint64_t a, uint64_t b) -> uint64_t {
  return (a / 2 + b / 2 + ((a % 2) * (b % 2)));
}

constexpr auto cmp_inverse(uint64_t estimate, uint64_t target, uint64_t (*func)(uint64_t)) -> int {
  return (func(estimate) >= target and func(estimate - 1) < target)
             ? 0
             : (func(estimate) > target ? 1 : -1);
}

constexpr auto
inverse_binary_search(uint64_t n, uint64_t low, uint64_t high, uint64_t (*func)(uint64_t))
    -> uint64_t {
  return cmp_inverse(midpoint(low, high), n, func) == 0
             ? midpoint(low, high)
             : cmp_inverse(midpoint(low, high), n, func) == 1
                   ? inverse_binary_search(n, low, midpoint(low, high) - 1, func)
                   : inverse_binary_search(n, midpoint(low, high) + 1, high, func);
}

/// func must be strictly increasing
/// returns smallest m such that func(m) >= n
constexpr auto inverse(uint64_t n, uint64_t (*func)(uint64_t)) {
  return n == 0 ? 0 : inverse_binary_search(n, 0, n, func);
}

constexpr auto _sqr(uint64_t n) -> uint64_t { return n * n; }
constexpr auto sqrt(uint64_t n) -> uint64_t { return inverse(n, _sqr); }
constexpr auto digits2_to_10(uint64_t n) -> uint64_t { return inverse(n, digits10_to_2); }

constexpr auto prec_to_nlimb(mpfr_prec_t prec) -> uint64_t {
  return round_up_to_multiple(mpfr_custom_get_size(prec), sizeof(mp_limb_t)) / sizeof(mp_limb_t);
}

template <typename T> struct remove_pointer;
template <typename T> struct remove_pointer<T*> { using type = T; };

static constexpr uint64_t limb_pack_size = 32 / sizeof(mp_limb_t);
static constexpr mp_limb_t zero_pack[limb_pack_size] = {};
static constexpr mp_limb_t pow2_mantissa_last = mp_limb_t{1}
                                                << mp_limb_t{sizeof(mp_limb_t) * CHAR_BIT - 1};

inline auto is_zero_pack(mp_limb_t const* ptr) -> bool {
  return std::memcmp(ptr, zero_pack, sizeof(zero_pack)) == 0;
}

inline constexpr auto prec_negate_if(mpfr_prec_t p, bool cond) -> mpfr_prec_t {
  return cond ? static_cast<mpfr_prec_t>(~static_cast<mpfr_uprec_t>(p)) : p;
}

inline constexpr auto prec_abs(mpfr_prec_t p) -> mpfr_prec_t { return prec_negate_if(p, p < 0); }

inline auto compute_actual_prec(mpfr_srcptr x) -> mpfr_prec_t {

  if (mpfr_custom_get_kind(x) != MPFR_REGULAR_KIND) {
    return 0;
  }

  auto const* xp = static_cast<mp_limb_t const*>(mpfr_custom_get_mantissa(x));

  // size of mantissa minus last block
  size_t size = prec_to_nlimb(mpfr_get_prec(x)) - 1;

  // count leading mantissa zeros
  mpfr_prec_t zero_bits = 0;

  size_t head = size / limb_pack_size * limb_pack_size;
  size_t end = size;

  size_t i = 0;
  for (; i < head; i += limb_pack_size) {
    if (is_zero_pack(xp + i)) {
      zero_bits += sizeof(zero_pack) * CHAR_BIT;
    } else {
      end = i + limb_pack_size;
      break;
    }
  }

  for (; i < end; ++i) {
    if (xp[i] == 0) {
      zero_bits += sizeof(mp_limb_t) * CHAR_BIT;
    } else {
      break;
    }
  }

  // count trailing zeros
  mp_limb_t last_limb = xp[i];
  if (last_limb == 0) {
    zero_bits += sizeof(mp_limb_t) * CHAR_BIT;
  } else {
#if HEDLEY_HAS_BUILTIN(__builtin_ctzll)
    zero_bits += __builtin_ctzll(last_limb);
#else
    while ((last_limb % 2) == 0) {
      last_limb /= 2;
      ++zero_bits;
    }
#endif
  }
  return static_cast<mpfr_prec_t>(round_up_to_multiple(
             static_cast<uint64_t>(mpfr_get_prec(x)), CHAR_BIT * sizeof(mp_limb_t))) //
         - zero_bits;
}

struct mpfr_cref_t {
  typename remove_pointer<mpfr_ptr>::type m;

  void flip_sign_if(bool cond) { cond ? (void)(MPFR_SIGN(&m) *= -1) : (void)0; }
  [[nodiscard]] auto pow2_exponent() const -> mpfr_exp_t { return mpfr_custom_get_exp(&m) - 1; }
  void set_pow2_exponent(mpfr_exp_t e) { mpfr_custom_get_exp(&m) = e + 1; }
};

struct mpfr_raii_setter_t /* NOLINT */ {
  /* can only set once
   * sign and exponent of m must be set
   * if precision of m is equal to actual_precision, compute actual_precision
   * otherwise, set actual_precision_ptr's value to actual_precision
   */
  typename remove_pointer<mpfr_ptr>::type m;
  mpfr_prec_t m_actual_precision;
  mpfr_exp_t* m_exponent_ptr{};
  mpfr_prec_t* m_actual_prec_sign_ptr{};

  mpfr_raii_setter_t(
      mpfr_prec_t precision,
      mp_limb_t* mantissa,
      mpfr_exp_t* exponent_ptr,
      mpfr_prec_t* actual_prec_sign_ptr)
      : m{},
        m_actual_precision{precision},
        m_exponent_ptr{exponent_ptr},
        m_actual_prec_sign_ptr{actual_prec_sign_ptr} {

    mpfr_custom_init_set( //
        &m,
        MPFR_NAN_KIND,
        0,
        precision,
        mantissa);

    // poison value
    MPFR_SIGN(&m) = std::numeric_limits<mpfr_sign_t>::max();
    mpfr_custom_get_exp(&m) = std::numeric_limits<mpfr_exp_t>::max();
  }
  mpfr_raii_setter_t(mpfr_raii_setter_t const&) = delete;
  mpfr_raii_setter_t(mpfr_raii_setter_t&&) = delete;
  auto operator=(mpfr_raii_setter_t const&) -> mpfr_raii_setter_t& = delete;
  auto operator=(mpfr_raii_setter_t &&) -> mpfr_raii_setter_t& = delete;

  ~mpfr_raii_setter_t() {
    // if exponent is poisoned, do nothing
    if (mpfr_custom_get_exp(&m) == std::numeric_limits<mpfr_exp_t>::max()) {
      return;
    }
    // exponent is set
    // if sign is poisoned, set to positive
    if (MPFR_SIGN(&m) == std::numeric_limits<mpfr_sign_t>::max()) {
      MPFR_SIGN(&m) = 1;
    }
    CXX_MPFR_ASSERT(MPFR_SIGN(&m) != std::numeric_limits<mpfr_sign_t>::max());
    CXX_MPFR_ASSERT(mpfr_custom_get_exp(&m) != std::numeric_limits<mpfr_exp_t>::max());

    if (mpfr_regular_p(&m)) {
      *m_exponent_ptr = mpfr_custom_get_exp(&m);
      *m_actual_prec_sign_ptr = prec_negate_if(
          (m_actual_precision == mpfr_get_prec(&m)) ? compute_actual_prec(&m) : m_actual_precision,
          mpfr_signbit(&m));
    } else {
      std::memset(
          mpfr_custom_get_significand(&m), 0, sizeof(mp_limb_t) * prec_to_nlimb(mpfr_get_prec(&m)));
      *m_actual_prec_sign_ptr = prec_negate_if(0, mpfr_signbit(&m));
      *m_exponent_ptr = mpfr_zero_p(&m) ? 0 : mpfr_get_exp(&m);
    }
  }

  void flip_sign_if(bool cond) { cond ? (void)(MPFR_SIGN(&m) *= -1) : (void)0; }
  [[nodiscard]] auto pow2_exponent() const -> mpfr_exp_t { return mpfr_custom_get_exp(&m) - 1; }
  void set_pow2_exponent(mpfr_exp_t e) { mpfr_custom_get_exp(&m) = e + 1; }
};

struct mpfr_raii_inout_t /* NOLINT */ {
  /* can only set once
   * sign and exponent of m must be set
   * if precision of m is equal to actual_precision, compute actual_precision
   * otherwise, set actual_precision_ptr's value to actual_precision
   */
  typename remove_pointer<mpfr_ptr>::type m;
  mpfr_exp_t* m_exponent_ptr{};
  mpfr_prec_t* m_actual_prec_sign_ptr{};

  mpfr_raii_inout_t(
      mpfr_prec_t precision,
      mp_limb_t* mantissa,
      mpfr_exp_t* exponent_ptr,
      mpfr_prec_t* actual_prec_sign_ptr)
      : m{
        precision,
        (*actual_prec_sign_ptr) < 0 ? -1 : 1,
        (*exponent_ptr),
        mantissa,
      },
        m_exponent_ptr{exponent_ptr},
        m_actual_prec_sign_ptr{actual_prec_sign_ptr} {}

  mpfr_raii_inout_t(mpfr_raii_inout_t const&) = delete;
  mpfr_raii_inout_t(mpfr_raii_inout_t&&) = delete;
  auto operator=(mpfr_raii_inout_t const&) -> mpfr_raii_inout_t& = delete;
  auto operator=(mpfr_raii_inout_t &&) -> mpfr_raii_inout_t& = delete;

  ~mpfr_raii_inout_t() {
    mpfr_exp_t new_exponent{};
    mpfr_prec_t new_actual_prec_sign{};
    if (mpfr_regular_p(&m)) {
      new_exponent = mpfr_custom_get_exp(&m);
      new_actual_prec_sign = prec_negate_if(compute_actual_prec(&m), mpfr_signbit(&m));
    } else {
      new_actual_prec_sign = prec_negate_if(0, mpfr_signbit(&m));
      new_exponent = mpfr_zero_p(&m) ? 0 : mpfr_get_exp(&m);
      typename remove_pointer<mpfr_ptr>::type original{
          mpfr_get_prec(&m),
          (*m_actual_prec_sign_ptr) < 0 ? -1 : 1,
          (*m_exponent_ptr),
          nullptr,
      };
      if (mpfr_regular_p(&original)) {
        std::memset(
            mpfr_custom_get_significand(&m),
            0,
            sizeof(mp_limb_t) * prec_to_nlimb(mpfr_get_prec(&m)));
      }
    }

    if (*m_exponent_ptr != new_exponent) {
      *m_exponent_ptr = new_exponent;
    }
    if (*m_actual_prec_sign_ptr != new_actual_prec_sign) {
      *m_actual_prec_sign_ptr = new_actual_prec_sign;
    }
  }

  void flip_sign_if(bool cond) { cond ? (void)(MPFR_SIGN(&m) *= -1) : (void)0; }
  [[nodiscard]] auto pow2_exponent() const -> mpfr_exp_t { return mpfr_custom_get_exp(&m) - 1; }
  void set_pow2_exponent(mpfr_exp_t e) { mpfr_custom_get_exp(&m) = e + 1; }
};

[[nodiscard]] HEDLEY_ALWAYS_INLINE auto mul_b_is_pow2(
    mpfr_exp_t* a_exp,
    mpfr_exp_t* a_prec_sign,
    mpfr_exp_t b_exponent,
    mpfr_exp_t b_actual_prec_sign,
    bool div) -> bool {
  typename _::remove_pointer<mpfr_ptr>::type ea_{0, 0, *a_exp, nullptr};

  mpfr_prec_t eb{div ? (1 - b_exponent) : b_exponent - 1};
  if (HEDLEY_LIKELY(
          mpfr_regular_p(&ea_) and                //
          (*a_exp + eb - 1) < mpfr_get_emax() and //
          (*a_exp + eb - 1) > mpfr_get_emin())) {

    *a_prec_sign = (b_actual_prec_sign < 0) ? -*a_prec_sign : *a_prec_sign;
    *a_exp += eb;
    return true;
  }
  return false;
}

inline void set_d(mpfr_raii_setter_t& out, double a) {
  mpfr_set_d(&out.m, a, MPFR_RNDN);
  if (mpfr_get_prec(&out.m) >= static_cast<mpfr_prec_t>(sizeof(double) * CHAR_BIT)) {
    auto p = out.m;

    size_t full_n_limb = prec_to_nlimb(mpfr_get_prec(&out.m));
    size_t actual_n_limb = prec_to_nlimb(sizeof(double) * CHAR_BIT);

    mpfr_custom_move(
        &p,
        (static_cast<mp_limb_t*>(mpfr_custom_get_significand(&p)) + (full_n_limb - actual_n_limb)));
    mpfr_get_prec(&p) = sizeof(double) * CHAR_BIT;

    out.m_actual_precision = compute_actual_prec(&p);
  }
}

inline void set_add(mpfr_raii_setter_t& out, mpfr_cref_t a, mpfr_cref_t b) {
  mpfr_add(&out.m, &a.m, &b.m, MPFR_RNDN);
}

inline void set_sub(mpfr_raii_setter_t& out, mpfr_cref_t a, mpfr_cref_t b) {
  MPFR_SIGN(&b.m) *= -1;
  set_add(out, a, b);
}

inline void set_mul(mpfr_raii_setter_t& out, mpfr_cref_t a, mpfr_cref_t b) {
  if ((mpfr_custom_get_significand(&a.m) == mpfr_custom_get_significand(&b.m)) and
      (mpfr_custom_get_exp(&a.m) == mpfr_custom_get_exp(&b.m)) and
      (mpfr_signbit(&a.m) == mpfr_signbit(&b.m)) and (mpfr_get_prec(&a.m) == mpfr_get_prec(&b.m))) {
    mpfr_sqr(&out.m, &a.m, MPFR_RNDN);
    return;
  }
  mpfr_mul(&out.m, &a.m, &b.m, MPFR_RNDN);
}

inline void set_div(mpfr_raii_setter_t& out, mpfr_cref_t a, mpfr_cref_t b) {
  mpfr_div(&out.m, &a.m, &b.m, MPFR_RNDN);
}

struct heap_str_t /* NOLINT(cppcoreguidelines-special-member-functions) */ {
  char* p;
  explicit heap_str_t(size_t n) : p{n > 0 ? new char[n] : nullptr} {}

  ~heap_str_t() { delete[] p; }
};

template <typename CharT, typename Traits>
auto has_flag(
    std::basic_ostream<CharT, Traits>& out, typename std::basic_ostream<CharT, Traits>::fmtflags f)
    -> bool {
  return (out.flags() & f) == f;
}

template <typename CharT, typename Traits>
void print_n(std::basic_ostream<CharT, Traits>& out, CharT c, size_t n) {
  constexpr size_t pack_size = 64;
  CharT buffer[pack_size];
  for (auto& e : buffer) {
    e = c;
  }
  while (n >= pack_size) {
    n -= pack_size;
    out.write(buffer, pack_size);
  }
  out.write(buffer, static_cast<std::streamsize>(n));
}

template <typename CharT, typename Traits>
inline void write_to_ostream(
    std::basic_ostream<CharT, Traits>& out,
    mpfr_cref_t x_,
    char* stack_buffer,
    size_t stack_bufsize) {
  using ostr = std::basic_ostream<CharT, Traits>;

  char format[128] = {};
  std::size_t pos = 0;
  format[pos++] = '%';
  format[pos++] = '.';

  bool hf = _::has_flag(out, ostr::scientific) and _::has_flag(out, ostr::fixed);
  long precision = hf ? (mpfr_get_prec(&x_.m) / 4) : static_cast<long>(out.precision());
  std::snprintf(format + pos, sizeof(format) - pos, "%ld", precision);
  pos = std::strlen(format);
  format[pos++] = 'R';

  bool u = has_flag(out, ostr::uppercase);
  if (hf) {
    format[pos++] = u ? 'A' : 'a';
  } else if (_::has_flag(out, ostr::scientific)) {
    format[pos++] = u ? 'E' : 'e';
  } else if (_::has_flag(out, ostr::fixed)) {
    format[pos++] = u ? 'F' : 'f';
  } else {
    format[pos++] = u ? 'G' : 'g';
  }
  format[pos++] = '\0';

  bool signbit = mpfr_signbit(&x_.m);
  MPFR_SIGN(&x_.m) = 1;

  auto const size_needed = static_cast<std::size_t>(mpfr_snprintf(nullptr, 0, format, &x_.m)) + 1;

  std::size_t zero_padding = 0;

  bool use_heap = size_needed > stack_bufsize;
  _::heap_str_t heap_buffer{use_heap ? size_needed : 0};

  char* ptr = use_heap ? heap_buffer.p : stack_buffer;

  mpfr_snprintf(ptr, size_needed, format, &x_.m);
  if (not hf and _::has_flag(out, ostr::showpoint) and out.precision() > 0) {
    char* dot_ptr = std::strchr(ptr, '.');
    if (dot_ptr == nullptr) {
      zero_padding = 1 + static_cast<std::size_t>(out.precision()) - size_needed;
    }
  }

  std::size_t n_padding = 0;
  if (static_cast<std::size_t>(out.width()) >=
      size_needed - 1 + ((signbit or _::has_flag(out, ostr::showpos)) ? 1 : 0)) {

    n_padding = static_cast<std::size_t>(out.width()) - size_needed + 1 -
                ((signbit or _::has_flag(out, ostr::showpos)) ? 1 : 0) - zero_padding;
  }

  out.width(0);
  if (_::has_flag(out, ostr::right)) {
    _::print_n(out, out.fill(), n_padding);
  }

  if (signbit) {
    out.put(out.widen('-'));
  } else if (_::has_flag(out, ostr::showpos)) {
    out.put(out.widen('+'));
  }

  if (_::has_flag(out, ostr::internal)) {
    _::print_n(out, out.fill(), n_padding);
  }

  out << ptr;

  if (zero_padding > 0) {
    out.put(out.widen('.'));
    _::print_n(out, '0', zero_padding - 1);
  }

  if (_::has_flag(out, ostr::left)) {
    _::print_n(out, out.fill(), n_padding);
  }
}

struct impl_access {

  template <precision_t P>
  [[nodiscard]] static auto mantissa_mut(mp_float_t<P>& x)
      -> mp_limb_t (&)[prec_to_nlimb(static_cast<std::uint64_t>(P))] {
    return x.m_mantissa;
  }
  template <precision_t P>
  [[nodiscard]] static auto mantissa_const(mp_float_t<P> const& x) -> mp_limb_t
      const (&)[prec_to_nlimb(static_cast<std::uint64_t>(P))] {
    return x.m_mantissa;
  }

  template <precision_t P>
  [[nodiscard]] static auto actual_prec_sign_mut(mp_float_t<P>& x) -> mpfr_prec_t& {
    return x.m_actual_prec_sign;
  }
  template <precision_t P>
  [[nodiscard]] static auto actual_prec_sign_const(mp_float_t<P> const& x) -> mpfr_prec_t {
    return x.m_actual_prec_sign;
  }

  template <precision_t P>[[nodiscard]] static auto exp_mut(mp_float_t<P>& x) -> mpfr_exp_t& {
    return x.m_exponent;
  }
  template <precision_t P>
  [[nodiscard]] static auto exp_const(mp_float_t<P> const& x) -> mpfr_exp_t {
    return x.m_exponent;
  }

  template <precision_t P>
  [[nodiscard]] static auto mpfr_cref(mp_float_t<P> const& x) -> mpfr_cref_t {
    mpfr_cref_t out{};
    mpfr_sign_t sign = (x.m_actual_prec_sign < 0) ? -1 : 1;

    mpfr_prec_t actual_prec = prec_abs(x.m_actual_prec_sign);

    if (actual_prec == 0) {
      mpfr_custom_init_set(&out.m, sign * MPFR_ZERO_KIND, x.m_exponent, 1, x.m_mantissa);
      return out;
    }
    constexpr size_t full_n_limb = prec_to_nlimb(mp_float_t<P>::precision);
    size_t actual_n_limb = prec_to_nlimb(actual_prec);
    out.m = {
        actual_prec,
        sign,
        x.m_exponent,
        const_cast<mp_limb_t*> // NOLINT(cppcoreguidelines-pro-type-const-cast)
        (x.m_mantissa + (full_n_limb - actual_n_limb)),
    };
    return out;
  }

  template <precision_t P>
  [[nodiscard]] static auto mpfr_setter(mp_float_t<P>& x) -> mpfr_raii_setter_t {
    return {
        mp_float_t<P>::precision,
        static_cast<mp_limb_t*>(x.m_mantissa),
        &x.m_exponent,
        &x.m_actual_prec_sign,
    };
  }

  template <precision_t P>
  [[nodiscard]] static auto mpfr_inout_setter(mp_float_t<P>& x) -> mpfr_raii_inout_t {
    return {
        mp_float_t<P>::precision,
        static_cast<mp_limb_t*>(x.m_mantissa),
        &x.m_exponent,
        &x.m_actual_prec_sign,
    };
  }
};

template <precision_t P>
auto apply_unary_op(mp_float_t<P> const& x, int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)) noexcept
    -> mp_float_t<P> {
  mp_float_t<P> out;
  {
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
    op(&g.m, &x_.m, MPFR_RNDN);
  }
  return out;
}

template <precision_t P1, precision_t P2>
auto apply_binary_op(
    mp_float_t<P1> const& x,
    mp_float_t<P2> const& y,
    int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)) noexcept
    -> mp_float_t<(P1 >= P2) ? P1 : P2> {
  mp_float_t<(P1 >= P2) ? P1 : P2> out;
  {
    _::mpfr_raii_setter_t&& g = _::impl_access::mpfr_setter(out);
    _::mpfr_cref_t x_ = _::impl_access::mpfr_cref(x);
    _::mpfr_cref_t y_ = _::impl_access::mpfr_cref(y);
    op(&g.m, &x_.m, &y_.m, MPFR_RNDN);
  }
  return out;
}

template <typename T> auto to_lvalue(T&& arg) -> T& { return arg; }
template <parameter_type T> struct into_mpfr;

template <> struct into_mpfr<in> {
  static auto get_pointer(mpfr_cref_t&& p) -> mpfr_srcptr { return &p.m; }
  template <precision_t P> static auto get_mpfr(mp_float_t<P> const& x) -> mpfr_cref_t {
    return impl_access::mpfr_cref(x);
  }
};

template <> struct into_mpfr<out> {
  static auto get_pointer(mpfr_raii_setter_t&& p) -> mpfr_ptr { return &p.m; }
  template <precision_t P> static auto get_mpfr(mp_float_t<P>& x) -> mpfr_raii_setter_t {
    return {
        static_cast<mpfr_prec_t>(P),
        static_cast<mp_limb_t*>(impl_access::mantissa_mut(x)),
        &impl_access::exp_mut(x),
        &impl_access::actual_prec_sign_mut(x),
    };
  }
};

template <> struct into_mpfr<inout> {
  static auto get_pointer(mpfr_raii_inout_t&& p) -> mpfr_ptr { return &p.m; }
  template <precision_t P> static auto get_mpfr(mp_float_t<P>& x) -> mpfr_raii_inout_t {
    return {
        static_cast<mpfr_prec_t>(P),
        static_cast<mp_limb_t*>(impl_access::mantissa_mut(x)),
        &impl_access::exp_mut(x),
        &impl_access::actual_prec_sign_mut(x),
    };
  }
};

} // namespace _

constexpr digits2::digits2(precision_t prec) noexcept : m_value{static_cast<std::uint64_t>(prec)} {}
constexpr digits2::operator precision_t() const noexcept {
  return static_cast<precision_t>(m_value);
}

constexpr digits10::digits10(precision_t prec) noexcept
    : m_value{_::digits2_to_10(static_cast<std::uint64_t>(prec))} {}
constexpr digits10::operator precision_t() const noexcept {
  return static_cast<precision_t>(_::digits10_to_2(m_value));
}

} // namespace mpfr

#endif /* end of include guard MPFR_HPP_NZTOL31N */
