#ifndef ENUMS_HPP_4PXFRLJ0
#define ENUMS_HPP_4PXFRLJ0

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <cstdint>
#include <iosfwd>
#include <utility>

#include <mpfr.h>

#ifndef CXX_MPFR_DEBUG
#define CXX_MPFR_DEBUG 0
#endif

#ifndef CXX_MPFR_SINGLE_HEADER
#include "mpfr/detail/hedley.h"
#endif

namespace mpfr {

/// The count of bits that make up the number's mantissa.\n
/// This must be a strictly positive value.
enum struct precision_t : mpfr_prec_t { reserved = 0 };

/// Types of the parameters passed to the mpfr interface layer.
///
/// See the documentation of :cpp:func:`mpfr::handle_as_mpfr_t` for details on the constraints of
/// each type.
enum parameter_type {
  /// Read-only argument.
  in,
  /// Writable argument.
  out,
  /// Readable/writable argument.
  inout,
};

/// Number of digits in base 2.
struct digits2 {
private:
  std::uint64_t m_value;

public:
  /// Constructs with the given number of digits.
  constexpr explicit digits2(std::uint64_t value) noexcept : m_value{value} {}
  /// Converts from precision enum.
  constexpr explicit digits2(precision_t prec) noexcept;
  /// Converts to precision enum.
  constexpr operator // NOLINT(hicpp-explicit-conversions)
      precision_t() const noexcept;
};

/// Number of digits in base 10.\n
/// The actual precision may be slightly larger since the conversion between base 2 and base 10 is
/// inexact.
struct digits10 {
private:
  std::uint64_t m_value;

public:
  /// Constructs with the given number of digits.
  constexpr explicit digits10(std::uint64_t value) noexcept : m_value{value} {}
  /// Converts from precision enum.
  constexpr explicit digits10(precision_t prec) noexcept;
  /// Converts to precision enum.
  constexpr operator // NOLINT(hicpp-explicit-conversions)
      precision_t() const noexcept;
};
} // namespace mpfr

#endif /* end of include guard ENUMS_HPP_4PXFRLJ0 */
