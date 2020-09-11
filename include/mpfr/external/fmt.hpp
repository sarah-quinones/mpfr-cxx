/*
 Formatting library for C++

 Copyright (c) 2012 - present, Victor Zverovich

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 --- Optional exception to the license ---

 As an exception, if, as a result of your compiling your source code, portions
 of this Software are embedded into a machine-executable object form of such
 source code, you may redistribute such embedded portions in such object form
 without including the above copyright and permission notices.
 */

#ifndef FMT_HPP_8STTEZY1
#define FMT_HPP_8STTEZY1

#include "mpfr/detail/mpfr.hpp"
#include "mpfr/detail/prologue.hpp"

#include <fmt/format.h>

namespace mpfr {
namespace _ {

enum struct float_format_e : unsigned char {
  general,
  exp,
  fixed,
  hex,
  binary,
};

struct mp_float_specs {
  int precision;
  float_format_e format;
  fmt::sign_t sign;
  bool upper = false;
  bool showpoint = false;
};

FMT_NORETURN FMT_API inline void on_error(const char* message) {
  FMT_THROW(fmt::format_error{message});
}

template <typename Error_Handler>
FMT_CONSTEXPR auto
parse_mp_float_type_specs(fmt::basic_format_specs<char> const& specs, Error_Handler&& eh)
    -> mp_float_specs {
  auto result = mp_float_specs{};
  result.showpoint = specs.alt;

  switch (specs.type) {
  case 0:
    result.format = float_format_e::general;
    result.showpoint |= specs.precision > 0;
    break;
  case 'G':
    result.upper = true;
    FMT_FALLTHROUGH;
  case 'g':
    result.format = float_format_e::general;
    break;
  case 'E':
    result.upper = true;
    FMT_FALLTHROUGH;
  case 'e':
    result.format = float_format_e::exp;
    result.showpoint |= specs.precision != 0;
    break;
  case 'F':
    result.upper = true;
    FMT_FALLTHROUGH;
  case 'f':
    result.format = float_format_e::fixed;
    result.showpoint |= specs.precision != 0;
    break;
  case 'A':
    result.upper = true;
    FMT_FALLTHROUGH;
  case 'a':
    result.format = float_format_e::hex;
    break;
  case 'b':
    result.format = float_format_e::binary;
    break;
  default:
    eh.on_error("invalid type specifier for mp_float_t");
    break;
  }
  return result;
}

inline char* format_to_buf(
    char* stack_buffer,
    size_t stack_bufsize,
    heap_str_t& heap_buffer,
    long& precision,
    long& width,
    bool& signbit,
    size_t& zero_padding,
    size_t& n_padding,
    size_t& left_padding,
    size_t& right_padding,
    fmt::basic_format_specs<char> const& specs_,
    mp_float_specs const& fspecs,
    mpfr_cref_t x) {
  char format[128] = {};
  std::size_t pos = 0;
  format[pos++] = '%';
  format[pos++] = '.';

  bool hex = fspecs.format == mpfr::_::float_format_e::hex;
  bool bin = fspecs.format == mpfr::_::float_format_e::binary;
  precision = hex ? mpfr_get_prec(&x.m) / 4
                  : bin ? mpfr_get_prec(&x.m) : static_cast<long>(specs_.precision);
  width = static_cast<long>(specs_.width);
  signbit = mpfr_signbit(&x.m);

  if (precision > 0) {
    std::snprintf(format + pos, sizeof(format) - pos, "%ld", precision);
  }
  pos = std::strlen(format);
  format[pos++] = 'R';

  bool u = fspecs.upper;
  if (fspecs.format == mpfr::_::float_format_e::hex) {
    format[pos++] = u ? 'A' : 'a';
  } else if (fspecs.format == mpfr::_::float_format_e::binary) {
    format[pos++] = 'b';
  } else if (fspecs.format == mpfr::_::float_format_e::exp) {
    format[pos++] = u ? 'E' : 'e';
  } else if (fspecs.format == mpfr::_::float_format_e::fixed) {
    format[pos++] = u ? 'F' : 'f';
  } else {
    format[pos++] = u ? 'G' : 'g';
  }
  format[pos++] = '\0';

  MPFR_SIGN(&x.m) = 1;
  auto const size_needed = static_cast<std::size_t>(mpfr_snprintf(nullptr, 0, format, &x.m)) + 1;

  bool use_heap = size_needed > stack_bufsize;
  heap_buffer.init(use_heap ? size_needed : 0);

  char* ptr = use_heap ? heap_buffer.p : stack_buffer;

  mpfr_snprintf(ptr, size_needed, format, &x.m);
  if (not hex and fspecs.showpoint and precision > 0) {
    char* dot_ptr = std::strchr(ptr, '.');
    if (dot_ptr == nullptr) {
      zero_padding = 1 + static_cast<std::size_t>(precision) - size_needed;
    }
  }

  if (static_cast<std::size_t>(width) >=
      size_needed - 1 + ((signbit or specs_.sign == fmt::sign_t::plus) ? 1 : 0)) {

    n_padding = static_cast<std::size_t>(width) - size_needed + 1 -
                ((signbit or specs_.sign == fmt::sign_t::plus) ? 1 : 0) - zero_padding;
  }

  left_padding = (specs_.align == fmt::align_t::none or specs_.align == fmt::align_t::right) //
                     ? n_padding
                     : (specs_.align == fmt::align_t::center) //
                           ? (n_padding / 2)
                           : 0;
  right_padding = (specs_.align == fmt::align_t::left) //
                      ? n_padding
                      : (specs_.align == fmt::align_t::center) //
                            ? (n_padding - left_padding)
                            : 0;
  return ptr;
}

} // namespace _
} // namespace mpfr

namespace fmt {

template <mpfr::precision_t P> struct formatter<mpfr::mp_float_t<P>> {
  FMT_CONSTEXPR formatter() = default;
  using T = mpfr::mp_float_t<P>;
  using Char = char;

  template <typename ParseContext>
  FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin()) {
    using handler_type = detail::dynamic_specs_handler<ParseContext>;
    auto type = detail::type_constant<T, Char>::value;
    detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);

    auto it = parse_format_specs(ctx.begin(), ctx.end(), handler);
    auto eh = ctx.error_handler();
    mpfr::_::parse_mp_float_type_specs(specs_, eh);
    return it;
  }

  template <typename Format_Context>
  auto format(mpfr::mp_float_t<P> const& value, Format_Context& ctx) -> decltype(ctx.out()) {

    constexpr std::size_t stack_bufsize =
        mpfr::_::digits2_to_10(static_cast<std::size_t>(mpfr::mp_float_t<P>::precision) + 64);
    char stack_buffer[stack_bufsize];

    auto eh = ctx.error_handler();

    detail::handle_dynamic_spec<detail::width_checker>(specs_.width, specs_.width_ref, ctx);
    detail::handle_dynamic_spec<detail::precision_checker>(
        specs_.precision, specs_.precision_ref, ctx);
    mpfr::_::mp_float_specs fspecs = mpfr::_::parse_mp_float_type_specs(specs_, eh);

    long precision = 0;
    long width = 0;
    std::size_t zero_padding = 0;
    std::size_t n_padding = 0;
    bool signbit = false;
    size_t left_padding = 0;
    size_t right_padding = 0;

    mpfr::_::heap_str_t heap_buffer{0};
    auto ptr = mpfr::_::format_to_buf(
        stack_buffer,
        stack_bufsize,
        heap_buffer,
        precision,
        width,
        signbit,
        zero_padding,
        n_padding,
        left_padding,
        right_padding,
        specs_,
        fspecs,
        mpfr::_::impl_access::mpfr_cref(value));

    char const* fmt_format_parts[] = {
        "{0:",
        "<{1}}{2}{0:",
        "<{3}}{4}{5}{0:0<{6}}{0:",
        "<{7}}",
    };
    std::size_t sizes[] = {
        std::strlen(fmt_format_parts[0]),
        std::strlen(fmt_format_parts[1]),
        std::strlen(fmt_format_parts[2]),
        std::strlen(fmt_format_parts[3]),
    };

    char const* fill_begin = specs_.fill.data();
    std::size_t fill_size = specs_.fill.size();
    char fmt_format[256]{};

    char* fmt_ptr = fmt_format;

    for (int i = 0; i < 4; ++i) {
      std::memcpy(fmt_ptr, fmt_format_parts[i], sizes[i]);
      fmt_ptr += sizes[i];

      if (i < 3) {
        if (fill_begin != nullptr) {
          std::memcpy(fmt_ptr, fill_begin, fill_size);
          fmt_ptr += fill_size;
        }
      } else {
        *fmt_ptr = '\0';
      }
    }

    return fmt::format_to(
        ctx.out(),
        fmt_format,
        "",
        left_padding,
        signbit ? "-" : specs_.sign == sign_t::plus ? "+" : "",
        specs_.align == align_t::numeric ? n_padding : 0,
        ptr,
        zero_padding > 0 ? "." : "",
        zero_padding > 0 ? (zero_padding - 1) : 0,
        right_padding);
  }

  detail::dynamic_format_specs<char> specs_;
};

} // namespace fmt

#include "mpfr/detail/epilogue.hpp"

#endif /* end of include guard FMT_HPP_8STTEZY1 */
