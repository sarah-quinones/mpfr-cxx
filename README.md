lightweight wrapper around [gnu mpfr](https://www.mpfr.org/). allowing usage of stack allocated multiprecision floats.  
[documentation](https://mpfr-cxx.readthedocs.io/en/latest/)

advantages over boost::multiprecision::mpfr
  - scalar types are [trivially copyable](https://en.cppreference.com/w/cpp/named_req/TriviallyCopyable).
  - a possible value representation of positive zero is all zero bits. i.e., for an object `x` of type `mp_float_t<_>`, `std::memset(&x, 0, sizeof(x))`, will set `x` to `+0.0`.
  - optimized for certain common operations (multiplication and division by powers of two).
  - better `ostream` formatting support.
  - smaller compilation overhead.
  - [fmtlib](https://github.com/fmtlib/fmt) support (requires v7). the format specification is interpreted as [standard format specification](https://en.cppreference.com/w/cpp/utility/format/formatter#Standard_format_specification), with the addition of a binary type specifier `b`.
 
# usage
put `./include/mpfr` in one of the compiler's include directories, or run `make_single_header.py` to generate a single header include file.

```cpp
#include "mpfr/mpfr.hpp"
#include <iomanip>
#include <iostream>
#include <fmt/core.h>

using mpfr::digits10;
using mpfr::digits2;
using mpfr::mp_float_t;
using scalar_t = mp_float_t<digits10{48}>; // same as mp_float_t<digits2{160}>

auto euler_mascheroni_constant() -> scalar_t {
  // advanced usage with mpfr
  scalar_t y;

  // binds a mpfr_ptr proxy (yp) to y, and sets y to the corresponding value
  mpfr::handle_as_mpfr_t(
      [](mpfr_ptr yp) {
        mpfr_const_euler(yp, MPFR_RNDN); // lim n->∞ sum(i in [1.. n]) 1/i
      },
      y);
  return y;
}

auto main() -> int {
  scalar_t x = 0;

  scalar_t n = 1'000'000;
  for (long i = 1; i <= n; ++i) {
    x += 1 / scalar_t{i};
  }

  auto y = euler_mascheroni_constant();

  // ostream support
  std::cout << std::setprecision(30) << std::left;
  std::cout << std::setw(50) << (x - (y + log(n) + 1 / (2 * n) - 1 / (12 * n * n))) << '\n';

  // optional fmtlib support
  fmt::print("{:^50.{}e}\n", x, 20);
  fmt::print("{:^50.{}A}\n", x, 20);

  // compile time format parsing, unicode fill character, binary format
  fmt::print("{:💖<{}b}\n", x, 200);
}
```

compiled with: `clang++ example.cpp -Iinclude -lfmt -lmpfr`
