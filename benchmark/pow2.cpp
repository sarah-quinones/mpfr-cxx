// c++ -std=c++17 -march=haswell -O3 -Iinclude

#define CXX_MPFR_DEBUG 0
#include "mpfr/mpfr.hpp"

#include <boost/multiprecision/mpfr.hpp>
#include "nanobench.h"

template <int N> using scalar_t = mpfr::mp_float_t<mpfr::digits10{N}>;
template <int N>
using bscalar_t = boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<N, boost::multiprecision::allocate_stack>,
    boost::multiprecision::et_off>;

void f(scalar_t<100>& a) { a = 128.0; }

auto g(scalar_t<100>& a, scalar_t<100>& b) { return a * b; }

auto main() -> int {

  auto bench = ankerl::nanobench::Bench();
  bench.minEpochIterations(100000);
  {
    using T = scalar_t<100>;
    T a = sqrt(T{2.0});
    T c;

    ankerl::nanobench::detail::doNotOptimizeAway(&a);
    ankerl::nanobench::detail::doNotOptimizeAway(&c);

    bench.run("100 digits10: mul 128", [&] {
      c = a * 128.0;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("100 digits10: div 128", [&] {
      c = a / 128.0;
      ankerl::nanobench::clobberMemory();
    });
  }
  {
    using T = bscalar_t<100>;
    T a = sqrt(T{2.0});
    T c;

    ankerl::nanobench::detail::doNotOptimizeAway(&a);
    ankerl::nanobench::detail::doNotOptimizeAway(&c);

    bench.run("100 digits10: boost mul 128", [&] {
      c = a * 128.0;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("100 digits10: boost div 128", [&] {
      c = a / 128.0;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("100 digits10: boost mul 128 int", [&] {
      c = a * 128;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("100 digits10: boost div 128 int", [&] {
      c = a / 128;
      ankerl::nanobench::clobberMemory();
    });
  }

  bench.run("", [] {});
  bench.run("", [] {});

  {
    using T = scalar_t<1000>;
    T a = sqrt(T{2.0});
    T c;

    ankerl::nanobench::detail::doNotOptimizeAway(&a);
    ankerl::nanobench::detail::doNotOptimizeAway(&c);

    bench.run("1000 digits10: mul 128", [&] {
      c = a * 128.0;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("1000 digits10: div 128", [&] {
      c = a / 128.0;
      ankerl::nanobench::clobberMemory();
    });
  }
  {
    using T = bscalar_t<1000>;
    T a = sqrt(T{2.0});
    T c;

    ankerl::nanobench::detail::doNotOptimizeAway(&a);
    ankerl::nanobench::detail::doNotOptimizeAway(&c);

    bench.run("1000 digits10: boost mul 128", [&] {
      c = a * 128.0;
      ankerl::nanobench::clobberMemory();
    });
    bench.minEpochIterations(1000)
        .run(
            "1000 digits10: boost div 128",
            [&] {
              c = a / 128.0;
              ankerl::nanobench::clobberMemory();
            })
        .minEpochIterations(100000);
    bench.run("1000 digits10: boost mul 128 int", [&] {
      c = a * 128;
      ankerl::nanobench::clobberMemory();
    });
    bench.run("1000 digits10: boost div 128 int", [&] {
      c = a / 128;
      ankerl::nanobench::clobberMemory();
    });
  }
}
