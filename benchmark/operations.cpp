#include "mpfr/mpfr.hpp"

#include <boost/multiprecision/mpfr.hpp>
#include "nanobench.h"

template <int N> using scalar_t = mpfr::mp_float_t<mpfr::digits10{N}>;

auto main() -> int {

  auto bench = ankerl::nanobench::Bench();

  bench.minEpochTime(std::chrono::milliseconds{100UL});

  using T1 = scalar_t<512>;
  using T2 = scalar_t<1024>;

  T1 a = sqrt(T1{2.0});
  T2 b = sqrt(T2{2.0});
  auto a_ = a + 1;
  auto b_ = b + 1;

  T1 c;
  T2 d;

  ankerl::nanobench::doNotOptimizeAway(&a);
  ankerl::nanobench::doNotOptimizeAway(&b);
  ankerl::nanobench::doNotOptimizeAway(&c);
  ankerl::nanobench::doNotOptimizeAway(&d);

  bench.run("copying same size 512", [&] {
    c = a;
    ankerl::nanobench::clobberMemory();
  });
  bench.run("copying same size 1024", [&] {
    d = b;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("copying 512 -> 1024", [&] {
    d = a;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("copying 1024 -> 512", [&] {
    c = b;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("add 512", [&] {
    c = a + a;
    ankerl::nanobench::clobberMemory();
  });

  a_ = 2;
  bench.run("add2 512", [&] {
    c = a_ + a_;
    ankerl::nanobench::clobberMemory();
  });
  a_ = a + 1;

  bench.run("sqr 512", [&] {
    c = a * a;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("mul 512", [&] {
    c = a * a_;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("div 512", [&] {
    c = a / a_;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("mul2 512", [&] {
    c = a * 2;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("add 1024", [&] {
    d = b + b;
    ankerl::nanobench::clobberMemory();
  });

  b_ = 2;
  bench.run("add2 1024", [&] {
    d = b_ + b_;
    ankerl::nanobench::clobberMemory();
  });
  b_ = b + 1;

  bench.run("sqr 1024", [&] {
    d = b * b;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("mul 1024", [&] {
    d = b * b_;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("div 1024", [&] {
    d = b / b_;
    ankerl::nanobench::clobberMemory();
  });

  bench.run("mul2 1024", [&] {
    d = b * 2;
    ankerl::nanobench::clobberMemory();
  });
}
