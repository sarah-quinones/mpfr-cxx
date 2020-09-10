#include "mpfr/mpfr.hpp"

#include <boost/multiprecision/mpfr.hpp>
#include "nanobench.h"

template <int N> using scalar_t = mpfr::mp_float_t<mpfr::digits10{N}>;
template <int N>
using bscalar_t = boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<N, boost::multiprecision::allocate_stack>,
    boost::multiprecision::et_off>;

template <typename T, typename Rhs>
void bench_mul_div(ankerl::nanobench::Bench& bench, Rhs x, std::string const& name) {
  T a = sqrt(T{2.0});
  T c{};

  ankerl::nanobench::doNotOptimizeAway(&a);
  ankerl::nanobench::doNotOptimizeAway(&c);

  bench.run("100 digits10: mul 128" + name, [&] {
    c = a * x;
    ankerl::nanobench::clobberMemory();
  });
  bench.run("100 digits10: div 128" + name, [&] {
    c = a / x;
    ankerl::nanobench::clobberMemory();
  });
}

auto main() -> int {

  auto bench = ankerl::nanobench::Bench();

  bench.minEpochTime(std::chrono::milliseconds{100UL});

  bench_mul_div<scalar_t<100>>(bench, 128.0, "d");
  bench_mul_div<scalar_t<100>>(bench, 128, "i");
  bench_mul_div<bscalar_t<100>>(bench, 128.0, "d boost");
  bench_mul_div<bscalar_t<100>>(bench, 128, "i boost");

  bench_mul_div<scalar_t<1000>>(bench, 128.0, "d");
  bench_mul_div<scalar_t<1000>>(bench, 128, "i");
  bench_mul_div<bscalar_t<1000>>(bench, 128.0, "d boost");
  bench_mul_div<bscalar_t<1000>>(bench, 128, "i boost");
}
