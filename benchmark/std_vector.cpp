#define CXX_MPFR_DEBUG 0
#include "mpfr/mpfr.hpp"

#include <vector>
#include <boost/multiprecision/mpfr.hpp>
#include "nanobench.h"

using scalar_t = mpfr::mp_float_t<mpfr::digits10{100}>;
using bscalar_t = boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<1000, boost::multiprecision::allocate_stack>,
    boost::multiprecision::et_off>;

template <typename T> T constant = sqrt(T{2});

template <typename T> void bench_create() {
  std::vector<T> v(1000);
  ankerl::nanobench::doNotOptimizeAway(v.data());
  ankerl::nanobench::clobberMemory();
}

template <typename T> void bench_resize() {
  std::vector<T> v;
  v.resize(1000);
  ankerl::nanobench::doNotOptimizeAway(v.data());
  ankerl::nanobench::clobberMemory();
}

template <typename T> void bench_resize_set() {
  std::vector<T> v;
  v.resize(1000, constant<T>);
  ankerl::nanobench::doNotOptimizeAway(v.data());
  ankerl::nanobench::clobberMemory();
}

template <typename T> void bench_grow() {
  std::vector<T> v;
  v.resize(1024, constant<T>);
  v.resize(v.size() * 2, constant<T>);
  v.resize(v.size() * 2, constant<T>);
  v.resize(v.size() * 2, constant<T>);
  ankerl::nanobench::doNotOptimizeAway(v.data());
  ankerl::nanobench::clobberMemory();
}

auto main() -> int {

  auto bench = ankerl::nanobench::Bench();
  bench.minEpochTime(std::chrono::nanoseconds{100'000'000ULL});

  bench.run("create", bench_create<scalar_t>);
  bench.run("create (boost)", bench_create<bscalar_t>);

  bench.run("resize", bench_resize<scalar_t>);
  bench.run("resize (boost)", bench_resize<bscalar_t>);

  bench.run("resize set", bench_resize_set<scalar_t>);
  bench.run("resize set (boost)", bench_resize_set<bscalar_t>);

  bench.run("grow", bench_grow<scalar_t>);
  bench.run("grow (boost)", bench_grow<bscalar_t>);
}
