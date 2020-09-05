#define CXX_MPFR_DEBUG 1
#include "doctest.h"
#include "mpfr/mpfr.hpp"
#include <iostream>

using namespace mpfr;
using scalar_t = mp_float_t<digits10{1000}>;

DOCTEST_TEST_CASE("in") {
  scalar_t x = 2;
  handle_as_mpfr_t<in>(
      [](mpfr_srcptr a) {
        mpfr_t b;
        mpfr_t c;

        mpfr_init2(b, 64);
        mpfr_init2(c, 64);

        mpfr_set_ui(b, 4, MPFR_RNDN);
        mpfr_add(c, a, a, MPFR_RNDN);

        DOCTEST_CHECK(mpfr_equal_p(c, b) != 0);

        mpfr_clear(b);
        mpfr_clear(c);
      },
      x);
}

DOCTEST_TEST_CASE("inout") {
  scalar_t x = 2;
  scalar_t y = 7;
  handle_as_mpfr_t<inout, in>(
      [](mpfr_ptr a, mpfr_srcptr b) {
        mpfr_t c;

        mpfr_init2(c, 64);

        mpfr_set_ui(c, 7, MPFR_RNDN);

        DOCTEST_CHECK(mpfr_equal_p(b, c) != 0);

        mpfr_add(a, a, c, MPFR_RNDN);
        mpfr_mul(a, a, b, MPFR_RNDN);

        mpfr_clear(c);
      },
      x,
      y);
  DOCTEST_CHECK(x == 63);
}

DOCTEST_TEST_CASE("out") {
  scalar_t x;
  scalar_t y = 7;
  handle_as_mpfr_t<out, in>([](mpfr_ptr a, mpfr_srcptr b) { mpfr_exp(a, b, MPFR_RNDN); }, x, y);
  DOCTEST_CHECK(x == exp(y));
}

DOCTEST_TEST_CASE("out-other") {
  using scalar512_t = mp_float_t<digits2{512}>;
  using scalar1024_t = mp_float_t<digits2{1024}>;
  scalar512_t x; // uninitialized
  scalar1024_t y = scalar1024_t{2.0 + 1.0 / 1024};
  handle_as_mpfr_t<out, in>(
      [](mpfr_ptr a, mpfr_srcptr b) {
        mpfr_set_str(a, "0.1111", 2, MPFR_RNDN); // set to 0b0.1111 == 0.9375
        mpfr_add(a, a, b, MPFR_RNDN);
      },
      x,
      y);
  DOCTEST_CHECK(x == scalar512_t{2 + 0.5 + 0.25 + 0.125 + 0.0625 + 1.0 / 1024});
}
