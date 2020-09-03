#include "doctest.h"
#include "mpfr/mpfr.hpp"
#include <iostream>

#define CXX_MPFR_DEBUG 1

using namespace mpfr;
using scalar_t = mp_float_t<digits10{1000}>;

DOCTEST_TEST_CASE("in") {
  scalar_t x = 2;
  apply_mpfr_fn<in>( //
      x,
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
      });
}

DOCTEST_TEST_CASE("inout") {
  scalar_t x = 2;
  scalar_t y = 7;
  apply_mpfr_fn<inout, in>( //
      x,
      y,
      [](mpfr_ptr a, mpfr_srcptr b) {
        mpfr_t c;

        mpfr_init2(c, 64);

        mpfr_set_ui(c, 7, MPFR_RNDN);

        DOCTEST_CHECK(mpfr_equal_p(b, c) != 0);

        mpfr_add(a, a, c, MPFR_RNDN);
        mpfr_mul(a, a, b, MPFR_RNDN);

        mpfr_clear(c);
      });
  DOCTEST_CHECK(x == 63);
}

DOCTEST_TEST_CASE("out") {
  scalar_t x;
  scalar_t y = 7;
  apply_mpfr_fn<out, in>( //
      x,
      y,
      [](mpfr_ptr a, mpfr_srcptr b) {
        mpfr_exp(a, b, MPFR_RNDN);
      });
  DOCTEST_CHECK(x == exp(y));
}
