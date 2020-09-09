#define MPFR_CXX_DEBUG 1

#include "doctest.h"
#include <iostream>
#include <fmt/format.h>
#include "mpfr/mpfr.hpp"
#include <cassert>

using namespace mpfr;
using scalar_t = mp_float_t<digits10{100}>;

DOCTEST_TEST_CASE("multiplication/division test") {
  scalar_t x = 1 - std::numeric_limits<scalar_t>::epsilon();
  scalar_t y = 0.5;
  scalar_t z = std::numeric_limits<scalar_t>::quiet_NaN();

  handle_as_mpfr_t(
      [](mpfr_ptr a, mpfr_ptr b) {
        mpfr_mul_2si(a, a, mpfr_get_emax(), MPFR_RNDN);
        mpfr_mul_2si(b, b, mpfr_get_emin(), MPFR_RNDN);
      },
      x,
      y);

  {
    DOCTEST_CHECK(x < std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(x * 2 == std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(2 * x == std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK((x / 2) + (x / 2) == (x / 2) * 2);
    DOCTEST_CHECK((x / 2) + (x / 2) == 2 * (x / 2));
    DOCTEST_CHECK(y > 0);
    DOCTEST_CHECK(y / 2 == 0);
    DOCTEST_CHECK(y + y == 2 * y);
    DOCTEST_CHECK(y + y == y * 2);
  }
  {
    DOCTEST_CHECK(-x > -std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(-x * 2 == -std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK((-2) * x == -std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(-(2 * x) == -std::numeric_limits<scalar_t>::infinity());

    DOCTEST_CHECK((-x / 2) + (-x / 2) == ((-x) / 2) * 2);
    DOCTEST_CHECK((-x / 2) + (-x / 2) == (-(x / 2)) * 2);
    DOCTEST_CHECK((-x / 2) + (-x / 2) == -((x / 2)) * 2);
    DOCTEST_CHECK((-x / 2) + (-x / 2) == 2 * ((-x) / 2));
    DOCTEST_CHECK((-x / 2) + (-x / 2) == 2 * (-(x / 2)));
    DOCTEST_CHECK((-x / 2) + (-x / 2) == 2 * -((x / 2)));

    DOCTEST_CHECK(-y < 0);
    DOCTEST_CHECK(-y / 2 == 0);
    DOCTEST_CHECK(-y + -y == 2 * -y);
    DOCTEST_CHECK(-y + -y == -y * 2);
  }
  {
    DOCTEST_CHECK(isnan(z));
    DOCTEST_CHECK(isnan(z / 2));
    DOCTEST_CHECK(isnan(z * 2));
  }
}

DOCTEST_TEST_CASE("math functions") {
  scalar_t const x = 1.312;

  DOCTEST_CHECK(fabs(cos(x) * cos(x) + sin(x) * sin(x) - 1) < 1e-99);
  DOCTEST_CHECK(fabs(cosh(x) * cosh(x) - sinh(x) * sinh(x) - 1) < 1e-99);
}
