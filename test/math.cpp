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
    DOCTEST_CHECK(isfinite(x));
    DOCTEST_CHECK(x * 2 == std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(2 * x == std::numeric_limits<scalar_t>::infinity());
    DOCTEST_CHECK(isinf(x * 2));
    DOCTEST_CHECK(isinf(2 * x));
    DOCTEST_CHECK(not signbit(x * 2));
    DOCTEST_CHECK(not signbit(2 * x));
    DOCTEST_CHECK(not isfinite(x * 2));
    DOCTEST_CHECK(not isfinite(2 * x));

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
    DOCTEST_CHECK(isinf((-2) * x));
    DOCTEST_CHECK(isinf(-(2 * x)));
    DOCTEST_CHECK(signbit((-2) * x));
    DOCTEST_CHECK(signbit(-(2 * x)));
    DOCTEST_CHECK(not isfinite((-2) * x));
    DOCTEST_CHECK(not isfinite(-(2 * x)));

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

    DOCTEST_CHECK(not isnan(x));
    DOCTEST_CHECK(not isnan(x / 2));
    DOCTEST_CHECK(not isnan(x * 2));

    DOCTEST_CHECK(not isnan(y));
    DOCTEST_CHECK(not isnan(y / 2));
    DOCTEST_CHECK(not isnan(y * 2));
  }
}

DOCTEST_TEST_CASE("math functions") {
  scalar_t const x = 1.312;
  scalar_t const y = 21.21922;

  DOCTEST_CHECK(fabs(x) == x);
  DOCTEST_CHECK(fabs(-x) == x);

  DOCTEST_CHECK(fabs(sin(mpfr::pi_c<scalar_t::precision>())) < 1e-99);
  DOCTEST_CHECK(fabs(cos(mpfr::pi_c<scalar_t::precision>()) + 1) < 1e-99);

  DOCTEST_CHECK(fabs(cos(x) * cos(x) + sin(x) * sin(x) - 1) < 1e-99);
  DOCTEST_CHECK(fabs(sin(x) / cos(x) - tan(x)) < 1e-99);
  DOCTEST_CHECK(fabs(cosh(x) * cosh(x) - sinh(x) * sinh(x) - 1) < 1e-99);
  DOCTEST_CHECK(fabs(sinh(x) / cosh(x) - tanh(x)) < 1e-99);

  DOCTEST_CHECK(fabs(asin(sin(x)) - x) < 1e-99);
  DOCTEST_CHECK(fabs(acos(cos(x)) - x) < 1e-99);

  {
    auto sc = sin_cos(x);
    DOCTEST_CHECK(sc.sin == sin(x));
    DOCTEST_CHECK(sc.cos == cos(x));
  }
  {
    auto sc = sinh_cosh(x);
    DOCTEST_CHECK(sc.sinh == sinh(x));
    DOCTEST_CHECK(sc.cosh == cosh(x));
  }

  DOCTEST_CHECK(pow(x, scalar_t{2}) == x * x);
  DOCTEST_CHECK(pow(x, scalar_t{3}) == x * x * x);
  DOCTEST_CHECK(pow(x, scalar_t{0.5}) == sqrt(x));

  DOCTEST_CHECK(nextabove(x) > x);
  DOCTEST_CHECK(nextbelow(x) < x);
  DOCTEST_CHECK(nextbelow(nextabove(x)) == x);

  DOCTEST_CHECK(isnan(nexttoward(x, std::numeric_limits<scalar_t>::quiet_NaN())));
  DOCTEST_CHECK(nexttoward(x, scalar_t{0}) == nextbelow(x));
  DOCTEST_CHECK(nexttoward(x, scalar_t{2}) == nextabove(x));
  DOCTEST_CHECK(nextabove(scalar_t{1}) == std::numeric_limits<scalar_t>::epsilon() + 1);

  DOCTEST_CHECK(atan2(y, x) == atan(y / x));
  DOCTEST_CHECK(atan2(-y, x) == -atan(y / x));
  DOCTEST_CHECK(fabs(sin(atan2(-y, -x) - atan(y / x))) < 1e-99);
  DOCTEST_CHECK(fabs(cos(atan2(-y, -x) - atan(y / x)) + 1) < 1e-99);
  DOCTEST_CHECK(fabs(sin(atan2(y, -x) + atan2(y, x))) < 1e-99);
  DOCTEST_CHECK(fabs(cos(atan2(y, -x) + atan2(y, x)) + 1) < 1e-99);

  DOCTEST_CHECK(ceil(x) == 2);
  DOCTEST_CHECK(floor(x) == 1);
  DOCTEST_CHECK(round(x) == 1);
  DOCTEST_CHECK(trunc(x) == 1);
  DOCTEST_CHECK(round(scalar_t{1.5}) == 2);
  DOCTEST_CHECK(trunc(scalar_t{1.5}) == 1);
  DOCTEST_CHECK(round(scalar_t{-1.5}) == -2);
  DOCTEST_CHECK(trunc(scalar_t{-1.5}) == -1);

  std::fesetround(FE_UPWARD);
  DOCTEST_CHECK(rint(x) == 2);
  std::fesetround(FE_DOWNWARD);
  DOCTEST_CHECK(rint(x) == 1);
  std::fesetround(FE_TOWARDZERO);
  DOCTEST_CHECK(rint(x) == 1);
  DOCTEST_CHECK(rint(scalar_t{1.5}) == 1);
  std::fesetround(FE_TONEAREST);
  DOCTEST_CHECK(rint(x) == 1);
  DOCTEST_CHECK(rint(scalar_t{1.5}) == 2);
}
