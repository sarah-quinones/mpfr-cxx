lightweight wrapper around [gnu mpfr](https://www.mpfr.org/). allowing usage of stack allocated multiprecision floats.  
benefits over boost/multiprecision/mpfr:
  - scalar types are [trivially copyable](https://en.cppreference.com/w/cpp/named_req/TriviallyCopyable).
  - a possible value representation of positive zero is all zero bits. i.e., for an object `x` of type `mp_float_t<_>`, `std::memset(&x, 0, sizeof(x))`, will set `x` to `+0.0`.
  - optimized for certain common operations (for instance, multiplication and division by powers of two).
  - better `ostream` formatting support.
  - smaller compilation overhead.
 
disadvantages:
  - no support for operations involving floats of different precisions.

# usage
run `make_single_header.py` to generate a single header include file.
