MPFR layer
==========

This function provides a (potentially unsafe) way to read and assign the
:code:`mp_float_t<_>` types through the MPFR C interface.

.. doxygenenum:: mpfr::parameter_type
.. doxygenfunction:: mpfr::handle_as_mpfr_t

Example usage:

.. code-block:: cpp

    using scalar512_t = mp_float_t<digits2{512}>;
    using scalar1024_t = mp_float_t<digits2{1024}>;
    scalar512_t x; // uninitialized
    scalar1024_t y = scalar1024_t{2.0 + 1.0/1024};
    handle_as_mpfr_t<out, in>(
      [](mpfr_ptr a, mpfr_srcptr b) {
        mpfr_set_str(a, "0.1111", 2, MPFR_RNDN); // set to 0b0.1111 == 0.9375
        mpfr_add(a, a, b, MPFR_RNDN);
      },
      x, y
    );
    assert(x == scalar512_t{2 + 0.5 + 0.25 + 0.125 + 0.0625 + 1.0/1024});
