MPFR layer
==========

This function provides a (potentially unsafe) way to read and assign the
:code:`mp_float_t<_>` types through the MPFR C interface.

.. doxygenfunction:: mpfr::handle_as_mpfr_t

Example usage:

.. literalinclude:: ../test/mpfr_layer.cpp
  :language: cpp
  :start-after: [mpfr-layer-test begin]
  :end-before: [mpfr-layer-test end]
  :lines: 3-6,9-,7
