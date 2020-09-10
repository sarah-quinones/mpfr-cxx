#ifdef MPFR_CXX_PROLOGUE
#error "epilogue missing"
#endif
#define MPFR_CXX_PROLOGUE

#ifndef MPFR_CXX_DEBUG
#define MPFR_CXX_DEBUG 0
#endif

#if MPFR_CXX_DEBUG == 1
#define MPFR_CXX_STRINGIZE2(...) #__VA_ARGS__
#define MPFR_CXX_STRINGIZE(...) MPFR_CXX_STRINGIZE2(__VA_ARGS__)
#define MPFR_CXX_ASSERT(...)                                                                       \
  (static_cast<bool>(__VA_ARGS__)                                                                  \
       ? (void)0                                                                                   \
       : ::mpfr::_::crash_with_message("assertion_failed at " __FILE__                             \
                                       ":" MPFR_CXX_STRINGIZE(__LINE__) "\n" #__VA_ARGS__))
#else
#define MPFR_CXX_ASSERT(...) ((void)0)
#endif

#define MPFR_CXX_NODISCARD HEDLEY_DIAGNOSTIC_DISABLE_CPP98_COMPAT_WRAP_(nodiscard)

#if HEDLEY_HAS_BUILTIN(__builtin_fabs) and HEDLEY_HAS_BUILTIN(__builtin_frexp) and                 \
    HEDLEY_HAS_BUILTIN(__builtin_fabsf) and HEDLEY_HAS_BUILTIN(__builtin_frexpf) and               \
    HEDLEY_HAS_BUILTIN(__builtin_fabsl) and HEDLEY_HAS_BUILTIN(__builtin_frexpl) and               \
    HEDLEY_HAS_BUILTIN(__builtin_signbit)

#define MPFR_CXX_HAS_MATH_BUILTINS 1
#define MPFR_CXX_FREXP __builtin_frexp
#define MPFR_CXX_FREXPF __builtin_frexpf
#define MPFR_CXX_FREXPL __builtin_frexpl

#define MPFR_CXX_FABS __builtin_fabs
#define MPFR_CXX_FABSF __builtin_fabsf
#define MPFR_CXX_FABSL __builtin_fabsl

#define MPFR_CXX_SIGNBIT __builtin_signbit

#else

#define MPFR_CXX_HAS_MATH_BUILTINS 0

#define MPFR_CXX_FREXP ::std::frexp
#define MPFR_CXX_FREXPF ::std::frexp
#define MPFR_CXX_FREXPL ::std::frexp

#define MPFR_CXX_FABS ::std::fabs
#define MPFR_CXX_FABSF ::std::fabs
#define MPFR_CXX_FABSL ::std::fabs
#define MPFR_CXX_SIGNBIT ::std::signbit

#endif
