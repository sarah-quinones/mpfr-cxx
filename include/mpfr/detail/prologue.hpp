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
    HEDLEY_HAS_BUILTIN(__builtin_signbit)

#define MPFR_CXX_HAS_MATH_BUILTINS 1

#else

#define MPFR_CXX_HAS_MATH_BUILTINS 0

#endif
