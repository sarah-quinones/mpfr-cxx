#ifdef CXX_MPFR_PROLOGUE
#error "epilogue missing"
#endif
#define CXX_MPFR_PROLOGUE

#ifndef CXX_MPFR_DEBUG
#define CXX_MPFR_DEBUG 0
#endif

#if CXX_MPFR_DEBUG == 1
#define CXX_MPFR_STRINGIZE2(...) #__VA_ARGS__
#define CXX_MPFR_STRINGIZE(...) CXX_MPFR_STRINGIZE2(__VA_ARGS__)
#define CXX_MPFR_ASSERT(...)                                                                       \
  (static_cast<bool>(__VA_ARGS__)                                                                  \
       ? (void)0                                                                                   \
       : ::mpfr::_::crash_with_message("assertion_failed at " __FILE__                             \
                                       ":" CXX_MPFR_STRINGIZE(__LINE__) "\n" #__VA_ARGS__))
#else
#define CXX_MPFR_ASSERT(...) ((void)0)
#endif
