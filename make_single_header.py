files = [
    "./include/mpfr/detail/hedley.h",
    "./include/mpfr/detail/mpfr.hpp",
    "./include/mpfr/mpfr.hpp",
]

with open("./mpfr.hpp", "w") as single_header:
    for f in files:
        with open(f) as header:
            single_header.write(header.read())
