from pathlib import Path
from io import TextIOWrapper


def append_header_recursively(lines, file_path: Path):
    with open(file_path) as file:
        for line in file.readlines():
            line: str
            if not line.startswith('#include "mpfr/'):
                lines.append(line)
            else:
                nested_header_path = Path(".") / "include"

                start = line.find('"') + 1
                while True:
                    end = line.find("/", start)
                    if end != -1:
                        # directory
                        nested_header_path = nested_header_path / line[start:end]
                        start = end + 1
                    else:
                        # remove trailing newline and end-quote
                        nested_header_path = nested_header_path / line[start:-2]
                        break

                append_header_recursively(lines, nested_header_path)


lines = []
main_header = Path("include") / "mpfr" / "mpfr.hpp"
append_header_recursively(lines, main_header)
print("".join(lines))
