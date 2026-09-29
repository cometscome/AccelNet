"""Check that every symbol declared in a public C header is dynamically visible.

This catches private BIND(C) symbols hidden by GCC PR126872, even for API
functions not exercised by the numerical tests. No function is invoked.
"""
import ctypes
from pathlib import Path
import re
import sys

library_path, header_path = sys.argv[1:]
header = re.sub(r"/\*.*?\*/|//[^\n]*", "", Path(header_path).read_text(), flags=re.S)
functions = set(re.findall(r"\b(accelnet_\w+)\s*\(", header))
variables = set(re.findall(
    r"^extern\s+(?:int|double|ACCELNET_BOOL)\s+(\w+)\s*;", header, flags=re.M))
if not functions:
    raise SystemExit(f"No C function declarations found in {header_path}")

library = ctypes.CDLL(library_path)
missing = []
for name in sorted(functions):
    try:
        getattr(library, name)
    except AttributeError:
        missing.append(name)
for name in sorted(variables):
    try:
        ctypes.c_byte.in_dll(library, name)
    except ValueError:
        missing.append(name)
if missing:
    raise SystemExit(f"Missing exports in {library_path}: {', '.join(missing)}")
print(f"Exported {len(functions)} functions and {len(variables)} variables from {library_path}")
