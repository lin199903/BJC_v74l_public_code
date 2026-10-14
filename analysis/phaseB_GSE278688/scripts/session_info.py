#!/usr/bin/env python
"""session_info.py -- exact interpreter and package versions for the frozen run."""
import sys, platform
print("=== Python ===")
print("version   :", sys.version.replace("\n", " "))
print("executable:", sys.executable)
print("platform  :", platform.platform())
print()
print("=== packages used ===")
for name in ("numpy", "pandas", "scipy", "anndata", "scanpy", "h5py", "statsmodels"):
    try:
        m = __import__(name)
        print(f"  {name:14s} {getattr(m, '__version__', 'unknown')}")
    except Exception as e:
        print(f"  {name:14s} NOT INSTALLED ({type(e).__name__})")
