"""Strict-JSON writer for 03-certified-supremum — re-export of `02-bounded-rule/io_util.py`, single implementation.

`02-bounded-rule/io_util.py` is the canonical writer (its `dump_json(..., strict=True)` refuses
non-finite floats instead of emitting bare `NaN`/`Infinity`).  03-certified-supremum imports it through
this shim so that imports resolve inside this directory and the two trees cannot drift.
"""

from __future__ import annotations

import importlib.util
import os

_W26 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "02-bounded-rule")
_SRC = os.path.join(_W26, "io_util.py")

_spec = importlib.util.spec_from_file_location("w26_io_util", _SRC)
if _spec is None or _spec.loader is None:  # pragma: no cover
    raise ImportError(f"cannot load {_SRC}")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

load_json = _mod.load_json
dump_json = _mod.dump_json
f64 = _mod.f64
i64 = _mod.i64
rnd = _mod.rnd
write_text = _mod.write_text

__all__ = ["load_json", "dump_json", "f64", "i64", "rnd", "write_text"]
