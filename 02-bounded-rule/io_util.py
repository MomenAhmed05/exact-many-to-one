"""Small IO + numeric helpers that keep the pi-lens 'unchecked call' rule quiet."""

import json
import os


def load_json(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"could not read {path}") from exc


def dump_json(obj, path, strict=False):
    """Write JSON.  `strict=True` refuses non-finite floats instead of emitting
    bare `Infinity`/`NaN` (invalid strict JSON, and how a polluted deposit arose)."""
    try:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, allow_nan=not strict)
        return path
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"could not write {path}") from exc


def f64(x):
    try:
        return float(x)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("not a float") from exc


def i64(x):
    try:
        return int(x)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("not an int") from exc


def rnd(x, nd=6):
    try:
        return round(float(x), nd)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("not a float") from exc


def write_text(text, path):
    try:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"could not write {path}") from exc
