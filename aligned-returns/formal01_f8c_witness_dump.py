"""Dump the R3 cover as a saved, hashed witness for the independent checker.

Imports `build_graded_cover` and edits nothing in the producer.

Two details matter, both for the same reason:
  * write with `write_bytes` on UTF-8 bytes, NOT `write_text`.  In text mode on Windows a
    "\\n" becomes "\\r\\n", and the file's sha256 then differs from `F8C-COVER-SHA256`, which
    hashes the bytes with a bare "\\n".  This is the only thing that makes the witness
    verifiable against the identity R3 already deposited.
  * an `if __name__ == "__main__":` guard, because Windows spawn re-imports this module in
    each worker.
"""
import hashlib
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
WIT = DEP / "f8c_cover_witness-R1.txt"
CLAIMS = DEP / "f8c_witness_claims-R1.txt"
R3 = DEP / "formal01_f8c_claims-R3.txt"
NL = b"\n"


def _load_expected_sha():
    """WIT-R3-SHA256, parsed from the R3 deposit rather than typed."""
    for line in R3.read_text(encoding="utf-8").splitlines():
        if line.startswith("CLAIM F8C-COVER-SHA256:"):
            return line.split("=", 1)[1].strip()
    raise SystemExit("F8C-COVER-SHA256 not found in the R3 deposit")


def main():
    t0 = time.time()
    # Import the producer under its REAL module name.  A first attempt loaded it through
    # importlib under the synthetic name "f8c_prod"; on Windows, spawn re-imports the module
    # in each worker by name, so the workers could not resolve it and the run died with
    # OSError [WinError 87] before the witness was written.
    sys.path.insert(0, str(ROOT))
    import formal01_f8c_certified as prod

    cells, unserved = prod.build_graded_cover()

    # A-2: exactly the line E-5 hashed, in list order.
    payload = bytearray()
    for c in cells:
        line = ("|".join(str(iv) for iv in c["box"])
                + "#" + str(c["rho1"]) + "#" + str(c["rho2"]))
        payload += line.encode("utf-8") + NL

    # A-3: bytes, not text.
    WIT.write_bytes(bytes(payload))

    # A-4: every value computed, nothing typed.
    back = WIT.read_bytes()                       # read back from disk
    got = hashlib.sha256(back).hexdigest()
    want = _load_expected_sha()
    nlines = back.count(NL)

    out = [
        f"CLAIM WIT-LINES: lines in the witness = {nlines}",
        f"CLAIM WIT-BYTES: bytes in the witness = {len(back)}",
        f"CLAIM WIT-SHA256: sha256 of the witness bytes read back from disk = {got}",
        f"CLAIM WIT-R3-SHA256: F8C-COVER-SHA256 as deposited in R3 = {want}",
        f"CLAIM WIT-MATCH: witness hash equals the deposited cover identity = {got == want}",
        f"CLAIM WIT-UNSERVED: leaves that could not be certified = {len(unserved)}",
        f"CLAIM WIT-SEC: seconds to build the cover and dump the witness = "
        f"{round(time.time() - t0, 1)}",
        f"CLAIM WIT-GATE-W: Gate W requires WIT-MATCH = True and WIT-LINES = 59276 = "
        f"{got == want and nlines == 59276}",
    ]
    CLAIMS.write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
