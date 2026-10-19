"""Spread of the certified critical constant c*(theta) across the declared box (Table 2).
Rigorous brackets from the rate engine (formal01_f8a_battery.cstar_bracket).
Writes deposits/cstar_spread_claims.txt."""
import pathlib
from fractions import Fraction as F
import formal01_f8a_battery as B
pts = [("corner (argmax)", (F(1,25), F(-3,2), F(3,2))),
       ("no loadings, equal var", (F(1), F(0), F(0))),
       ("same sign 1,1", (F(1), F(1), F(1))),
       ("same sign 1.5,1.5, lam 25", (F(25), F(3,2), F(3,2))),
       ("opposite 1,-1", (F(1), F(1), F(-1))),
       ("benchmark unloaded 0,1", (F(1), F(0), F(1)))]
rows = []
for name, th in pts:
    lo, hi = B.cstar_bracket(th)
    rows.append("CLAIM SPREAD-%s: cstar bracket at theta=%s = [%.6f, %.6f]" % (name, tuple(float(x) for x in th), float(lo), float(hi)))
out = pathlib.Path(__file__).resolve().parent / "deposits" / "cstar_spread_claims.txt"
out.write_text("\n".join(rows) + "\n", encoding="utf-8")
print("\n".join(rows))
