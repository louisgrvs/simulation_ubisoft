"""Trace les résultats de projection_demo.
Gauche : résidu ||r||_∞ du gradient conjugué au fil des itérations.
Droite : nombre d'itérations en fonction de n (pente attendue ≈ 1, car κ(A) ~ n²
et CG demande O(√κ) itérations).
Usage : python tools/plot_projection.py <dossier_des_sorties>"""
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
res = np.genfromtxt(d / "projection_residuals.csv", delimiter=",", names=True)
summ = np.genfromtxt(d / "projection_summary.csv", delimiter=",", names=True)

plt.style.use("dark_background")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
for n in np.unique(res["n"]):
    m = res["n"] == n
    a1.semilogy(res["iteration"][m], res["residual"][m], label=f"n = {int(n)}")
a1.axhline(1e-8, ls="--", c="0.6", lw=1, label="tolérance")
a1.set_xlabel("itération"); a1.set_ylabel(r"$\|r\|_\infty$ = max |div| (s$^{-1}$)")
a1.set_title("Convergence du gradient conjugué"); a1.legend()

n, it = summ["n"], summ["iterations"]
a2.loglog(n, it, "o-", label="mesuré")
slope = np.polyfit(np.log(n), np.log(it), 1)[0]
a2.loglog(n, it[0] * n / n[0], "--", c="0.6", label="pente 1")
a2.set_xlabel("n (grille n × n)"); a2.set_ylabel("itérations")
a2.set_title(f"Itérations de CG : pente mesurée {slope:.2f}"); a2.legend()

fig.tight_layout()
fig.savefig(d / "projection.png", dpi=130)
print(f"pente mesurée : {slope:.2f} -> {d / 'projection.png'}")
