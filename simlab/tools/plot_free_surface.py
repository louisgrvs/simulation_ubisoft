"""Trace les mesures de free_surface_demo (jour 8 : surface libre).

À gauche  : cuve au repos, max |u_p| au fil du temps (vitesses parasites),
            comparé à sqrt(g H), la vitesse caractéristique des ondes.
Au milieu : profil de pression au milieu de la cuve, contre la référence
            discrète ρ g (y_s − y) — y_s = centre de la première cellule d'air.
À droite  : cuve agitée, aire fluide / aire initiale, pour 32² et 64²,
            avec (traits pleins) et sans (pointillés) comblement des bulles.

Usage : python tools/plot_free_surface.py [<dossier_csv>] [<image.png>]
        (par défaut : out_free_surface, et l'image free_surface_j8.png dans ce dossier)
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "out_free_surface")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else d / "free_surface_j8.png"

G = 9.81


def read(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


plt.style.use("dark_background")
colors = {32: "#7fd1ff", 64: "#ffb86b"}
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 5))

# --- 1) Vitesses parasites ---------------------------------------------------------------
rest = read(d / "rest.csv")
H = rest["surface_height"][0]
a1.semilogy(rest["t"], np.maximum(rest["max_speed"], 1e-18), color="#b6f07a", lw=2, label="max |u_p|")
a1.semilogy(rest["t"], np.maximum(rest["rms_speed"], 1e-18), color="#b6f07a", lw=1, ls="--", label="rms |u_p|")
a1.axhline(np.sqrt(G * H), color="#ff6e9c", lw=1, ls=":", label=f"√(gH) = {np.sqrt(G * H):.2f} m/s")
a1.set_ylim(1e-14, 10)
a1.set_xlabel("t (s)")
a1.set_ylabel("vitesse (m/s)")
a1.set_title("Cuve au repos, 64 × 64 : vitesses parasites\n"
             f"(aire fluide constante : {rest['fluid_cells'].min()} à {rest['fluid_cells'].max()} cellules)")
a1.legend()

# --- 2) Pression hydrostatique -------------------------------------------------------------
prof = read(d / "pressure.csv")
a2.plot(prof["p_ref"], prof["y"], color="#ff6e9c", lw=6, alpha=0.35, label="ρ g (y_s − y)")
a2.plot(prof["p"], prof["y"], "o", color="#7fd1ff", ms=4, label="p calculée")
err = np.max(np.abs(prof["p"] - prof["p_ref"])) / np.max(prof["p_ref"])
a2.set_xlabel("p (Pa)")
a2.set_ylabel("y (m)")
a2.set_title(f"Profil de pression au milieu de la cuve\n(écart relatif max : {err:.1e})")
a2.legend()

# --- 3) Volume de la cuve agitée -------------------------------------------------------------
for n in (32, 64):
    for bubbles, ls in ((8, "-"), (0, "--")):
        path = d / f"stirred_n{n}_bulles{bubbles}.csv"
        if not path.exists():
            continue
        r = read(path)
        lab = f"{n}², bulles comblées" if bubbles else f"{n}², sans comblement"
        a3.plot(r["t"], r["fluid_cells"] / r["fluid_cells"][0], color=colors[n], lw=2 if bubbles else 1.2,
                ls=ls, label=lab)
a3.axhline(1, color="white", lw=0.6, alpha=0.5)
a3.set_xlabel("t (s)")
a3.set_ylabel("aire fluide / aire initiale")
a3.set_title("Cuve agitée : volume de fluide\n(une bulle à p = 0 implose et fait fondre le volume)")
a3.legend(fontsize=9)

fig.tight_layout()
fig.savefig(out, dpi=130)
print(f"Figure écrite : {out}")
