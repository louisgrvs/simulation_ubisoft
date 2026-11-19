"""Trace les mesures de flip_demo (jours 9 et 10) : PIC contre FLIP (alpha = 0,95, 0,99, 1) contre APIC.

Ligne du haut, rupture de barrage (dam.csv) :
    (E_c + E_p) / E0   énergie mécanique totale : ce qui manque a été dissipé
    E_c / E0           énergie cinétique seule (la chute la convertit, l'impact la détruit)
    bruit relatif      <sigma²> / <|u|²> : variance des vitesses dans chaque cellule
Ligne du bas :
    tourbillon en boîte pleine (vortex.csv) : E/E0, puis bruit relatif
    front de la rupture de barrage, en variables de Martin & Moyce :
    Z = x_front / a  contre  T = t sqrt(2 g / a),  a = largeur initiale de la colonne

Usage : python tools/plot_flip.py <sorties_flip> [largeur_colonne_a] [g]
        (par défaut : a = 0,25 m, g = 9,81 m/s²)
"""
import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "sorties_flip")
a_col = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
g = float(sys.argv[3]) if len(sys.argv) > 3 else 9.81


def read(path):
    """CSV -> {schéma: {colonne: tableau numpy}}, dans l'ordre du fichier."""
    out = OrderedDict()
    if not path.exists():
        return out
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            s = out.setdefault(row["schema"], {})
            for k, v in row.items():
                if k != "schema":
                    s.setdefault(k, []).append(float(v))
    return OrderedDict((s, {k: np.array(v) for k, v in cols.items()}) for s, cols in out.items())


plt.style.use("dark_background")
colors = {"PIC": "#7fd1ff", "FLIP 0.95": "#b6f07a", "FLIP 0.99": "#ffb86b", "FLIP 1": "#ff6e9c",
          "APIC": "#e0e0e0"}
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
(a1, a2, a3), (a4, a5, a6) = axes

# --- Rupture de barrage ----------------------------------------------------------------
dam = read(d / "dam.csv")
if dam:
    for s, r in dam.items():
        c = colors.get(s, "white")
        E0 = r["Etot"][0]
        a1.plot(r["t"], r["Etot/E0"], color=c, lw=2, label=s)
        a2.plot(r["t"], r["Ec"] / E0, color=c, lw=2, label=s)
        a3.semilogy(r["t"][1:], np.maximum(r["bruit_relatif"][1:], 1e-12), color=c, lw=1.5, label=s)
        T = r["t"] * np.sqrt(2 * g / a_col)
        a6.plot(T, r["front"] / a_col, color=c, lw=2, label=s)
    a1.set_title("Rupture de barrage : énergie mécanique (E_c + E_p) / E0")
    a2.set_title("Rupture de barrage : énergie cinétique E_c / E0")
    a3.set_title("Rupture de barrage : bruit <σ²> / <|u|²>")
    a6.set_title("Front de l'eau (variables de Martin & Moyce)")
    a6.set_xlabel("T = t √(2g/a)")
    a6.set_ylabel("Z = x_front / a")
    a6.axhline(1 / a_col, color="gray", lw=0.8, ls=":")   # paroi de droite (Lx = 1)
    for ax in (a1, a2, a3):
        ax.set_xlabel("t (s)")
else:
    for ax in (a1, a2, a3, a6):
        ax.text(0.5, 0.5, "dam.csv absent ou vide\n(surface libre du jour 8 nécessaire)",
                ha="center", va="center", transform=ax.transAxes, color="gray")

# --- Tourbillon en boîte pleine ----------------------------------------------------------
vortex = read(d / "vortex.csv")
for s, r in vortex.items():
    c = colors.get(s, "white")
    a4.plot(r["t"], r["E/E0"], color=c, lw=2, label=s)
    a5.semilogy(r["t"][1:], np.maximum(r["bruit_relatif"][1:], 1e-12), color=c, lw=1.5, label=s)
a4.set_title("Tourbillon en boîte pleine : E / E0 (idéal : 1)")
a5.set_title("Tourbillon en boîte pleine : bruit <σ²> / <|u|²>")
for ax in (a4, a5):
    ax.set_xlabel("t (s)")

for ax in axes.flat:
    ax.grid(alpha=0.25)
    if ax.lines:
        ax.legend(fontsize=9)

fig.suptitle("Jours 9 et 10 : PIC, FLIP et APIC (même grille, même pas de temps, même graine)", fontsize=14)
fig.tight_layout()
fig.savefig(d / "flip.png", dpi=130)
print(f"figure : {d / 'flip.png'}")
