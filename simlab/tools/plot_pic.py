"""Trace les mesures de pic_transfer_demo (jour 7).

En haut à gauche  : L/L0 (traits pleins) et E/E0 (pointillés) au fil des
                    allers-retours P2G -> G2P, pour n = 32, 64, 128.
En haut à droite  : L/L0 selon le désordre de l'ensemencement et le nombre
                    de particules par cellule (n = 64).
En bas à gauche   : erreur de vitesse après UN aller-retour, selon r/R,
                    réseau aléatoire contre réseau régulier.
En bas à droite   : tourbillon cellulaire en boîte pleine, solveur PIC complet ;
                    si l'on donne le dossier de smoke_demo (jour 5), les courbes
                    de Stable Fluids sont superposées pour comparaison.

Usage : python tools/plot_pic.py <sorties_pic> [<sorties_smoke_jour5>]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "out_pic")
smoke = Path(sys.argv[2]) if len(sys.argv) > 2 else None


def read(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def series(data):
    return list(dict.fromkeys(data["serie"]))


plt.style.use("dark_background")
colors = ["#7fd1ff", "#ffb86b", "#ff6e9c", "#b6f07a"]
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
(a1, a2), (a3, a4) = axes

# --- 1) Résolution -----------------------------------------------------------------
res = read(d / "transfers_n.csv")
first_loss, ns = [], []
for s, c in zip(series(res), colors):
    r = res[res["serie"] == s]
    a1.plot(r["k"], r["L_ratio"], color=c, lw=2, label=f"{s} : L/L0")
    a1.plot(r["k"], r["E_ratio"], color=c, lw=1, ls="--")
    ns.append(int(s.split("=")[1]))
    first_loss.append(1 - r["L_ratio"][1])
slope = np.polyfit(np.log(ns), np.log(first_loss), 1)[0]
a1.set_xlabel("allers-retours P2G -> G2P")
a1.set_ylabel("rapport à la valeur initiale")
a1.set_title(f"Disque en rotation, transferts PIC seuls\n"
             f"(pointillés : E/E0 ; perte de L par aller-retour en n^{slope:.2f})")
a1.legend()

# --- 2) Ensemencement ----------------------------------------------------------------
jit = read(d / "transfers_jitter.csv")
ppc = read(d / "transfers_ppc.csv")
for s, c in zip(series(jit), colors):
    r = jit[jit["serie"] == s]
    a2.plot(r["k"], r["L_ratio"], color=c, lw=2, label=s)
for s, c in zip(series(ppc), ["#d9d9d9", "#9be3c4", "#c9a7ff"]):
    r = ppc[ppc["serie"] == s]
    a2.plot(r["k"], r["L_ratio"], color=c, lw=1, ls=":", label=s)
a2.set_xlabel("allers-retours P2G -> G2P")
a2.set_ylabel("L / L0")
a2.set_title("n = 64 : désordre (traits pleins) et particules par cellule (pointillés)")
a2.legend(fontsize=8)

# --- 3) Localisation de l'erreur -----------------------------------------------------
rad = read(d / "radial.csv")
for j, c, lab in [(1.0, colors[0], "aléatoire (jitter 1)"), (0.0, colors[1], "réseau régulier (jitter 0)")]:
    r = rad[rad["jitter"] == j]
    # Les erreurs nulles (réseau régulier, intérieur) sont posées à 1e-16 pour rester visibles.
    a3.scatter(r["r_over_R"], np.maximum(r["err_over_omegaR"], 1e-16), s=3, color=c, alpha=0.5, label=lab)
a3.set_yscale("log")
a3.set_ylim(3e-17, 1e-1)
a3.set_xlabel("r / R")
a3.set_ylabel("|u' - u| / (ω R) après un aller-retour")
a3.set_title("Où PIC se trompe : bord du disque, et bruit si le réseau est aléatoire")
a3.legend(markerscale=4)

# --- 4) Solveur complet ----------------------------------------------------------------
vor = read(d / "vortex_pic.csv")
for s, c in zip(series(vor), colors):
    r = vor[vor["serie"] == s]
    a4.plot(r["t"], r["E_ratio"], color=c, lw=2, label=f"PIC, {s}")
if smoke is not None and (smoke / "energy_n.csv").exists():
    st = read(smoke / "energy_n.csv")
    for s, c in zip(series(vor), colors):
        r = st[st["serie"] == s]
        if len(r):
            a4.plot(r["t"], r["E_ratio"], color=c, lw=1, ls="--", label=f"Stable Fluids, {s}")
a4.set_xlabel("t (s)")
a4.set_ylabel("E(t) / E(0)")
a4.set_title("Tourbillon cellulaire en boîte pleine, dt = 0,01 s, ν = 0")
a4.legend()

fig.tight_layout()
out = d / "pic_transfers.png"
fig.savefig(out, dpi=120)
print(f"perte de L par aller-retour en n^{slope:.2f} -> {out}")
