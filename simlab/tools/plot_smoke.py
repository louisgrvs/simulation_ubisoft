"""Trace les mesures de smoke_demo (jour 5).

Gauche  : E(t)/E(0) pour plusieurs viscosités ν, avec en pointillés le facteur
          visqueux exact exp(-2 ν k² t) multiplié par la courbe ν = 0
          (les deux dissipations, physique et numérique, se composent).
Milieu  : viscosité numérique effective en fonction du CFL (ν = 0).
Droite  : viscosité numérique effective en fonction de n (ν = 0, dt fixé).

Usage : python tools/plot_smoke.py <dossier_des_sorties>
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "out_smoke")


def read(name):
    return np.genfromtxt(d / name, delimiter=",", names=True, dtype=None, encoding="utf-8")


nu_runs = read("energy_nu.csv")
summary = read("summary.csv")

plt.style.use("dark_background")
colors = ["#7fd1ff", "#ffb86b", "#ff6e9c", "#b6f07a"]
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.4))

# --- 1) Viscosité --------------------------------------------------------------
series = list(dict.fromkeys(nu_runs["serie"]))
base = nu_runs[nu_runs["serie"] == series[0]]          # ν = 0 : dissipation numérique seule
for s, c in zip(series, colors):
    r = nu_runs[nu_runs["serie"] == s]
    a1.plot(r["t"], r["E_ratio"], color=c, lw=2, label=s.replace("nu", "ν"))
    if s != series[0]:
        a1.plot(r["t"], base["E_ratio"] * r["exact_viscous"], color=c, ls="--", lw=1)
a1.set_xlabel("t (s)")
a1.set_ylabel("E(t) / E(0)")
a1.set_title("Tourbillon stationnaire, n = 64, dt = 0,01 s\n(pointillés : courbe ν = 0 × exp(-2νk²t))")
a1.legend()

# --- 2) Pas de temps -------------------------------------------------------------
p = summary[summary["serie"] == "pas"]
a2.loglog(p["cfl"], p["nu_num"], "o-", color=colors[1])
for cfl, nu, dt in zip(p["cfl"], p["nu_num"], p["dt"]):
    a2.annotate(f"dt={dt:g}", (cfl, nu), textcoords="offset points", xytext=(4, 6), fontsize=8, color="0.75")
a2.set_xlabel("CFL = max|u| dt / h")
a2.set_ylabel("ν numérique effective (m²/s)")
a2.set_title("Dissipation numérique selon le pas (ν = 0)")

# --- 3) Résolution ---------------------------------------------------------------
r = summary[summary["serie"] == "resolution"]
a3.loglog(r["n"], r["nu_num"], "o-", color=colors[0], label="mesuré")
slope = np.polyfit(np.log(r["n"]), np.log(r["nu_num"]), 1)[0]
a3.loglog(r["n"], r["nu_num"][0] * r["n"][0] / r["n"], "--", color="0.6", label="pente -1")
a3.set_xlabel("n (grille n × n)")
a3.set_ylabel("ν numérique effective (m²/s)")
a3.set_title(f"Dissipation numérique selon h : pente {slope:.2f}")
a3.legend()

fig.tight_layout()
out = d / "smoke_energy.png"
fig.savefig(out, dpi=130)
print(f"pente en n : {slope:.2f} -> {out}")
