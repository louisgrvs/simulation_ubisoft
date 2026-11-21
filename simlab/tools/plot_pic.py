"""Trace les mesures de apic_demo (jour 11) : PIC contre FLIP (alpha = 0,95 et 1) contre APIC.

apic.png, disque en rotation :
    ligne du haut   transferts seuls (disk_transfers.csv) : 1 - L/L0 (échelle log) et E/E0
                    en fonction du nombre d'allers-retours ; puis simulation complète :
                    Lg/Lg0, moment cinétique de la grille après P2G
    ligne du bas    simulation complète (disk_free.csv) : E/E0, R_g/R_g0 (forme), bruit relatif.
                    La bande grise commence au quart de tour : au-delà, le disque se défait
                    (instabilité de Rayleigh-Taylor) et la conservation n'est plus mesurable.

stability.png, balayage en pas de temps (stability.csv) :
    E_final / E0 en fonction du CFL maximal atteint, pour la rupture de barrage et le
    disque libre. Trait plein : avec sous-pas d'advection ; tirets : un seul déplacement
    par pas. Croix : calcul explosé (valeur non finie, vitesse > 10 fois l'échelle ou
    E > 1,05 E0).

Usage : python tools/plot_apic.py <sorties_apic> [omega]
        (par défaut : sorties_apic, omega = 2 rad/s)
"""
import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "sorties_apic")
omega = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
t_quarter = 0.25 * 2 * np.pi / omega   # le disque libre tient à peu près jusque-là

# Le dossier est créé par apic_demo, dans le dossier d'où on l'a lancé (pas forcément simlab/).
if not d.is_dir() or not any((d / f).exists() for f in ("disk_transfers.csv", "disk_free.csv", "stability.csv")):
    sys.exit(f"Aucun CSV de apic_demo dans « {d.resolve()} ».\n"
             "Lance apic_demo depuis ce dossier, ou donne au script le chemin du dossier qu'il a créé\n"
             "(par exemple out\\build\\x64-Release\\sorties_apic si la démo a tourné depuis Visual Studio).")


def read(path, key="schema"):
    """CSV -> {clé: {colonne: tableau numpy}}, dans l'ordre du fichier."""
    out = OrderedDict()
    if not path.exists():
        return out
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            s = out.setdefault(row[key], {})
            for k, v in row.items():
                if k != key:
                    try:
                        s.setdefault(k, []).append(float(v))
                    except ValueError:
                        s.setdefault(k, []).append(v)
    return OrderedDict((s, {k: np.array(v) for k, v in cols.items()}) for s, cols in out.items())


def missing(ax, name):
    ax.text(0.5, 0.5, f"{name} absent ou vide", ha="center", va="center", transform=ax.transAxes,
            color="#8a8a85")


# Palette catégorielle validée pour le fond sombre (ordre fixe : PIC, FLIP 0,95, FLIP 1, APIC).
SURFACE = "#1a1a19"
INK, INK2 = "#f0efec", "#c3c2b7"
colors = OrderedDict([("PIC", "#3987e5"), ("FLIP 0.95", "#d95926"), ("FLIP 1", "#199e70"),
                      ("APIC", "#c98500")])
styles = {"PIC": "-", "FLIP 0.95": "--", "FLIP 1": "-.", "APIC": "-"}   # identité pas seulement par la couleur

plt.style.use("dark_background")
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": "#4a4a46", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlecolor": INK, "grid.color": "#3a3a37", "legend.frameon": False,
})


def finish(ax, xlabel):
    ax.set_xlabel(xlabel)
    ax.grid(alpha=0.6, lw=0.6)
    if ax.lines:
        ax.legend(fontsize=9)


def shade_unstable(ax):
    """Bande grise après le quart de tour : le disque n'est plus un disque."""
    ax.axvspan(t_quarter, ax.get_xlim()[1], color="#8a8a85", alpha=0.12, lw=0)
    ax.axvline(t_quarter, color="#8a8a85", lw=0.8, ls=":")


# =============================================================================
#  Figure 1 : le disque en rotation
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
(a1, a2, a3), (a4, a5, a6) = axes

# --- Transferts seuls ----------------------------------------------------------------
tr = read(d / "disk_transfers.csv")
if tr:
    for s, r in tr.items():
        c, ls = colors.get(s, INK), styles.get(s, "-")
        k = r["k"][1:]
        # FLIP 1 et APIC sont au bruit machine : on les plaque à 1e-16 pour les voir sur l'échelle log.
        a1.semilogy(k, np.maximum(1 - r["L/L0"][1:], 1e-16), color=c, ls=ls, lw=2, label=s)
        a2.plot(r["k"], r["E/E0"], color=c, ls=ls, lw=2, label=s)
    a1.axhline(1e-16, color="#8a8a85", lw=0.8, ls=":")
    a1.text(1, 2e-16, "précision machine (FLIP 1, APIC)", color=INK2, fontsize=8, va="bottom")
    a1.set_title("Transferts seuls : perte de moment cinétique 1 − L/L0")
    a2.set_title("Transferts seuls : énergie E/E0")
    finish(a1, "allers-retours P2G + G2P")
    finish(a2, "allers-retours P2G + G2P")
else:
    missing(a1, "disk_transfers.csv")
    missing(a2, "disk_transfers.csv")

# --- Simulation complète -------------------------------------------------------------
fr = read(d / "disk_free.csv")
if fr:
    for s, r in fr.items():
        c, ls = colors.get(s, INK), styles.get(s, "-")
        a3.plot(r["t"], r["Lg/Lg0"], color=c, ls=ls, lw=2, label=s)
        a4.plot(r["t"], r["E/E0"], color=c, ls=ls, lw=2, label=s)
        a5.plot(r["t"], r["Rg/Rg0"], color=c, ls=ls, lw=2, label=s)
        a6.semilogy(r["t"][1:], np.maximum(r["bruit_relatif"][1:], 1e-12), color=c, ls=ls, lw=1.5, label=s)
    a3.set_title("Disque libre : moment cinétique de la grille Lg/Lg0")
    a4.set_title("Disque libre : énergie E/E0")
    a5.set_title("Disque libre : forme R_g / R_g(0) (1 = disque intact)")
    a6.set_title("Disque libre : bruit <σ²> / <|u|²>")
    for ax in (a3, a4, a5, a6):
        finish(ax, "t (s)")
        shade_unstable(ax)
    a3.text(t_quarter, 0.02, "  ¼ de tour", transform=a3.get_xaxis_transform(), color=INK2, fontsize=8)
else:
    for ax in (a3, a4, a5, a6):
        missing(ax, "disk_free.csv")

fig.suptitle("Jour 11 : disque en rotation, PIC, FLIP et APIC (même grille, même pas, même graine)", fontsize=14)
fig.tight_layout()
fig.savefig(d / "apic.png", dpi=130)
print(f"figure : {d / 'apic.png'}")

# =============================================================================
#  Figure 2 : stabilité selon le pas de temps
# =============================================================================
st = read(d / "stability.csv", key="scene")
fig, axes = plt.subplots(1, 2, figsize=(14, 6.2))
titles = {"barrage": "Rupture de barrage (1 s)", "disque_libre": "Disque libre (¼ de tour)"}
Y_BLOWN = 1.12   # les calculs explosés sont rangés dans une bande au-dessus de E0 (leur E réel n'a pas de sens)
for ax, scene in zip(axes, ["barrage", "disque_libre"]):
    if scene not in st:
        missing(ax, "stability.csv")
        continue
    r = st[scene]
    for i, (s, c) in enumerate(colors.items()):
        for sub, ls in ((1.0, "-"), (0.0, "--")):
            m = (r["schema"] == s) & (r["sous_pas"] == sub)
            if not m.any():
                continue
            order = np.argsort(r["cfl_max"][m])
            x, y, blown = r["cfl_max"][m][order], r["E_final/E0"][m][order], r["explose"][m][order] > 0
            label = f"{s}, {'avec' if sub else 'sans'} sous-pas"
            ax.plot(x[~blown], y[~blown], color=c, ls=ls, lw=2, marker="o", ms=5, label=label)
            if blown.any():   # décalés verticalement par schéma pour ne pas se recouvrir
                ax.plot(x[blown], np.full(blown.sum(), Y_BLOWN + 0.035 * i), ls="none", marker="X", ms=11,
                        color=c, markeredgecolor=SURFACE, markeredgewidth=1.5, clip_on=False)
    ax.plot([], [], ls="none", marker="X", ms=10, color=INK2, label="explosé")
    ax.axhline(1, color="#8a8a85", lw=0.8, ls=":")
    ax.axhspan(1.07, 1.3, color="#8a8a85", alpha=0.12, lw=0)
    ax.text(0.01, 1.085, "explosés", transform=ax.get_yaxis_transform(), color=INK2, fontsize=8)
    ax.set_ylim(0, 1.3)
    ax.set_xscale("log")
    ax.set_ylabel("E_final / E0")
    ax.set_title(titles[scene] + " : énergie finale selon le CFL atteint")
    ax.set_xlabel("CFL max = max|u_p| dt / h")
    ax.grid(alpha=0.6, lw=0.6)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=9, handlelength=3.5)
fig.suptitle("Jour 11 : stabilité de chaque schéma selon le pas de temps "
             "(trait plein : avec sous-pas ; tirets : un seul déplacement par pas)", fontsize=13)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(d / "stability.png", dpi=130)
print(f"figure : {d / 'stability.png'}")