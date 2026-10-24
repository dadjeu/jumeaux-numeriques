"""Prédictions du jumeau recalé, résumées en une figure.

1. Fréquence propre en fonction de la masse ajoutée en bout.
2. Réponse en régime forcé (force de 0,15 N en bout) : flèche et contrainte maximale
   en fonction de la fréquence d'excitation, avec les zones d'alerte.

Sortie : resultats/predictions.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import parametres as P
from fe_poutre import modes_propres, evaluer, masse_modale, contrainte_flexion

RACINE = Path(__file__).resolve().parent.parent
FORCE = 0.15


def main():
    rec = json.loads((RACINE / "resultats" / "parametres_recales.json").read_text(encoding="utf-8"))
    E, zeta = rec["E_effectif_Pa"], rec["amortissement"]

    masses = np.linspace(0, 150, 31)
    f_masse = [modes_propres(E=E, masse_bout=P.masse_bloc + m / 1000, n_modes=1)[0][0] for m in masses]

    freqs, modes, le = modes_propres(E=E, n_modes=1)
    f1, q = freqs[0], modes[0]
    k = (2 * np.pi * f1) ** 2 * masse_modale(q, E, P.masse_bloc)
    _, courbure = evaluer(q, le, np.linspace(0, P.L, 400))
    sigma_par_m = np.abs(contrainte_flexion(E, courbure)).max()
    fe = np.linspace(0.2, 12, 1000)
    r = fe / f1
    fleche = FORCE / k / np.sqrt((1 - r**2) ** 2 + (2 * zeta * r) ** 2)
    sigma = sigma_par_m * fleche / 1e6

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    ax[0].plot(masses, f_masse, "o-", ms=3)
    ax[0].set(xlabel="masse ajoutée en bout [g]", ylabel="fréquence propre prédite [Hz]",
              title="Fréquence propre et masse ajoutée")
    ax[0].grid(alpha=0.3)

    ax[1].semilogy(fe, sigma, color="C0")
    ax[1].axhline(P.sigma_e / 1e6, color="C3", ls="--", label=f"limite élastique {P.sigma_e/1e6:.0f} MPa")
    ax[1].axhline(P.sigma_e / P.coef_securite / 1e6, color="C1", ls=":", label="limite / 1,5")
    ax[1].axvspan(f1 * (1 - P.seuil_resonance), f1 * (1 + P.seuil_resonance), color="C3", alpha=0.12,
                  label=f"zone de résonance ({f1:.2f} Hz ± 10 %)")
    danger = fe[sigma > P.sigma_e / 1e6]
    if danger.size:
        ax[1].set_title(f"Plastification prédite entre {danger.min():.2f} et {danger.max():.2f} Hz")
    ax[1].set(xlabel="fréquence d'excitation [Hz]", ylabel="contrainte max à l'encastrement [MPa]")
    ax[1].legend(fontsize=8, loc="upper right")
    ax[1].grid(alpha=0.3, which="both")
    fig.savefig(RACINE / "resultats" / "predictions.png", dpi=130)

    print(f"f propre sans masse ajoutée : {f_masse[0]:.2f} Hz ; avec 100 g : "
          f"{np.interp(100, masses, f_masse):.2f} Hz")
    if danger.size:
        print(f"Plastification prédite pour une excitation entre {danger.min():.2f} et {danger.max():.2f} Hz")


if __name__ == "__main__":
    main()
