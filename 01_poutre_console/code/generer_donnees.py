"""Génère le jeu de données SYNTHÉTIQUE de l'essai de lâcher.

On simule le « vrai » système, dont les paramètres diffèrent un peu des valeurs
nominales (encastrement imparfait, amortissement inconnu), puis on ajoute ce que
feraient de vrais capteurs : bruit, décalage de zéro, échantillonnage à 1 kHz.

Le jumeau ne connaît PAS ces paramètres : il doit les retrouver par recalage.
Ils sont enregistrés à part (donnees/parametres_reels.json) pour vérifier le recalage.

Sorties :
    donnees/essai_lacher.csv
    donnees/parametres_reels.json
"""
import json
from pathlib import Path

import numpy as np

import parametres as P
from fe_poutre import modes_propres, evaluer

DOSSIER = Path(__file__).resolve().parent.parent / "donnees"

# --- « Vérité » cachée du système réel ---
REEL = {
    "E_effectif_Pa": 196e9,       # rigidité un peu plus faible : encastrement imparfait
    "amortissement_mode1": 0.012,
    "amortissement_mode2": 0.020,
    "fleche_initiale_m": 0.0292,  # l'opérateur a écarté de 29,2 mm au lieu de 30
    "amplitude_mode2_m": 0.0002,  # petite part du 2e mode excitée au lâcher
    "masse_bout_kg": 0.0486,      # colle et vis ajoutent 1,5 g au bloc (inconnu du jumeau)
}

# --- Caractéristiques des capteurs ---
BRUIT_ACCELERO = 0.08      # écart-type [m/s2]
DECALAGE_ACCELERO = 0.03   # décalage de zéro [m/s2]
BRUIT_JAUGE = 2.0          # écart-type [microdéformations]


def simuler(t, reel):
    f, modes, le = modes_propres(E=reel["E_effectif_Pa"], masse_bout=reel["masse_bout_kg"], n_modes=2)
    amplitudes = [reel["fleche_initiale_m"], reel["amplitude_mode2_m"]]
    zetas = [reel["amortissement_mode1"], reel["amortissement_mode2"]]
    acc = np.zeros_like(t)
    deformation = np.zeros_like(t)
    for fj, qj, Aj, zj in zip(f, modes, amplitudes, zetas):
        w = 2 * np.pi * fj
        wd = w * np.sqrt(1 - zj**2)
        x = Aj * np.exp(-zj * w * t) * np.cos(wd * t)     # flèche en bout du mode j
        acc += -w**2 * x                                   # accélération (approx. faible amortissement)
        _, courbure = evaluer(qj, le, P.position_jauge)
        deformation += (P.h / 2) * courbure * x            # déformation en peau à la jauge
    return acc, deformation * 1e6, f


def main():
    rng = np.random.default_rng(2026)
    t = np.arange(0, P.duree_essai, 1 / P.frequence_echantillonnage)
    acc, micro_def, f = simuler(t, REEL)
    acc_mesure = acc + DECALAGE_ACCELERO + rng.normal(0, BRUIT_ACCELERO, t.size)
    def_mesure = micro_def + rng.normal(0, BRUIT_JAUGE, t.size)

    DOSSIER.mkdir(exist_ok=True)
    entete = "temps_s,acceleration_bout_m_s2,deformation_jauge_microdef"
    np.savetxt(DOSSIER / "essai_lacher.csv", np.column_stack([t, acc_mesure, def_mesure]),
               delimiter=",", header=entete, comments="", fmt=["%.4f", "%.4f", "%.2f"])
    infos = dict(REEL, frequences_reelles_Hz=[round(v, 4) for v in f],
                 note="Données synthétiques générées par code/generer_donnees.py")
    (DOSSIER / "parametres_reels.json").write_text(json.dumps(infos, indent=2), encoding="utf-8")
    print(f"{t.size} échantillons écrits dans {DOSSIER / 'essai_lacher.csv'}")


if __name__ == "__main__":
    main()
