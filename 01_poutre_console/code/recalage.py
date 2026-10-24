"""Recale le jumeau numérique sur les données de l'essai de lâcher.

Étape 1 : FFT de l'accélération pour repérer la fréquence propre.
Étape 2 : moindres carrés sur un modèle oscillant amorti (fréquence, amortissement,
          amplitude, phase, décalage), à partir de t = 0,5 s pour laisser s'éteindre le 2e mode.
Étape 3 : on cherche le module d'Young effectif qui donne cette fréquence dans le
          modèle éléments finis (f proportionnelle à racine de E).

Sorties :
    resultats/parametres_recales.json
    resultats/recalage.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

import parametres as P
from fe_poutre import modes_propres

RACINE = Path(__file__).resolve().parent.parent
T_DEBUT = 0.5


def charger():
    d = np.loadtxt(RACINE / "donnees" / "essai_lacher.csv", delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1], d[:, 2]


def modele(p, t):
    f, zeta, A, phi, c = p
    w = 2 * np.pi * f
    return A * np.exp(-zeta * w * t) * np.cos(w * np.sqrt(1 - zeta**2) * t + phi) + c


def recaler(t, acc):
    # 1) FFT
    spectre = np.abs(np.fft.rfft(acc - acc.mean()))
    freqs = np.fft.rfftfreq(t.size, t[1] - t[0])
    f_fft = freqs[np.argmax(spectre[1:]) + 1]
    # 2) moindres carrés
    m = t >= T_DEBUT
    p0 = [f_fft, 0.01, acc[m].max(), 0.0, 0.0]
    res = least_squares(lambda p: modele(p, t[m]) - acc[m], p0,
                        bounds=([0.5 * f_fft, 0, 0, -np.pi, -5], [1.5 * f_fft, 0.2, 1e3, np.pi, 5]))
    f1, zeta, A, phi, c = res.x
    # 3) module d'Young effectif
    f_nom = modes_propres(E=P.E_nominal)[0][0]
    E_eff = P.E_nominal * (f1 / f_nom) ** 2
    residu = acc[m] - modele(res.x, t[m])
    return dict(f_fft_Hz=f_fft, f1_Hz=f1, amortissement=zeta, decalage_capteur=c,
                E_effectif_Pa=E_eff, f_nominale_Hz=f_nom,
                ecart_type_residu=float(residu.std())), res.x, (freqs, spectre)


def main():
    t, acc, _ = charger()
    resultat, p, (freqs, spectre) = recaler(t, acc)
    f2 = modes_propres(E=resultat["E_effectif_Pa"], n_modes=2)[0][1]

    reel_path = RACINE / "donnees" / "parametres_reels.json"
    if reel_path.exists():   # vérification seulement : le recalage ne s'en sert pas
        reel = json.loads(reel_path.read_text(encoding="utf-8"))
        # Note : le jumeau ignore la masse de colle et de vis ; E effectif l'absorbe en partie.
        resultat["erreur_E_pourcent"] = 100 * abs(resultat["E_effectif_Pa"] - reel["E_effectif_Pa"]) / reel["E_effectif_Pa"]
        resultat["erreur_amortissement_pourcent"] = 100 * abs(resultat["amortissement"] - reel["amortissement_mode1"]) / reel["amortissement_mode1"]

    (RACINE / "resultats").mkdir(exist_ok=True)
    (RACINE / "resultats" / "parametres_recales.json").write_text(
        json.dumps({k: float(v) for k, v in resultat.items()}, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), layout="constrained")
    ax[0].plot(t, acc, lw=0.6, color="0.6", label="mesure (accéléromètre)")
    ax[0].plot(t[t >= T_DEBUT], modele(p, t[t >= T_DEBUT]), lw=1.2, color="C3", label="jumeau recalé")
    ax[0].set(xlabel="temps [s]", ylabel="accélération en bout [m/s²]",
              title=f"Recalage : f1 = {resultat['f1_Hz']:.3f} Hz, amortissement = {resultat['amortissement']*100:.2f} %")
    ax[0].legend(loc="upper right")
    ax[1].semilogy(freqs, spectre, lw=0.8)
    ax[1].axvline(resultat["f1_Hz"], color="C3", ls="--", label=f"mode 1 : {resultat['f1_Hz']:.2f} Hz")
    ax[1].axvline(f2, color="C2", ls="--", label=f"mode 2 prédit : {f2:.1f} Hz")
    ax[1].set(xlim=(0, 80), ylim=(1, None), xlabel="fréquence [Hz]", ylabel="amplitude FFT")
    ax[1].legend()
    fig.savefig(RACINE / "resultats" / "recalage.png", dpi=130)

    print(f"Fréquence propre mesurée : {resultat['f1_Hz']:.3f} Hz (nominale : {resultat['f_nominale_Hz']:.3f} Hz)")
    print(f"Amortissement recalé     : {resultat['amortissement']*100:.2f} %")
    print(f"Module d'Young effectif  : {resultat['E_effectif_Pa']/1e9:.1f} GPa (nominal : {P.E_nominal/1e9:.0f} GPa)")
    if "erreur_E_pourcent" in resultat:
        print(f"Vérification : erreur sur E = {resultat['erreur_E_pourcent']:.2f} %, "
              f"sur l'amortissement = {resultat['erreur_amortissement_pourcent']:.1f} %")


if __name__ == "__main__":
    main()
