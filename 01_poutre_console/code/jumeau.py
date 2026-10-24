"""Jumeau numérique 3D de la poutre console.

Lance :  python code/jumeau.py
         python code/jumeau.py --capture     (enregistre resultats/jumeau.gif sans ouvrir de fenêtre)

À l'écran :
  - à gauche, le modèle Fusion (cao/poutre.stl) qui vibre, coloré selon la contrainte de flexion ;
  - à droite, l'accélération mesurée pendant l'essai et celle du jumeau, avec le curseur du temps ;
  - deux curseurs : masse ajoutée en bout, fréquence d'une excitation sinusoïdale en bout ;
  - un panneau de prédictions et d'alertes (résonance, plastification).

Le mouvement est montré au ralenti (x 4) pour qu'on puisse le suivre.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pyvista as pv

import parametres as P
from fe_poutre import modes_propres, evaluer, masse_modale, contrainte_flexion

RACINE = Path(__file__).resolve().parent.parent
RALENTI = 4.0
FORCE_EXCITATION = 0.15     # amplitude de la force sinusoïdale en bout [N]


# ---------------------------------------------------------------- géométrie
def charger_geometrie():
    """Charge le STL Fusion (mm) ; à défaut, construit la même géométrie."""
    chemin = RACINE / "cao" / "poutre.stl"
    if chemin.exists():
        maillage = pv.read(str(chemin))
        source = "cao/poutre.stl (Fusion)"
    else:
        lame = pv.Box(bounds=(0, 300, 0, 30, 0, 1))
        bloc = pv.Box(bounds=(280, 300, 0, 30, 1, 11))
        maillage = lame.merge(bloc)
        source = "géométrie de secours (cao/poutre.stl introuvable)"
    # Raffiner pour que la déformée soit lisse
    maillage = maillage.triangulate().subdivide_adaptive(max_edge_len=3.0)
    return maillage, source


def reperer_axes(maillage):
    """Trouve l'axe de la lame, l'axe vertical et le côté encastré, quelle que soit l'orientation."""
    pts = maillage.points
    etendues = pts.max(0) - pts.min(0)
    axe_long = int(np.argmax(etendues))
    autres = [a for a in range(3) if a != axe_long]
    axe_vertical = min(autres, key=lambda a: etendues[a])
    # le côté encastré est celui SANS le bloc : les points hauts (bloc) sont du côté du bout libre
    hauts = pts[:, axe_vertical] > pts[:, axe_vertical].min() + 0.6 * etendues[axe_vertical]
    s_min, s_max = pts[:, axe_long].min(), pts[:, axe_long].max()
    bout_libre_en_max = pts[hauts, axe_long].mean() > (s_min + s_max) / 2
    if bout_libre_en_max:
        s = pts[:, axe_long] - s_min
    else:
        s = s_max - pts[:, axe_long]
    return axe_long, axe_vertical, s / 1000.0 * (P.L / (etendues[axe_long] / 1000.0))


# ---------------------------------------------------------------- physique
class Jumeau:
    def __init__(self):
        rec = json.loads((RACINE / "resultats" / "parametres_recales.json").read_text(encoding="utf-8"))
        self.E = rec["E_effectif_Pa"]
        self.zeta = rec["amortissement"]
        self.f_mesure = rec["f1_Hz"]
        d = np.loadtxt(RACINE / "donnees" / "essai_lacher.csv", delimiter=",", skiprows=1)
        self.t_essai, self.acc_essai = d[:, 0], d[:, 1]
        self.fleche_lacher = self._fleche_initiale()
        self.masse_ajoutee = 0.0
        self.f_excitation = 0.0
        self.mettre_a_jour()

    def _fleche_initiale(self):
        """Flèche initiale de l'essai, déduite de l'amplitude d'accélération recalée."""
        w = 2 * np.pi * self.f_mesure
        m = self.t_essai >= 0.5
        a = self.acc_essai[m] - self.acc_essai[m].mean()
        amp = np.abs(a[:200]).max() / np.exp(-self.zeta * w * 0.5)
        return amp / w**2

    def mettre_a_jour(self):
        mb = P.masse_bloc + self.masse_ajoutee / 1000.0
        freqs, modes, self.le = modes_propres(E=self.E, masse_bout=mb, n_modes=1)
        self.f1 = freqs[0]
        self.mode = modes[0]
        self.k_modal = (2 * np.pi * self.f1) ** 2 * masse_modale(self.mode, self.E, mb)
        s = np.linspace(0, P.L, 400)
        _, courbure = evaluer(self.mode, self.le, s)
        self.sigma_par_metre = np.abs(contrainte_flexion(self.E, courbure)).max()  # [Pa / m de flèche]
        if self.f_excitation > 0:
            r = self.f_excitation / self.f1
            self.amplitude = FORCE_EXCITATION / self.k_modal / np.sqrt((1 - r**2) ** 2 + (2 * self.zeta * r) ** 2)
        else:
            self.amplitude = self.fleche_lacher
        self.sigma_max = self.sigma_par_metre * self.amplitude

    def fleche_bout(self, t):
        """Flèche en bout [m] à l'instant t (lâcher libre ou régime forcé établi)."""
        w = 2 * np.pi * self.f1
        if self.f_excitation > 0:
            return self.amplitude * np.cos(2 * np.pi * self.f_excitation * t)
        tt = t % P.duree_essai
        return self.fleche_lacher * np.exp(-self.zeta * w * tt) * np.cos(w * np.sqrt(1 - self.zeta**2) * tt)

    def alertes(self):
        msgs = []
        if self.f_excitation > 0 and abs(self.f_excitation / self.f1 - 1) < P.seuil_resonance:
            msgs.append(f"RÉSONANCE : excitation à {self.f_excitation:.2f} Hz, f propre = {self.f1:.2f} Hz")
        if self.sigma_max > P.sigma_e:
            msgs.append(f"PLASTIFICATION : sigma max = {self.sigma_max/1e6:.0f} MPa > {P.sigma_e/1e6:.0f} MPa")
        elif self.sigma_max > P.sigma_e / P.coef_securite:
            msgs.append(f"MARGE INSUFFISANTE : sigma max = {self.sigma_max/1e6:.0f} MPa "
                        f"> {P.sigma_e/P.coef_securite/1e6:.0f} MPa (coef. 1,5)")
        return msgs


# ---------------------------------------------------------------- affichage
def construire(capture=False, masse=0.0, excitation=0.0):
    maillage, source = charger_geometrie()
    axe_long, axe_vert, s = reperer_axes(maillage)
    j = Jumeau()
    pts0 = maillage.points.copy()

    def champs():
        w_mode, courbure = evaluer(j.mode, j.le, s)
        sigma_unit = np.abs(contrainte_flexion(j.E, courbure)) / 1e6   # MPa par m de flèche
        return w_mode, sigma_unit

    etat = {"t": 0.0, "w_mode": None, "sigma_unit": None, "ymax": 45.0}
    etat["w_mode"], etat["sigma_unit"] = champs()
    maillage.point_data["Contrainte [MPa]"] = np.zeros(maillage.n_points)

    pl = pv.Plotter(shape=(1, 2), col_weights=[3, 2], window_size=(1500, 780), off_screen=capture)
    pl.set_background("white")

    # --- vue 3D ---
    pl.subplot(0, 0)
    pl.add_mesh(maillage, scalars="Contrainte [MPa]", cmap="jet", clim=[0, P.sigma_e / 1e6],
                smooth_shading=True,
                scalar_bar_args={"title": "MPa", "color": "black", "vertical": True,
                                 "position_x": 0.90, "position_y": 0.45, "height": 0.35, "width": 0.07})
    socle = pv.Box(bounds=_socle(pts0, axe_long, axe_vert))
    pl.add_mesh(socle, color="dimgray")
    pl.add_text(f"Jumeau numérique - poutre console\nModèle : {source}", font_size=10, color="black")
    info = pl.add_text("", position="lower_right", font_size=9, color="black")
    alerte = pl.add_text("", position="lower_left", font_size=11, color="red")
    # caméra : vue de face, légèrement plongeante, la lame à l'horizontale
    centre = pts0.mean(0)
    devant = np.zeros(3)
    devant[3 - axe_long - axe_vert] = -1.0
    haut = np.zeros(3)
    haut[axe_vert] = 1.0
    pl.camera_position = [tuple(centre + 520 * devant + 200 * haut), tuple(centre + 20 * haut), tuple(haut)]
    pl.reset_camera()
    pl.camera.zoom(0.8)

    # --- courbes ---
    pl.subplot(0, 1)
    graphe = pv.Chart2D(x_label="temps [s]", y_label="accélération en bout [m/s²]")
    graphe.title = "Essai (gris) et jumeau (rouge)"
    graphe.line(j.t_essai, j.acc_essai, color="gray", width=1.0)
    w1 = 2 * np.pi * j.f1
    courbe_jumeau = graphe.line(j.t_essai, -w1**2 * j.fleche_bout(j.t_essai), color="red", width=1.5)
    curseur = graphe.line([0, 0], [-40, 40], color="blue", width=2.0)
    graphe.y_range = [-45, 45]
    pl.add_chart(graphe)

    def rafraichir_textes():
        info.SetText(1, (f"Paramètres recalés :\n  E effectif = {j.E/1e9:.1f} GPa\n"
                         f"  amortissement = {j.zeta*100:.2f} %\n\n"
                         f"Prédictions :\n  f propre = {j.f1:.2f} Hz\n"
                         f"  flèche en bout = {j.amplitude*1000:.1f} mm\n"
                         f"  contrainte max = {j.sigma_max/1e6:.0f} MPa\n"
                         f"  masse ajoutée = {j.masse_ajoutee:.0f} g\n"
                         f"  excitation = {'aucune (lâcher)' if j.f_excitation == 0 else f'{j.f_excitation:.2f} Hz'}"))
        msgs = j.alertes()
        alerte.SetText(0, "\n".join(msgs) if msgs else "")
        tt = j.t_essai
        acc_jumeau = -(2 * np.pi * (j.f_excitation or j.f1)) ** 2 * j.fleche_bout(tt)
        courbe_jumeau.update(tt, acc_jumeau)
        etat["ymax"] = max(45.0, 1.1 * float(np.abs(acc_jumeau).max()))
        graphe.y_range = [-etat["ymax"], etat["ymax"]]
        graphe.title = ("Essai (gris) et jumeau en lâcher libre (rouge)" if j.f_excitation == 0
                        else f"Essai (gris) et prédiction en régime forcé à {j.f_excitation:.2f} Hz (rouge)")

    def sur_masse(v):
        j.masse_ajoutee = float(v)
        j.mettre_a_jour()
        etat["w_mode"], etat["sigma_unit"] = champs()
        rafraichir_textes()

    def sur_excitation(v):
        j.f_excitation = 0.0 if v < 0.5 else float(v)
        j.mettre_a_jour()
        rafraichir_textes()

    pl.subplot(0, 0)
    pl.add_slider_widget(sur_masse, [0, 150], value=masse, title="Masse ajoutée en bout [g]",
                         pointa=(0.03, 0.82), pointb=(0.35, 0.82), style="modern", color="black",
                         title_height=0.02, fmt="%.0f")
    pl.add_slider_widget(sur_excitation, [0, 12], value=excitation, title="Excitation en bout [Hz] (0 = lâcher)",
                         pointa=(0.03, 0.67), pointb=(0.35, 0.67), style="modern", color="black",
                         title_height=0.02, fmt="%.2f")
    sur_masse(masse)
    sur_excitation(excitation)

    def pas(dt):
        etat["t"] += dt / RALENTI
        q = j.fleche_bout(etat["t"])
        nouveaux = pts0.copy()
        nouveaux[:, axe_vert] += etat["w_mode"] * q * 1000.0          # m -> mm
        maillage.points = nouveaux
        maillage.point_data["Contrainte [MPa]"] = etat["sigma_unit"] * abs(q)
        tc = etat["t"] % P.duree_essai
        curseur.update([tc, tc], [-etat["ymax"], etat["ymax"]])
        return q

    return pl, pas


def _socle(pts, axe_long, axe_vert):
    """Petit bloc gris qui figure l'encastrement, côté racine."""
    mn, mx = pts.min(0), pts.max(0)
    b = [mn[0], mx[0], mn[1], mx[1], mn[2], mx[2]]
    racine_en_min = True
    hauts = pts[:, axe_vert] > mn[axe_vert] + 0.6 * (mx[axe_vert] - mn[axe_vert])
    if pts[hauts, axe_long].mean() < (mn[axe_long] + mx[axe_long]) / 2:
        racine_en_min = False
    if racine_en_min:
        b[2 * axe_long], b[2 * axe_long + 1] = mn[axe_long] - 25, mn[axe_long]
    else:
        b[2 * axe_long], b[2 * axe_long + 1] = mx[axe_long], mx[axe_long] + 25
    b[2 * axe_vert], b[2 * axe_vert + 1] = mn[axe_vert] - 40, mn[axe_vert] + 40
    autre = 3 - axe_long - axe_vert
    b[2 * autre] -= 10
    b[2 * autre + 1] += 10
    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", action="store_true", help="enregistre resultats/jumeau.gif")
    ap.add_argument("--masse", type=float, default=0.0, help="masse ajoutée en bout [g] (capture)")
    ap.add_argument("--excitation", type=float, default=0.0, help="fréquence d'excitation [Hz] (capture)")
    ap.add_argument("--nom", default="jumeau", help="nom des fichiers de capture")
    args = ap.parse_args()
    pl, pas = construire(capture=args.capture, masse=args.masse, excitation=args.excitation)
    if args.capture:
        sortie = RACINE / "resultats" / f"{args.nom}.gif"
        pl.open_gif(str(sortie), fps=20)
        q_max = 0.0
        for _ in range(120):
            q = pas(1 / 20)
            pl.write_frame()
            if abs(q) >= q_max:          # image fixe à l'instant de flèche maximale
                q_max = abs(q)
                pl.screenshot(str(RACINE / "resultats" / f"{args.nom}.png"))
        pl.close()
        print(f"Capture enregistrée : {sortie}")
    else:
        pl.add_timer_event(max_steps=10**7, duration=33, callback=lambda _: (pas(0.033), pl.render()))
        pl.show()


if __name__ == "__main__":
    main()
