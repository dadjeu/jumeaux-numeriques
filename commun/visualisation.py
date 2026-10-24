"""Outils communs : chargement du modèle 3D Fusion et affichage coloré d'une grandeur."""
from pathlib import Path

import numpy as np
import pyvista as pv


def charger_stl(chemin):
    """Charge un export STL de Fusion et renvoie un maillage PyVista (unités : mm)."""
    chemin = Path(chemin)
    if not chemin.exists():
        raise FileNotFoundError(f"Modèle introuvable : {chemin} — exporte-le depuis Fusion en STL.")
    return pv.read(str(chemin))


def afficher_champ(maillage, valeurs, nom="grandeur", unite="", cmap="coolwarm"):
    """Affiche le maillage coloré par une grandeur définie en chaque point."""
    maillage = maillage.copy()
    maillage.point_data[nom] = np.asarray(valeurs)
    p = pv.Plotter()
    p.add_mesh(maillage, scalars=nom, cmap=cmap,
               scalar_bar_args={"title": f"{nom} {unite}".strip()})
    p.add_axes()
    p.show()
