# Jumeaux numériques de systèmes mécaniques

Série de 10 jumeaux numériques, un par semaine : un modèle 3D conçu dans **Autodesk Fusion**,
un essai de référence, un jeu de données d'essai, le recalage du modèle, une **simulation 3D animée**
en Python, des prédictions et la détection de cas critiques.

Auteur : Junior Teddy Dadjeu Temgoua — élève ingénieur, ENSTA Paris (majeure mécanique).

> Les jeux de données sont **synthétiques** (générés à partir de la physique, avec bruit de mesure),
> sauf mention « données mesurées » dans le README du projet.

## Projets

| N° | Système | Essai de référence | Grandeur visualisée | Cas critique surveillé | Statut |
| --- | --- | --- | --- | --- | --- |
| 01 | [Poutre console en acier](01_poutre_console/) | Vibration libre (lâcher) | Déformée animée, contrainte de Von Mises | Résonance, plastification | À faire |
| 02 | [Cuve cylindrique avec orifice de vidange](02_cuve_vidange/) | Vidange chronométrée | Niveau d'eau, vitesse de jet | Prise d'air (vortex) en fin de vidange | À faire |
| 03 | [Presse hydraulique de l'atelier](03_presse_hydraulique/) | Compression d'une éprouvette | Pression, contrainte dans le bâti et la tige | Surpression, flambement, fuite | À faire |
| 04 | [Dissipateur thermique à ailettes](04_dissipateur_thermique/) | Chauffe à puissance constante | Champ de température 3D | Surchauffe du composant | À faire |
| 05 | [Échangeur tubulaire eau chaude / eau froide](05_echangeur_tubulaire/) | Essai à deux débits | Température le long des tubes | Encrassement | À faire |
| 06 | [Conduite avec coude et vanne](06_conduite_vanne/) | Perte de charge en fonction du débit | Pression et vitesse dans la conduite | Cavitation, coup de bélier | À faire |
| 07 | [Rotor avec disque sur deux paliers](07_rotor_paliers/) | Montée en vitesse | Orbite de l'arbre, amplitude de vibration | Vitesse critique, balourd | À faire |
| 08 | [Pompe centrifuge](08_pompe_centrifuge/) | Courbe hauteur-débit | Hauteur, rendement, NPSH | Cavitation, hors plage | À faire |
| 09 | [Pile de pont dans un courant](09_pile_de_pont/) | Essai en canal à plusieurs vitesses | Sillage, force, vibration | Accrochage fréquentiel (lock-in) | À faire |
| 10 | [Petite éolienne tripale](10_eolienne/) | Essai en soufflerie | Puissance, rotation, efforts sur les pales | Survitesse, décrochage, surcharge du mât | À faire |

## Structure de chaque projet

```
NN_nom_du_projet/
├── cao/          modèle Fusion (.f3d) et exports (.stl pour Python, .step pour COMSOL / Ansys)
├── donnees/      jeu de données d'essai (.csv)
├── code/         scripts Python : generer_donnees.py, recalage.py, jumeau.py
├── resultats/    figures, captures et vidéo de la simulation
└── README.md     objectif, modèle, résultats chiffrés, limites
```

Le dossier `commun/` contient les fonctions partagées par tous les projets (lecture du STL, affichage 3D, alertes).

## Installation

```
python -m venv .venv
.venv\Scripts\activate        (Windows)   |   source .venv/bin/activate   (Mac / Linux)
pip install -r requirements.txt
```

Puis, dans un projet : `python code/jumeau.py`
