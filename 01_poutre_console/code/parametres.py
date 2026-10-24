"""Paramètres nominaux du système : lame d'acier encastrée avec masse en bout.

Unités SI partout (m, kg, s, Pa), sauf indication contraire.
Les dimensions correspondent au modèle Fusion de cao/poutre.stl.
"""

# --- Géométrie de la lame (partie flexible) ---
L = 0.300        # longueur libre [m] (encastrement en x = 0)
b = 0.030        # largeur [m]
h = 0.001        # épaisseur [m]

# --- Masse en bout (bloc d'acier 20 x 30 x 10 mm posé sur la lame) ---
bloc_longueur = 0.020
bloc_largeur = 0.030
bloc_hauteur = 0.010

# --- Matériau : acier de construction S235 ---
E_nominal = 210e9      # module d'Young [Pa]
rho = 7850.0           # masse volumique [kg/m3]
sigma_e = 235e6        # limite d'élasticité [Pa]

# --- Grandeurs dérivées ---
I = b * h**3 / 12                                          # moment quadratique [m4]
masse_lame = rho * L * b * h                               # [kg]
masse_bloc = rho * bloc_longueur * bloc_largeur * bloc_hauteur  # [kg]

# --- Essai de référence : lâcher depuis une flèche imposée en bout ---
fleche_initiale = 0.030      # [m]
duree_essai = 5.0            # [s]
frequence_echantillonnage = 1000.0   # [Hz]
position_jauge = 0.010       # jauge de déformation à 10 mm de l'encastrement [m]

# --- Seuils d'alerte ---
seuil_resonance = 0.10       # alerte si |f_excitation / f_propre - 1| < 10 %
coef_securite = 1.5          # alerte plastification si sigma_max > sigma_e / 1,5
