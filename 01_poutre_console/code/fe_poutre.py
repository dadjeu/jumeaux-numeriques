"""Modèle éléments finis d'une poutre d'Euler-Bernoulli encastrée avec masse en bout.

Même démarche que le TP de méthode des éléments finis (matrices élémentaires,
assemblage, résolution), étendue aux éléments de poutre d'Hermite :
2 degrés de liberté par nœud (flèche w, rotation theta).
"""
import numpy as np
from scipy.linalg import eigh

import parametres as P


def matrices_elementaires(EI, rhoA, le):
    """Matrices de rigidité et de masse (cohérente) d'un élément de poutre d'Hermite."""
    k = EI / le**3 * np.array([
        [12, 6 * le, -12, 6 * le],
        [6 * le, 4 * le**2, -6 * le, 2 * le**2],
        [-12, -6 * le, 12, -6 * le],
        [6 * le, 2 * le**2, -6 * le, 4 * le**2],
    ])
    m = rhoA * le / 420 * np.array([
        [156, 22 * le, 54, -13 * le],
        [22 * le, 4 * le**2, 13 * le, -3 * le**2],
        [54, 13 * le, 156, -22 * le],
        [-13 * le, -3 * le**2, -22 * le, 4 * le**2],
    ])
    return k, m


def assembler(E, masse_bout, n_elem=40):
    """Assemble K et M globales ; encastrement en x = 0, masse ponctuelle au dernier nœud."""
    le = P.L / n_elem
    n_ddl = 2 * (n_elem + 1)
    K = np.zeros((n_ddl, n_ddl))
    M = np.zeros((n_ddl, n_ddl))
    k_e, m_e = matrices_elementaires(E * P.I, P.rho * P.b * P.h, le)
    for e in range(n_elem):
        ddl = slice(2 * e, 2 * e + 4)
        K[ddl, ddl] += k_e
        M[ddl, ddl] += m_e
    M[-2, -2] += masse_bout          # masse ponctuelle sur la flèche du dernier nœud
    libres = np.arange(2, n_ddl)      # on retire w et theta du nœud encastré
    return K[np.ix_(libres, libres)], M[np.ix_(libres, libres)], le


def modes_propres(E=P.E_nominal, masse_bout=P.masse_bloc, n_modes=3, n_elem=40):
    """Fréquences propres [Hz] et déformées modales (normalisées à 1 en bout)."""
    K, M, le = assembler(E, masse_bout, n_elem)
    valeurs, vecteurs = eigh(K, M, subset_by_index=[0, n_modes - 1])
    freqs = np.sqrt(valeurs) / (2 * np.pi)
    modes = []
    for j in range(n_modes):
        q = np.concatenate([[0.0, 0.0], vecteurs[:, j]])   # on remet l'encastrement
        q /= q[-2]                                           # flèche en bout = 1
        modes.append(q)
    return freqs, modes, le


def evaluer(q, le, x):
    """Flèche w(x) et courbure w''(x) d'une solution nodale q (interpolation d'Hermite)."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, P.L)
    n_elem = (len(q) // 2) - 1
    e = np.minimum((x / le).astype(int), n_elem - 1)
    xi = x / le - e                                    # coordonnée locale dans [0, 1]
    w1, t1, w2, t2 = (q[2 * e], q[2 * e + 1], q[2 * e + 2], q[2 * e + 3])
    N1 = 1 - 3 * xi**2 + 2 * xi**3
    N2 = le * (xi - 2 * xi**2 + xi**3)
    N3 = 3 * xi**2 - 2 * xi**3
    N4 = le * (-xi**2 + xi**3)
    w = N1 * w1 + N2 * t1 + N3 * w2 + N4 * t2
    d2N1 = (-6 + 12 * xi) / le**2
    d2N2 = (-4 + 6 * xi) / le
    d2N3 = (6 - 12 * xi) / le**2
    d2N4 = (-2 + 6 * xi) / le
    courbure = d2N1 * w1 + d2N2 * t1 + d2N3 * w2 + d2N4 * t2
    return w, courbure


def masse_modale(q, E=P.E_nominal, masse_bout=P.masse_bloc, n_elem=40):
    """Masse généralisée q^T M q du mode q (normalisé à 1 en bout) [kg]."""
    _, M, _ = assembler(E, masse_bout, n_elem)
    ql = q[2:]
    return float(ql @ M @ ql)


def contrainte_flexion(E, courbure):
    """Contrainte normale de flexion en peau de la lame : sigma = E * (h/2) * w''."""
    return E * (P.h / 2) * courbure


if __name__ == "__main__":
    f, modes, le = modes_propres()
    k = 3 * P.E_nominal * P.I / P.L**3
    m_eff = P.masse_bloc + 0.2357 * P.masse_lame
    f_analytique = np.sqrt(k / m_eff) / (2 * np.pi)
    print(f"Masse lame = {P.masse_lame*1e3:.1f} g, masse bloc = {P.masse_bloc*1e3:.1f} g")
    print(f"Fréquences propres EF : {', '.join(f'{v:.2f}' for v in f)} Hz")
    print(f"Mode 1 analytique (Rayleigh) : {f_analytique:.2f} Hz "
          f"-> écart {abs(f[0] - f_analytique) / f_analytique * 100:.2f} %")
