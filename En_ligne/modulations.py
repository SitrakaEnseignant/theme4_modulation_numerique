"""Description unique des modulations numériques du Thème 4.
Utilisée par outils/generer_grc.py (graphes GNU Radio) et par En_ligne/executer_grc_ci.py (mesures).
Les mêmes valeurs figurent dans Simulink/construire_et_verifier_modeles.m."""
import numpy as np

FS = 32000.0        # cadence d'échantillonnage (Hz)
SPS = 32            # échantillons par symbole
RS = FS / SPS       # rapidité de modulation : 1 000 symboles par seconde
FC = 4000.0         # fréquence de la porteuse (Hz)
ALPHA = 0.35        # facteur de retombée du filtre en racine de cosinus surélevé
NTAPS = 11 * SPS + 1  # 353 coefficients : retard de 176 échantillons par filtre, 11 symboles pour la paire
DF = 1000.0         # excursion de la FSK binaire (Hz) : tons à fc - 1000 et fc + 1000 Hz, indice h = 2 (comme la leçon Learn SDR)
PHI_OL = 0.7        # erreur de phase de l'oscillateur local du récepteur (rad)
DF_OL = 1.0         # erreur de fréquence de l'oscillateur local du récepteur (Hz)
W_COSTAS = '2*math.pi/12800'   # bande de boucle de Costas, par échantillon (cadence complète)
W_RECEPTEUR = '2*math.pi/400'  # bande de boucle du récepteur de constellation, par symbole
CONVERGENCE = 300   # symboles ignorés pendant l'accrochage des boucles


def gray(k):
    return k ^ (k >> 1)


def _psk(M, rot=0.0):
    # le point d'indice k porte l'étiquette binaire gray(k) : deux points voisins diffèrent d'un bit
    return [complex(np.round(np.cos(rot + 2 * np.pi * k / M), 12), np.round(np.sin(rot + 2 * np.pi * k / M), 12)) for k in range(M)], \
           [gray(k) for k in range(M)]


def _qam16():
    niv = [-1.0, -1 / 3, 1 / 3, 1.0]          # niveaux -3, -1, 1, 3 divisés par 3
    pts, lab = [], []
    for i in range(4):
        for q in range(4):
            pts.append(complex(niv[i], niv[q])); lab.append((gray(i) << 2) | gray(q))
    return pts, lab


MODS = {
    'ook':   dict(titre='OOK (tout ou rien)', M=2, pts=[0j, 1 + 0j], lab=[0, 1], mise_en_forme='rect',
                  detection='enveloppe', ebn0=[6, 8, 10, 12, 14], compare='TEB'),
    'ask4':  dict(titre='ASK à 4 niveaux', M=4, pts=[0j, 1 / 3 + 0j, 2 / 3 + 0j, 1 + 0j], lab=[gray(k) for k in range(4)],
                  mise_en_forme='rect', detection='enveloppe', ebn0=[10, 12, 14, 16, 18], compare='TES'),
    'bpsk':  dict(titre='BPSK', M=2, pts=[1 + 0j, -1 + 0j], lab=[0, 1], mise_en_forme='rect', detection='costas', ordre=2,
                  ebn0=[0, 2, 4, 6, 8], compare='TES'),
    'qpsk':  dict(titre='QPSK', M=4, mise_en_forme='rect', detection='costas', ordre=4, ebn0=[0, 2, 4, 6, 8], compare='TES'),
    'psk8':  dict(titre='8-PSK', M=8, mise_en_forme='rrc', detection='costas', ordre=8, ebn0=[4, 6, 8, 10, 12], compare='TES'),
    'qam16': dict(titre='16-QAM', M=16, mise_en_forme='rrc', detection='recepteur', ordre=4, ebn0=[4, 6, 8, 10, 12], compare='TES'),
    'fsk2':  dict(titre='FSK binaire', M=2, mise_en_forme='rect', detection='energie', ebn0=[4, 6, 8, 10, 12], compare='TEB'),
}
MODS['qpsk']['pts'], MODS['qpsk']['lab'] = _psk(4, np.pi / 4)
MODS['psk8']['pts'], MODS['psk8']['lab'] = _psk(8, np.pi / 8)
MODS['qam16']['pts'], MODS['qam16']['lab'] = _qam16()
MODS['fsk2']['pts'], MODS['fsk2']['lab'] = [-1 + 0j, 1 + 0j], [0, 1]
for _k, _m in MODS.items():
    _m['nom'] = _k
    _m['bits'] = int(np.log2(_m['M']))
    _m['amp_moy'] = float(np.mean(np.abs(_m['pts'])))
    _m.setdefault('ordre', 1)


def table(pts):
    """liste Python des points, écrite pour un champ de GRC"""
    return '[' + ', '.join(f'complex({p.real!r}, {p.imag!r})' for p in pts) + ']'
