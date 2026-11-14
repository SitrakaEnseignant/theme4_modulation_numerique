"""Probabilités d'erreur théoriques (canal à bruit blanc additif gaussien), en fonction de Eb/N0 linéaire."""
import numpy as np
from scipy.special import erfc
from scipy.integrate import quad
from scipy.stats import ncx2
Q = lambda x: 0.5 * erfc(np.asarray(x) / np.sqrt(2))

def ser_mpsk(M, g):
    """SER exact de la M-PSK cohérente (formule de Craig), g = Eb/N0"""
    gs = g * np.log2(M); a = np.sin(np.pi / M) ** 2
    return quad(lambda p: np.exp(-gs * a / np.sin(p) ** 2), 0, (M - 1) * np.pi / M)[0] / np.pi

def ser_qam(M, g):
    gs = g * np.log2(M); P = 2 * (1 - 1 / np.sqrt(M)) * Q(np.sqrt(3 * gs / (M - 1)))
    return 1 - (1 - P) ** 2

def ser_ask_unipolaire(M, g):
    """M-ASK unipolaire {0, 1, ..., M-1}·d, détection cohérente, seuils à mi-distance"""
    pts = np.arange(M); E = np.mean(pts ** 2); gs = g * np.log2(M)
    x = (1 / 2) * np.sqrt(2 * gs / E)          # d/(2 sigma) avec d = 1
    return 2 * (M - 1) / M * Q(x)

def ber_ook_enveloppe(g):
    """OOK {0, A}, détection d'enveloppe, seuil A/2 : P = 1/2 [exp(-g/2) + 1 - Q1(2 sqrt g, sqrt g)]"""
    a2, b2 = 4 * g, g
    q1 = ncx2.sf(b2, 2, a2)
    return 0.5 * (np.exp(-b2 / 2) + 1 - q1)

def ber_fsk_noncoherent(g): return 0.5 * np.exp(-g / 2)
def ber_fsk_coherent(g): return Q(np.sqrt(g))
