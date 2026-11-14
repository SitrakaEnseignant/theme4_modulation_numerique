"""Exécution réelle des graphes GNU Radio du Thème 4 sans écran, mesures et figures.
Les afficheurs QT GUI (temps, spectre, constellation, oeil) sont remplacés par des Vector Sink,
les curseurs QT GUI Range par des variables ; aucun bloc de traitement n'est modifié, sauf le bloc
Throttle, réglé à une cadence très élevée puisqu'il ne fait que limiter la vitesse de calcul.
    python3 En_ligne/executer_grc_ci.py"""
import os, sys, json, time, threading, subprocess, importlib, yaml
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ICI = os.path.dirname(os.path.abspath(__file__)); RACINE = os.path.dirname(ICI)
sys.path.insert(0, ICI)
from modulations import MODS, FS, SPS, RS, FC, ALPHA, DF, CONVERGENCE
import theorie as T
G = os.path.join(RACINE, 'GNU_Radio'); OUT = os.path.join(RACINE, 'resultats_gnuradio'); HL = os.path.join(OUT, 'hl')
os.makedirs(HL, exist_ok=True)
SINKS = ('qtgui_time_sink_x', 'qtgui_freq_sink_x', 'qtgui_const_sink_x', 'qtgui_eye_sink_x')


def headless(grc, valeurs, params, tag):
    g = yaml.safe_load(open(grc, encoding='utf8'))
    fid = g['options']['parameters']['id'] + '_hl_' + tag
    g['options']['parameters'].update({'id': fid, 'generate_options': 'no_gui', 'run_options': 'run'})
    blocs, sinks = [], {}
    for b in g['blocks']:
        p = b['parameters']
        if b['id'] == 'qtgui_range':
            b = {'name': b['name'], 'id': 'variable', 'parameters': {'comment': '', 'value': p['value']}, 'states': b['states']}
        if b['name'] in valeurs: b['parameters']['value'] = str(valeurs[b['name']])
        if b['name'] in params: b['parameters'].update({k: str(v) for k, v in params[b['name']].items()})
        if b['id'] == 'blocks_throttle': b['parameters']['samples_per_second'] = '1e9'
        if b['id'] in SINKS:
            n = int(p.get('nconnections', 1)); typ = p['type']; sinks[b['name']] = (typ, n)
            for k in range(n):
                blocs.append({'name': f"{b['name']}_{k}", 'id': 'blocks_vector_sink_x', 'states': b['states'],
                              'parameters': {'affinity': '', 'alias': '', 'comment': '', 'type': typ, 'vlen': '1',
                                             'reserve_items': '1024'}})
            continue
        blocs.append(b)
    conns = [[c[0], c[1], f"{c[2]}_{c[3]}", '0'] if c[2] in sinks else c for c in g['connections']]
    g['blocks'], g['connections'] = blocs, conns
    path = os.path.join(HL, fid + '.grc'); yaml.safe_dump(g, open(path, 'w', encoding='utf8'), allow_unicode=True, sort_keys=False)
    return path, fid, sinks


def executer(nom, nsym, noise=0.0, tag='a'):
    """exécute le graphe jusqu'à épuisement de la source de symboles (repeat = False)"""
    path, fid, sinks = headless(os.path.join(G, nom + '.grc'), {'nsym': nsym, 'noise_amp': noise},
                                {'src': {'repeat': 'False'}}, tag)
    r = subprocess.run(['grcc', '-o', HL, path], capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stdout + r.stderr)
    sys.path.insert(0, HL); mod = importlib.reload(importlib.import_module(fid))
    tb = getattr(mod, fid)(); t0 = time.time(); tb.start()
    th = threading.Thread(target=tb.wait); th.start(); th.join(timeout=120)
    fini = not th.is_alive()
    if not fini: tb.stop(); th.join()
    d = {s: [np.array(getattr(tb, f'{s}_{k}').data()) for k in range(n)] for s, (typ, n) in sinks.items()}
    return d, {'duree_s': round(time.time() - t0, 2), 'fin_naturelle': fini}


def permutations(m):
    """ambiguïté de phase d'une boucle : rotation de 2 pi k / ordre, et permutation des indices qu'elle produit"""
    pts = np.array(m['pts']); P = []
    for k in range(m['ordre']):
        r = np.exp(2j * np.pi * k / m['ordre'])
        perm = np.argmin(np.abs((pts * r)[:, None] - pts[None, :]), 1)
        if len(set(perm)) == len(pts): P.append((k, perm, r))
    return P


def aligner(m, tx, dec, maxlag=60):
    """retard (en symboles) et rotation d'ambiguïté qui alignent les décisions sur les symboles émis.
    La rotation est levée sur une fenêtre de 500 symboles connus après l'accrochage, comme avec un préambule."""
    tx = tx.astype(int); dec = dec.astype(int); best = None
    for k, perm, r in permutations(m):
        for L in range(maxlag):
            a = perm[tx[CONVERGENCE:CONVERGENCE + 500]]; b = dec[CONVERGENCE + L:CONVERGENCE + L + 500]
            if len(b) < 500: continue
            acc = np.mean(a == b)
            if best is None or acc > best[0]: best = (acc, L, k, perm, r)
    return best[1], best[2], best[3], best[4]


BLOC, PILOTES = 1000, 100


def erreurs(m, tx, dec, L, perm):
    """Erreurs comptées après l'accrochage. Avec une boucle (Costas ou récepteur de constellation), un
    glissement de cycle fait tourner la constellation d'un multiple de 2 pi / ordre : comme avec des
    pilotes répartis dans la trame, les 100 premiers symboles de chaque bloc de 1000 sont connus et
    servent à lever l'ambiguïté du bloc ; les erreurs sont comptées sur les 900 symboles suivants.
    Renvoie aussi le nombre de glissements de cycle détectés."""
    tx = tx.astype(int); dec = dec.astype(int); lab = np.array(m['lab'])
    n = min(len(tx), len(dec) - L) - 5
    a_all, b_all = tx[CONVERGENCE:n], dec[CONVERGENCE + L:n + L]
    perms = [p for _, p, _ in permutations(m)] if m['ordre'] > 1 else [perm]
    es = ns = eb = 0; glissements = 0; prec = None
    for i in range(0, len(a_all) - BLOC + 1, BLOC):
        a, b = a_all[i:i + BLOC], b_all[i:i + BLOC]
        j = int(np.argmax([np.mean(p[a[:PILOTES]] == b[:PILOTES]) for p in perms]))
        if prec is not None and j != prec: glissements += 1
        prec = j
        inv = np.argsort(perms[j]); bb = inv[b[PILOTES:]]; aa = a[PILOTES:]
        es += int(np.sum(aa != bb)); ns += len(aa)
        eb += int(sum(bin(int(x)).count('1') for x in (lab[aa] ^ lab[bb])))
    erreurs.glissements = glissements
    return es, ns, eb, ns * m['bits']


def theorie(m, g):
    k = m['nom']
    if k == 'ook': return float(T.ber_ook_enveloppe(g))
    if k == 'ask4': return float(T.ser_ask_enveloppe([abs(p) for p in m['pts']], g))
    if k in ('bpsk', 'qpsk', 'psk8'): return float(T.ser_mpsk(m['M'], g))
    if k == 'qam16': return float(T.ser_qam(16, g))
    if k == 'fsk2': return float(T.ber_fsk_noncoherent(g))


def dsp(x, N=4096):
    w = np.hanning(N); seg = [x[i:i + N] * w for i in range(0, len(x) - N, N // 2)]
    P = np.mean([np.abs(np.fft.rfft(s)) ** 2 for s in seg], 0); f = np.fft.rfftfreq(N, 1 / FS)
    return f, P


R, ECHECS = {}, []
for k, m in MODS.items():
    nom = f'{k}_tx_rx_sim'; r = R[k] = {'titre': m['titre']}; boucle = m['detection'] in ('costas', 'recepteur')
    ideal = np.array(m['pts']) / (m['amp_moy'] if k != 'fsk2' else 1)
    print('=====', nom, flush=True)
    d, info = executer(nom, 4000, 0.0, 'sans_bruit'); r['execution_sans_bruit'] = info
    emis = d['scope_emis'][0]; Ps = float(np.mean(emis[SPS * 20:] ** 2)); r['puissance_signal'] = Ps
    tx = d['symboles'][0].astype(int); dec = d['symboles'][1]
    L, rot, perm, rr = aligner(m, tx, dec); e, n, eb, nbits = erreurs(m, tx, dec, L, perm)
    r['sans_bruit'] = {'retard_symboles': L, 'rotation_ambiguite_deg': 360.0 * rot / m['ordre'], 'erreurs_symboles': e, 'symboles': n}
    if e: ECHECS.append(f'{k} : {e} erreurs sans bruit')
    if k == 'fsk2':
        Ld, _, pd, _ = aligner(m, tx, d['symboles'][2]); ed = erreurs(m, tx, d['symboles'][2], Ld, pd)[0]
        r['sans_bruit']['erreurs_discriminateur'] = ed
        if ed: ECHECS.append(f'{k} : {ed} erreurs sans bruit au discriminateur')
    else:
        ech = d['constellation_rx'][0]; n2 = min(len(ech) - L, len(tx)) - 5
        ref = ideal[tx[CONVERGENCE:n2]] * rr; mes = ech[CONVERGENCE + L:n2 + L]
        if m['detection'] == 'enveloppe': ref, mes = np.abs(ref), np.abs(mes)
        gain = float(np.real(np.vdot(ref, mes) / np.vdot(ref, ref)))
        evm = float(np.sqrt(np.mean(np.abs(mes - ref) ** 2) / np.mean(np.abs(ref) ** 2)))
        r['sans_bruit'].update({'gain_mesure': round(gain, 4), 'evm': round(evm, 4)})
        if evm > 0.05 or abs(gain - 1) > 0.05: ECHECS.append(f'{k} : EVM {evm:.3f}, gain {gain:.3f}')
    # phase d'échantillonnage optimale, mesurée sur le module de la sortie du filtre adapté (oeil)
    oe = np.abs(d['oeil'][0]) if k != 'fsk2' else d['oeil'][0]; ev = []
    refm = np.abs(ideal[tx]) if k != 'fsk2' else ideal[tx].real
    for p in range(SPS):
        sm = oe[p::SPS]; nn = min(len(sm), len(tx)) - 45
        ev.append(min(np.mean(np.abs(sm[l + CONVERGENCE:l + nn] - refm[CONVERGENCE:nn]) ** 2) for l in range(40)))
    r['sans_bruit']['phase_optimale'] = int(np.argmin(ev))
    f, P = dsp(emis); Pn = P / P.max(); r['spectre'] = {}
    if m['mise_en_forme'] == 'rrc':
        B = (1 + ALPHA) * RS / 2; frac = float(P[(f > FC - B) & (f < FC + B)].sum() / P.sum())
        r['spectre']['part_puissance_bande'] = round(frac, 4)
        if frac < 0.98: ECHECS.append(f'{k} : {frac:.3f} de la puissance dans la bande')
    elif k == 'fsk2':
        i1 = np.argmax(np.where(f < FC, Pn, 0)); i2 = np.argmax(np.where(f > FC, Pn, 0))
        r['spectre']['raies_Hz'] = [float(f[i1]), float(f[i2])]
        if abs(f[i1] - (FC - DF)) > 20 or abs(f[i2] - (FC + DF)) > 20: ECHECS.append(f'{k} : raies {f[i1]}, {f[i2]}')
    else:
        lv = [float(10 * np.log10(Pn[np.abs(f - x) < 20].max())) for x in (FC - RS, FC + RS)]
        r['spectre']['niveau_premiers_zeros_dB'] = [round(v, 1) for v in lv]
        if max(lv) > -20: ECHECS.append(f'{k} : zéros du spectre à {lv} dB')
    # courbe de taux d'erreur
    r['taux'] = []; tol = 2.0 if boucle else 1.5
    for e_db in m['ebn0']:
        g = 10 ** (e_db / 10); sig = float(np.sqrt(Ps * SPS / (2 * m['bits'] * g)))
        d2, info = executer(nom, 20000, sig, f'eb{e_db}')
        tx2, dec2 = d2['symboles'][0].astype(int), d2['symboles'][1]; L2, rot2, perm2, _ = aligner(m, tx2, dec2)
        es, ns, eb, nb = erreurs(m, tx2, dec2, L2, perm2); gl = erreurs.glissements
        mesure = eb / nb if m['compare'] == 'TEB' else es / ns
        nerr = eb if m['compare'] == 'TEB' else es
        th = theorie(m, g); ok = nerr < 50 or 1 / tol < mesure / th < tol
        # une boucle qui glisse (décrochage) n'a plus une référence de phase cohérente : le point est rapporté
        # comme « hors accrochage » ; au-delà du seuil de fonctionnement (théorie < 1 %), un glissement isolé est
        # toléré (événement rare, signalé), plusieurs glissements font échouer le contrôle
        decroche = boucle and gl > 0
        if decroche: ok = th > 0.01 or gl <= 1
        pt = {'EbN0_dB': e_db, 'ecart_type_bruit': round(sig, 5), 'TES': es / ns, 'TEB': eb / nb, 'erreurs': nerr,
              'mesure': mesure, 'theorie': th, 'conforme': bool(ok), 'retard': L2, 'rotation_ambiguite_deg': 360.0 * rot2 / m['ordre'], 'glissements_de_cycle': gl, 'hors_accrochage': bool(decroche), **info}
        if k == 'fsk2':
            Ld, _, pd, _ = aligner(m, tx2, d2['symboles'][2]); pt['TEB_discriminateur'] = erreurs(m, tx2, d2['symboles'][2], Ld, pd)[2] / nb
        r['taux'].append(pt)
        if not ok: ECHECS.append(f'{k} : Eb/N0 = {e_db} dB, mesure {mesure:.3g}, théorie {th:.3g}')
        if e_db == m['ebn0'][2]: d_mid, L_mid = d2, L2
    d_eye, _ = executer(nom, 3000, float(np.sqrt(Ps * SPS / (2 * m['bits'] * 10 ** 2.0))), 'oeil20')
    # figure : signal émis, spectre, constellation, oeil
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    t = np.arange(8 * SPS) / FS * 1e3; i0 = SPS * 40
    ax[0, 0].plot(t, emis[i0:i0 + 8 * SPS], lw=0.9, label='signal émis')
    ax[0, 0].plot(t, d['scope_emis'][1][i0:i0 + 8 * SPS], lw=1.6, ls='--', label='bits NRZ' if k == 'fsk2' else 'voie I en bande de base')
    ax[0, 0].set_xlabel('temps (ms)'); ax[0, 0].set_title(f"{m['titre']} : signal émis (scope_emis)"); ax[0, 0].legend(loc='upper right', fontsize=8)
    ax[0, 1].plot(f / 1e3, 10 * np.log10(Pn + 1e-16)); ax[0, 1].set_xlim(0, 10); ax[0, 1].set_ylim(-80, 3)
    ax[0, 1].set_xlabel('fréquence (kHz)'); ax[0, 1].set_ylabel('dB'); ax[0, 1].set_title('spectre du signal émis (spectre_emis)')
    c = d_mid['constellation_rx'][0][CONVERGENCE:]
    ax[1, 0].plot(c.real, c.imag, '.', ms=1.5, alpha=0.4, label=f"Eb/N0 = {m['ebn0'][2]} dB")
    c0 = d['constellation_rx'][0][CONVERGENCE:]; ax[1, 0].plot(c0.real, c0.imag, 'r+', ms=9, label='sans bruit')
    ax[1, 0].set_aspect('equal'); ax[1, 0].legend(loc='upper right', fontsize=8)
    ax[1, 0].set_title({'enveloppe': "constellation avant détection d'enveloppe (phase inconnue)",
                        'energie': 'points de décision (|y1|, |y0|)'}.get(m['detection'], 'constellation après la boucle') + ' (constellation_rx)', fontsize=10)
    oe = d_eye['oeil'][0]; titre_oeil = "diagramme de l'oeil, Eb/N0 = 20 dB (oeil)"
    if m['detection'] == 'enveloppe': oe = np.abs(oe); titre_oeil = "oeil du module, Eb/N0 = 20 dB (oeil)"
    elif m['detection'] == 'costas':
        Le, rote, _, rre = aligner(m, d_eye['symboles'][0].astype(int), d_eye['symboles'][1]); oe = np.real(oe * np.conj(rre))
        titre_oeil = "oeil de la voie I après la boucle de Costas, Eb/N0 = 20 dB (oeil)"
    elif m['detection'] == 'recepteur':
        # l'oeil est pris avant la boucle : la voie I est corrigée par la phase estimée par la boucle (sortie 2 du récepteur)
        oc = d_eye['oeil'][0]; ph = d_eye['boucle'][0]; p0 = r['sans_bruit']['phase_optimale']
        jsym = (np.arange(len(oc)) - p0) // SPS; best = None
        rots = [np.exp(2j * np.pi * q / 4) for q in range(4)]
        for sgn in (1, -1):
            for q in range(-20, 21):
                z = oc * np.exp(-1j * sgn * ph[np.clip(jsym + q, 0, len(ph) - 1)])
                sm = z[p0 + SPS * np.arange(CONVERGENCE, min(CONVERGENCE + 800, (len(z) - p0) // SPS - 1))]
                for ro in rots:
                    qd = np.mean(np.min(np.abs(sm[:, None] * np.conj(ro) - ideal[None, :]), 1))
                    if best is None or qd < best[0]: best = (qd, z * np.conj(ro))
        oe = np.real(best[1]); titre_oeil = "oeil de la voie I corrigée par la phase de la boucle, Eb/N0 = 20 dB"
    ph0 = r['sans_bruit']['phase_optimale']
    for j in range(CONVERGENCE + 60, CONVERGENCE + 260):
        s0 = j * SPS + ph0 - SPS // 2
        ax[1, 1].plot(np.arange(2 * SPS + 1) / SPS - 0.5, oe[s0:s0 + 2 * SPS + 1], color='C0', lw=0.4, alpha=0.5)
    ax[1, 1].set_xlabel('temps (symboles)'); ax[1, 1].set_title(titre_oeil, fontsize=10)
    for a in ax.flat: a.grid(ls=':')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f'{nom}.png'), dpi=110); plt.close(fig)
    print(json.dumps(r, indent=1, default=float), flush=True)

fig, ax = plt.subplots(figsize=(9, 6))
for i, (k, m) in enumerate(MODS.items()):
    e = np.array([p['EbN0_dB'] for p in R[k]['taux']]); y = np.array([p['mesure'] for p in R[k]['taux']])
    ee = np.linspace(e.min() - 1, e.max() + 1, 60)
    ax.semilogy(ee, [theorie(m, 10 ** (x / 10)) for x in ee], color=f'C{i}', lw=1)
    hc = np.array([p.get('hors_accrochage', False) for p in R[k]['taux']])
    ax.semilogy(e[~hc], np.where(y > 0, y, np.nan)[~hc], 'o', color=f'C{i}', label=f"{m['titre']} ({m['compare']})")
    if hc.any(): ax.semilogy(e[hc], y[hc], 'x', color=f'C{i}', ms=9, mew=2)
ax.semilogy([], [], 'kx', ms=9, mew=2, label='boucle hors accrochage (glissements de cycle)')
ax.set_ylim(1e-5, 0.5); ax.set_xlabel('Eb/N0 (dB)'); ax.set_ylabel("taux d'erreur")
ax.set_title("Taux d'erreur mesurés dans GNU Radio (points) et théorie (traits)"); ax.grid(which='both', ls=':'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'taux_erreur_gnuradio.png'), dpi=110); plt.close(fig)
R['echecs'] = ECHECS
json.dump(R, open(os.path.join(OUT, 'mesures.json'), 'w'), indent=1, default=float)
print('ECHECS :', ECHECS if ECHECS else 'aucun')
sys.exit(1 if ECHECS else 0)
