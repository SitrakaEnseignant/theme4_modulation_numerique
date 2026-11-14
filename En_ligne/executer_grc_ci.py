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
from modulations import MODS, FS, SPS, RS, FC, ALPHA, DF
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


def aligner(tx, dec, maxlag=60):
    """retard (en symboles) entre les décisions et les symboles émis"""
    n = min(len(tx), len(dec)) - maxlag - 10
    best = max(range(maxlag), key=lambda L: np.mean(dec[L:L + n] == tx[:n]))
    return best


def erreurs(m, tx, dec, L):
    n = min(len(tx), len(dec) - L) - 20
    a, b = tx[20:n].astype(int), dec[20 + L:n + L].astype(int)
    lab = np.array(m['lab'])
    nb = int(sum(bin(int(x)).count('1') for x in (lab[a] ^ lab[b])))
    return int(np.sum(a != b)), len(a), nb, len(a) * m['bits']


def theorie(m, g):
    k = m['nom']
    if k == 'ook': return float(T.ber_ook_enveloppe(g))
    if k == 'ask4': return float(T.ser_ask_unipolaire(4, g))
    if k in ('bpsk', 'qpsk', 'psk8'): return float(T.ser_mpsk(m['M'], g))
    if k == 'qam16': return float(T.ser_qam(16, g))
    if k == 'fsk2': return float(T.ber_fsk_noncoherent(g))


def dsp(x, N=4096):
    w = np.hanning(N); seg = [x[i:i + N] * w for i in range(0, len(x) - N, N // 2)]
    P = np.mean([np.abs(np.fft.rfft(s)) ** 2 for s in seg], 0); f = np.fft.rfftfreq(N, 1 / FS)
    return f, P


R, ECHECS = {}, []
for k, m in MODS.items():
    nom = f'{k}_tx_rx_sim'; r = R[k] = {'titre': m['titre']}
    print('=====', nom, flush=True)
    d, info = executer(nom, 4000, 0.0, 'sans_bruit'); r['execution_sans_bruit'] = info
    emis = d['scope_emis'][0]; Ps = float(np.mean(emis[SPS * 20:] ** 2)); r['puissance_signal'] = Ps
    tx, dec = d['symboles'][0], d['symboles'][1]
    L = aligner(tx, dec); e, n, eb, nbits = erreurs(m, tx, dec, L)
    r['sans_bruit'] = {'retard_symboles': L, 'erreurs_symboles': e, 'symboles': n}
    if e: ECHECS.append(f'{k} : {e} erreurs sans bruit')
    ech = d['constellation_rx'][0]
    if k != 'fsk2':
        ideal = np.array(m['pts']) / m['amp_moy']
        n2 = min(len(ech) - L, len(tx)) - 20
        ref = ideal[tx[20:n2].astype(int)]; mes = ech[20 + L:n2 + L]
        gain = float(np.real(np.vdot(ref, mes) / np.vdot(ref, ref)))
        evm = float(np.sqrt(np.mean(np.abs(mes - ref) ** 2) / np.mean(np.abs(ref) ** 2)))
        r['sans_bruit'].update({'gain_mesure': round(gain, 4), 'evm': round(evm, 4)})
        if evm > 0.05 or abs(gain - 1) > 0.05: ECHECS.append(f'{k} : EVM {evm:.3f}, gain {gain:.3f}')
        # meilleure phase d'échantillonnage mesurée sur la sortie du filtre adapté (oeil)
        oe = d['oeil'][0]; ev = []
        for p in range(SPS):
            s = oe[p::SPS]; nn = min(len(s), len(tx) + 40) - 40
            ev.append(min(np.mean(np.abs(s[l + 20:l + nn] - ideal[tx[20:nn].astype(int)]) ** 2) for l in range(30)))
        r['sans_bruit']['phase_optimale'] = int(np.argmin(ev)); r['sans_bruit']['evm_par_phase_min'] = float(np.sqrt(min(ev)))
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
    r['taux'] = []
    for e_db in m['ebn0']:
        g = 10 ** (e_db / 10); sig = float(np.sqrt(Ps * SPS / (2 * m['bits'] * g)))
        d2, info = executer(nom, 20000, sig, f'eb{e_db}')
        tx2, dec2 = d2['symboles'][0], d2['symboles'][1]; L2 = aligner(tx2, dec2)
        es, ns, eb, nb = erreurs(m, tx2, dec2, L2)
        mesure = eb / nb if m['compare'] == 'TEB' else es / ns
        nerr = eb if m['compare'] == 'TEB' else es
        th = theorie(m, g); ok = nerr < 50 or 1 / 1.5 < mesure / th < 1.5
        r['taux'].append({'EbN0_dB': e_db, 'ecart_type_bruit': round(sig, 5), 'TES': es / ns, 'TEB': eb / nb,
                          'erreurs': nerr, 'mesure': mesure, 'theorie': th, 'conforme': bool(ok), 'retard': L2, **info})
        if not ok: ECHECS.append(f'{k} : Eb/N0 = {e_db} dB, mesure {mesure:.3g}, théorie {th:.3g}')
        if e_db == m['ebn0'][2]: d_mid = d2
    d_eye, _ = executer(nom, 3000, float(np.sqrt(Ps * SPS / (2 * m['bits'] * 10 ** 2.0))), 'oeil20')
    # figure : signal émis, spectre, constellation, oeil
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    t = np.arange(8 * SPS) / FS * 1e3; i0 = SPS * 40
    ax[0, 0].plot(t, emis[i0:i0 + 8 * SPS], lw=0.9, label='signal émis')
    ax[0, 0].plot(t, d['scope_emis'][1][i0:i0 + 8 * SPS], lw=1.6, ls='--',
                  label='bits NRZ' if k == 'fsk2' else 'voie I en bande de base')
    ax[0, 0].set_xlabel('temps (ms)'); ax[0, 0].set_title(f"{m['titre']} : signal émis (scope_emis)"); ax[0, 0].legend(loc='upper right', fontsize=8)
    ax[0, 1].plot(f / 1e3, 10 * np.log10(Pn + 1e-16)); ax[0, 1].set_xlim(0, 10); ax[0, 1].set_ylim(-80, 3)
    ax[0, 1].set_xlabel('fréquence (kHz)'); ax[0, 1].set_ylabel('dB'); ax[0, 1].set_title('spectre du signal émis (spectre_emis)')
    c = d_mid['constellation_rx'][0][50:]
    ax[1, 0].plot(c.real, c.imag, '.', ms=1.5, alpha=0.4, label=f"Eb/N0 = {m['ebn0'][2]} dB")
    c0 = d['constellation_rx'][0][50:]; ax[1, 0].plot(c0.real, c0.imag, 'r+', ms=9, label='sans bruit')
    ax[1, 0].set_aspect('equal'); ax[1, 0].legend(loc='upper right', fontsize=8)
    ax[1, 0].set_title('points de décision (|y1|, |y0|) (constellation_rx)' if k == 'fsk2' else 'constellation (constellation_rx)')
    oe = d_eye['oeil'][0]; oe = oe.real if np.iscomplexobj(oe) else oe
    ph = r['sans_bruit'].get('phase_optimale', SPS - 1)
    for j in range(60, 260):
        s0 = j * SPS + ph - SPS // 2
        ax[1, 1].plot(np.arange(2 * SPS + 1) / SPS - 0.5, oe[s0:s0 + 2 * SPS + 1], color='C0', lw=0.4, alpha=0.5)
    ax[1, 1].set_xlabel('temps (symboles)'); ax[1, 1].set_title("diagramme de l'oeil, Eb/N0 = 20 dB (oeil)")
    for a in ax.flat: a.grid(ls=':')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f'{nom}.png'), dpi=110); plt.close(fig)
    print(json.dumps(r, indent=1, default=float), flush=True)

fig, ax = plt.subplots(figsize=(9, 6))
for i, (k, m) in enumerate(MODS.items()):
    e = np.array([p['EbN0_dB'] for p in R[k]['taux']]); y = np.array([p['mesure'] for p in R[k]['taux']])
    ee = np.linspace(e.min() - 1, e.max() + 1, 60)
    ax.semilogy(ee, [theorie(m, 10 ** (x / 10)) for x in ee], color=f'C{i}', lw=1)
    ax.semilogy(e, np.where(y > 0, y, np.nan), 'o', color=f'C{i}', label=f"{m['titre']} ({m['compare']})")
ax.set_ylim(1e-5, 0.5); ax.set_xlabel('Eb/N0 (dB)'); ax.set_ylabel("taux d'erreur")
ax.set_title("Taux d'erreur mesurés dans GNU Radio (points) et théorie (traits)"); ax.grid(which='both', ls=':'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'taux_erreur_gnuradio.png'), dpi=110); plt.close(fig)
R['echecs'] = ECHECS
json.dump(R, open(os.path.join(OUT, 'mesures.json'), 'w'), indent=1, default=float)
print('ECHECS :', ECHECS if ECHECS else 'aucun')
sys.exit(1 if ECHECS else 0)
