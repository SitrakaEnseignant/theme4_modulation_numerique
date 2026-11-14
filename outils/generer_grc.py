"""Génère les graphes de flux GNU Radio Companion 3.10 du Thème 4 (format YAML), en simulation pure :
aucune source ni aucun puits matériel. Les paramètres de chaque bloc reprennent ceux de fichiers .grc
de GNU Radio 3.10 (outils/modeles_blocs.json) ; seules les valeurs utiles sont modifiées.
    python3 outils/generer_grc.py"""
import os, sys, json, copy, yaml
import numpy as np
ICI = os.path.dirname(os.path.abspath(__file__)); RACINE = os.path.dirname(ICI)
sys.path.insert(0, os.path.join(RACINE, 'En_ligne'))
from modulations import MODS, FS, SPS, FC, ALPHA, NTAPS, DF, PHI_OL, DF_OL, W_COSTAS, W_RECEPTEUR, table
MODELES = json.load(open(os.path.join(ICI, 'modeles_blocs.json'), encoding='utf8'))
OUT = os.path.join(RACINE, 'GNU_Radio')


def rrc_gr(gain, fs, rs, alpha, ntaps):
    """réplique de firdes.root_raised_cosine de GNU Radio 3.10 (somme des coefficients = gain)"""
    ntaps |= 1; spb = fs / rs; h = np.zeros(ntaps); scale = 0.0
    for i in range(ntaps):
        x = i - ntaps // 2; x1 = np.pi * x / spb; x2 = 4 * alpha * x / spb; x3 = x2 * x2 - 1
        if abs(x3) >= 1e-6:
            num = np.cos((1 + alpha) * x1) + (np.sin((1 - alpha) * x1) / (4 * alpha * x / spb) if i != ntaps // 2
                                               else (1 - alpha) * np.pi / (4 * alpha))
            den = x3 * np.pi
        else:
            x3b = (1 - alpha) * x1; x2b = (1 + alpha) * x1
            num = (np.sin(x2b) * (1 + alpha) * np.pi - np.cos(x3b) * ((1 - alpha) * np.pi * spb) / (4 * alpha * x)
                   + np.sin(x3b) * spb * spb / (4 * alpha * x * x))
            den = -32 * np.pi * alpha * alpha * x / spb
        h[i] = 4 * alpha * num / den; scale += h[i]
    return h * gain / scale


class Q(str): pass
class L(str): pass


class FG:
    def __init__(s, fid, titre, desc):
        s.fid, s.titre, s.desc, s.blocks, s.conns = fid, titre, desc, [], []

    def blk(s, nom, bid, x, y, etat='enabled', **p):
        base = copy.deepcopy(MODELES.get(bid, {}))
        for k in [k for k in base if k.startswith(('label', 'color', 'alpha', 'width', 'style', 'marker'))
                  and k[-1].isdigit() and int(''.join(c for c in k if c.isdigit()) or 0) > 10]:
            base.pop(k)
        base.update({k: str(v) for k, v in p.items()})
        if 'comment' not in base: base['comment'] = ''
        base = {k: Q(v) for k, v in base.items()}
        s.blocks.append({'name': nom, 'id': bid, 'parameters': base,
                         'states': {'bus_sink': False, 'bus_source': False, 'bus_structure': None,
                                    'coordinate': [int(x), int(y)], 'rotation': 0, 'state': etat}})
        return nom

    def c(s, a, b, pa=0, pb=0):
        s.conns.append([a, str(pa), b, str(pb)])

    def save(s):
        d = {'options': {'parameters': {'author': 'Thème 4 - modulations numériques (simulation pure)', 'catch_exceptions': 'True',
                 'category': '[GRC Hier Blocks]', 'cmake_opt': '', 'comment': '', 'copyright': '', 'description': s.desc,
                 'gen_cmake': 'On', 'gen_linking': 'dynamic', 'generate_options': 'qt_gui', 'hier_block_src_path': '.:',
                 'id': s.fid, 'max_nouts': '0', 'output_language': 'python', 'placement': '(0,0)', 'qt_qss_theme': '',
                 'realtime_scheduling': '', 'run': 'True', 'run_command': '{python} -u {filename}', 'run_options': 'prompt',
                 'sizing_mode': 'fixed', 'thread_safe_setters': '', 'title': s.titre, 'window_size': '(1200,1000)'},
                 'states': {'bus_sink': False, 'bus_source': False, 'bus_structure': None, 'coordinate': [8, 8],
                            'rotation': 0, 'state': 'enabled'}},
             'blocks': s.blocks, 'connections': s.conns,
             'metadata': {'file_format': 1, 'grc_version': '3.10.9.2'}}

        class D(yaml.SafeDumper): pass
        D.add_representer(Q, lambda dp, v: dp.represent_scalar('tag:yaml.org,2002:str', str(v), style="'"))
        with open(os.path.join(OUT, s.fid + '.grc'), 'w', encoding='utf8') as f:
            yaml.dump(d, f, Dumper=D, sort_keys=False, allow_unicode=True, default_flow_style=False, width=1000)


def rng(g, nom, val, a, b, pas, lab, hint, x, y=110):
    return g.blk(nom, 'qtgui_range', x, y, value=val, start=a, stop=b, step=pas, label=lab, gui_hint=hint,
                 widget='counter_slider', rangeType='float', orient='QtCore.Qt.Horizontal', min_len='200', outputmsgname='value')


def time_sink(g, nom, typ, n, labels, hint, x, y, size='512', ymin='-2', ymax='2'):
    p = dict(type=typ, name=f'"{nom}"', size=size, srate='samp_rate', nconnections=n, ymin=ymin, ymax=ymax, gui_hint=hint,
             autoscale='False', grid='True', legend='True', update_time='0.10', entags='False', ctrlpanel='False', stemplot='False')
    for i, l in enumerate(labels, 1): p[f'label{i}'] = f'"{l}"'
    return g.blk(nom, 'qtgui_time_sink_x', x, y, **p)


def freq_sink(g, nom, typ, hint, x, y):
    return g.blk(nom, 'qtgui_freq_sink_x', x, y, type=typ, name=f'"{nom}"', fftsize='4096', bw='samp_rate', fc='0',
                 nconnections='1', average='0.1', ymin='-100', ymax='10', gui_hint=hint, autoscale='False', grid='True',
                 legend='True', freqhalf='True', label1='"signal émis"', wintype='window.WIN_BLACKMAN_hARRIS')


def const_sink(g, nom, hint, x, y, lim='2.5'):
    return g.blk(nom, 'qtgui_const_sink_x', x, y, type='complex', name=f'"{nom}"', size='1024', nconnections='1',
                 xmin=f'-{lim}', xmax=lim, ymin=f'-{lim}', ymax=lim, gui_hint=hint, autoscale='False', grid='True', label1='"échantillons"')


def eye_sink(g, nom, typ, hint, x, y, lim='3'):
    return g.blk(nom, 'qtgui_eye_sink_x', x, y, type=typ, samp_per_symbol='sps', size='1024', srate='samp_rate',
                 nconnections='1', ymin=f'-{lim}', ymax=lim, gui_hint=hint, autoscale='False', grid='True', entags='False',
                 label1='"diagramme de l\'oeil"', ylabel='Amplitude', ctrlpanel='False')


def variables(g, m, noise='0.0'):
    g.blk('imports', 'import', 1300, 12, imports='import math\nimport numpy as np')
    g.blk('samp_rate', 'variable', 16, 12, value=repr(FS)); g.blk('sps', 'variable', 150, 12, value=str(SPS))
    g.blk('nsym', 'variable', 260, 12, value='20000')
    rng(g, 'fc', repr(FC), '2000', '8000', '100', 'Fréquence de la porteuse fc (Hz)', '0,0,1,1', 16)
    rng(g, 'noise_amp', noise, '0', '1', '0.005', 'Écart type du bruit du canal', '0,1,1,1', 230)
    rng(g, 'phi_ol', repr(PHI_OL), '-3.14', '3.14', '0.01', "Erreur de phase de l'OL du récepteur (rad)", '0,2,1,1', 440)
    rng(g, 'df_ol', repr(DF_OL), '-50', '50', '0.5', "Erreur de fréquence de l'OL du récepteur (Hz)", '0,3,1,1', 650)


def recepteur_ol(g, x, y, freq='-(fc+df_ol)', nom='lo_rx'):
    return g.blk(nom, 'analog_sig_source_x', x, y, type='complex', samp_rate='samp_rate', waveform='analog.GR_COS_WAVE',
                 freq=freq, amp='2', offset='0', phase='-phi_ol', showports='False')


def lineaire(m):
    nom, M, rrc, det = m['nom'], m['M'], m['mise_en_forme'] == 'rrc', m['detection']
    g = FG(f'{nom}_tx_rx_sim', f"{m['titre']} : émetteur, canal et récepteur (simulation)",
           f"{m['titre']}, {M} symboles, {'racine de cosinus surélevé' if rrc else 'impulsions rectangulaires'}, porteuse fc, "
           "bruit blanc gaussien, oscillateur local du récepteur décalé en phase et en fréquence, filtre adapté, "
           + {'enveloppe': "détection d'enveloppe.", 'costas': 'boucle de Costas, échantillonnage et décision.',
              'recepteur': 'échantillonnage et récepteur de constellation (boucle de phase pilotée par les décisions).'}[det])
    variables(g, m)
    g.blk('constellation', 'variable_constellation', 860, 12, const_points=table(m['pts']), sym_map=str(list(range(M))),
          normalization='digital.constellation.AMPLITUDE_NORMALIZATION', rot_sym=str(m['ordre']), dims='1', precision='8',
          npwr='1.0', soft_dec_lut='None', type='calcdist')
    if rrc:
        g.blk('alpha', 'variable', 1060, 12, value=repr(ALPHA)); g.blk('ntaps', 'variable', 1160, 12, value=str(NTAPS))
        h = rrc_gr(SPS, FS, FS / SPS, ALPHA, NTAPS); hr = rrc_gr(1.0, FS, FS / SPS, ALPHA, NTAPS)
        gain_mf = float(np.convolve(h, hr).max()); decal = 1
    else:
        gain_mf, decal = 1.0, 0
    g.blk('gain_ech', 'variable', 360, 12, value=repr(round(1.0 / (gain_mf * m['amp_moy']), 6)))
    g.blk('decalage', 'variable', 500, 12, value=str(decal))
    if det != 'enveloppe':
        g.blk('w_boucle', 'variable', 620, 12, value=W_COSTAS if det == 'costas' else W_RECEPTEUR)
    y1, y2, y3 = 220, 400, 580
    g.blk('src', 'analog_random_source_x', 16, y1, type='byte', min='0', max=str(M), num_samps='nsym', repeat='True')
    g.blk('symboles_tx', 'digital_chunks_to_symbols_xx', 200, y1, in_type='byte', out_type='complex',
          symbol_table=table(m['pts']), dimension='1', num_ports='1')
    if rrc:
        g.blk('mise_en_forme', 'root_raised_cosine_filter', 400, y1, type='interp_fir_filter_ccf', interp='sps', decim='1',
              gain='sps', samp_rate='samp_rate', sym_rate='samp_rate/sps', alpha='alpha', ntaps='ntaps', comment='')
    else:
        g.blk('mise_en_forme', 'blocks_repeat', 400, y1, type='complex', interp='sps', vlen='1')
    g.blk('throttle', 'blocks_throttle', 600, y1, type='complex', samples_per_second='samp_rate', ignoretag='True', vlen='1')
    g.blk('lo_tx', 'analog_sig_source_x', 600, y1 + 90, type='complex', samp_rate='samp_rate', waveform='analog.GR_COS_WAVE',
          freq='fc', amp='1', offset='0', phase='0', showports='False')
    g.blk('mix_tx', 'blocks_multiply_xx', 800, y1, type='complex', num_inputs='2', vlen='1')
    g.blk('emis_reel', 'blocks_complex_to_real', 960, y1, vlen='1')
    g.blk('bande_base_I', 'blocks_complex_to_real', 800, y1 + 110, vlen='1')
    g.blk('bruit', 'analog_noise_source_x', 960, y1 + 90, type='float', noise_type='analog.GR_GAUSSIAN', amp='noise_amp', seed='1')
    g.blk('canal', 'blocks_add_xx', 1150, y1, type='float', num_inputs='2', vlen='1')
    time_sink(g, 'scope_emis', 'float', 2, ['signal émis', 'voie I en bande de base'], '1,0,1,2', 1150, y1 - 110, size='384')
    freq_sink(g, 'spectre_emis', 'float', '1,2,1,2', 1350, y1 - 110)
    g.blk('zero', 'analog_const_source_x', 16, y2 + 80, type='float', const='0')
    g.blk('to_cplx', 'blocks_float_to_complex', 200, y2, vlen='1')
    recepteur_ol(g, 200, y2 + 100)
    g.blk('mix_rx', 'blocks_multiply_xx', 400, y2, type='complex', num_inputs='2', vlen='1')
    if rrc:
        g.blk('filtre_adapte', 'root_raised_cosine_filter', 560, y2, type='fir_filter_ccf', interp='1', decim='1', gain='1',
              samp_rate='samp_rate', sym_rate='samp_rate/sps', alpha='alpha', ntaps='ntaps', comment='')
    else:
        g.blk('filtre_adapte', 'blocks_moving_average_xx', 560, y2, type='complex', length='sps', scale='1.0/sps',
              max_iter='4000', vlen='1')
    g.blk('normalisation', 'blocks_multiply_const_vxx', 760, y2, type='complex', const='gain_ech', vlen='1')
    g.blk('tx_float', 'blocks_char_to_float', 200, y1 - 110, scale='1', vlen='1')
    C = [('src', 'symboles_tx'), ('symboles_tx', 'mise_en_forme'), ('mise_en_forme', 'throttle'),
         ('lo_tx', 'mix_tx', 0, 1), ('throttle', 'mix_tx', 0, 0), ('mix_tx', 'emis_reel'), ('throttle', 'bande_base_I'),
         ('emis_reel', 'canal', 0, 0), ('bruit', 'canal', 0, 1), ('emis_reel', 'scope_emis', 0, 0),
         ('bande_base_I', 'scope_emis', 0, 1), ('emis_reel', 'spectre_emis'), ('canal', 'to_cplx', 0, 0),
         ('zero', 'to_cplx', 0, 1), ('to_cplx', 'mix_rx', 0, 0), ('lo_rx', 'mix_rx', 0, 1), ('mix_rx', 'filtre_adapte'),
         ('filtre_adapte', 'normalisation'), ('src', 'tx_float'), ('tx_float', 'symboles', 0, 0)]
    if det == 'enveloppe':
        g.blk('module', 'blocks_complex_to_mag', 960, y2, vlen='1')
        eye_sink(g, 'oeil', 'float', '2,0,1,2', 1150, y2 + 110)
        g.blk('decalage_ech', 'blocks_skiphead', 1150, y2, type='float', num_items='decalage', vlen='1')
        g.blk('echantillonnage', 'blocks_keep_one_in_n', 1350, y2, type='float', n='sps', vlen='1')
        g.blk('decalage_c', 'blocks_skiphead', 960, y3, type='complex', num_items='decalage', vlen='1')
        g.blk('echantillonnage_c', 'blocks_keep_one_in_n', 1150, y3, type='complex', n='sps', vlen='1')
        const_sink(g, 'constellation_rx', '2,2,1,2', 1350, y3)
        C += [('normalisation', 'module'), ('module', 'oeil'), ('module', 'decalage_ech'), ('decalage_ech', 'echantillonnage'),
              ('normalisation', 'decalage_c'), ('decalage_c', 'echantillonnage_c'), ('echantillonnage_c', 'constellation_rx')]
        if nom == 'ook':
            g.blk('decision', 'blocks_threshold_ff', 1550, y2, low='1.0', high='1.0', init='0'); dec_out = 'decision'
            C += [('echantillonnage', 'decision')]
        else:
            g.blk('zero_m', 'analog_const_source_x', 1350, y2 + 90, type='float', const='0')
            g.blk('module_c', 'blocks_float_to_complex', 1550, y2, vlen='1')
            g.blk('decision', 'digital_constellation_decoder_cb', 1700, y2, constellation='constellation')
            g.blk('dec_float', 'blocks_char_to_float', 1850, y2, scale='1', vlen='1'); dec_out = 'dec_float'
            C += [('echantillonnage', 'module_c', 0, 0), ('zero_m', 'module_c', 0, 1), ('module_c', 'decision'), ('decision', 'dec_float')]
    elif det == 'costas':
        g.blk('costas', 'digital_costas_loop_cc', 960, y2, w='w_boucle', order=str(m['ordre']), use_snr='False')
        eye_sink(g, 'oeil', 'complex', '2,0,1,2', 1150, y2 + 110)
        g.blk('decalage_ech', 'blocks_skiphead', 1150, y2, type='complex', num_items='decalage', vlen='1')
        g.blk('echantillonnage', 'blocks_keep_one_in_n', 1350, y2, type='complex', n='sps', vlen='1')
        const_sink(g, 'constellation_rx', '2,2,1,2', 1550, y2 + 110)
        g.blk('decision', 'digital_constellation_decoder_cb', 1550, y2, constellation='constellation')
        g.blk('dec_float', 'blocks_char_to_float', 1700, y2, scale='1', vlen='1'); dec_out = 'dec_float'
        C += [('normalisation', 'costas'), ('costas', 'oeil'), ('costas', 'decalage_ech'), ('decalage_ech', 'echantillonnage'),
              ('echantillonnage', 'constellation_rx'), ('echantillonnage', 'decision'), ('decision', 'dec_float')]
    else:
        eye_sink(g, 'oeil', 'complex', '2,0,1,2', 960, y2 + 110)
        g.blk('decalage_ech', 'blocks_skiphead', 960, y2, type='complex', num_items='decalage', vlen='1')
        g.blk('echantillonnage', 'blocks_keep_one_in_n', 1150, y2, type='complex', n='sps', vlen='1')
        g.blk('recepteur', 'digital_constellation_receiver_cb', 1350, y2, constellation='constellation', loop_bw='w_boucle',
              fmin='-0.1', fmax='0.1')
        g.blk('dec_float', 'blocks_char_to_float', 1600, y2, scale='1', vlen='1'); dec_out = 'dec_float'
        const_sink(g, 'constellation_rx', '2,2,1,2', 1600, y2 + 110)
        time_sink(g, 'boucle', 'float', 2, ['phase estimée (rad)', 'fréquence estimée (rad/symbole)'], '3,2,1,2', 1600, y3,
                  size='1024', ymin='-4', ymax='4')
        g.blk('erreur_phase', 'blocks_null_sink', 1600, y3 + 120, type='float', num_inputs='1', vlen='1', bus_structure_sink='[[0,],]')
        C += [('normalisation', 'oeil'), ('normalisation', 'decalage_ech'), ('decalage_ech', 'echantillonnage'),
              ('echantillonnage', 'recepteur'), ('recepteur', 'dec_float', 0, 0), ('recepteur', 'erreur_phase', 1, 0),
              ('recepteur', 'boucle', 2, 0), ('recepteur', 'boucle', 3, 1), ('recepteur', 'constellation_rx', 4, 0)]
    time_sink(g, 'symboles', 'float', 2, ['symbole émis', 'symbole décidé'], '3,0,1,2', 1850, y3, size='64',
              ymin='-1', ymax=str(M))
    C += [(dec_out, 'symboles', 0, 1)]
    for c in C: g.c(*c)
    g.save()
    return g.fid


def fsk(m):
    g = FG('fsk2_tx_rx_sim', 'FSK binaire : modulateur de fréquence, canal, détection non cohérente et discriminateur (simulation)',
           'FSK binaire à phase continue, tons fc - 1000 Hz et fc + 1000 Hz (indice 2), bruit blanc gaussien, oscillateur local '
           'décalé, détection non cohérente par deux corrélateurs, et discriminateur de fréquence suivi d\'un comparateur.')
    variables(g, m)
    g.blk('excursion', 'variable', 860, 12, value=repr(DF))
    y1, y2, y3, y4 = 220, 400, 580, 760
    g.blk('src', 'analog_random_source_x', 16, y1, type='byte', min='0', max='2', num_samps='nsym', repeat='True')
    g.blk('nrz', 'digital_chunks_to_symbols_xx', 200, y1, in_type='byte', out_type='float', symbol_table='[-1.0, 1.0]',
          dimension='1', num_ports='1')
    g.blk('mise_en_forme', 'blocks_repeat', 400, y1, type='float', interp='sps', vlen='1')
    g.blk('fsk_mod', 'analog_frequency_modulator_fc', 560, y1, sensitivity='2*math.pi*excursion/samp_rate')
    g.blk('throttle', 'blocks_throttle', 760, y1, type='complex', samples_per_second='samp_rate', ignoretag='True', vlen='1')
    g.blk('lo_tx', 'analog_sig_source_x', 760, y1 + 90, type='complex', samp_rate='samp_rate', waveform='analog.GR_COS_WAVE',
          freq='fc', amp='1', offset='0', phase='0', showports='False')
    g.blk('mix_tx', 'blocks_multiply_xx', 960, y1, type='complex', num_inputs='2', vlen='1')
    g.blk('emis_reel', 'blocks_complex_to_real', 1100, y1, vlen='1')
    g.blk('bruit', 'analog_noise_source_x', 1100, y1 + 90, type='float', noise_type='analog.GR_GAUSSIAN', amp='noise_amp', seed='1')
    g.blk('canal', 'blocks_add_xx', 1280, y1, type='float', num_inputs='2', vlen='1')
    time_sink(g, 'scope_emis', 'float', 2, ['signal émis', 'bits NRZ'], '1,0,1,2', 1280, y1 - 110, size='384')
    freq_sink(g, 'spectre_emis', 'float', '1,2,1,2', 1480, y1 - 110)
    g.blk('zero', 'analog_const_source_x', 16, y2 + 80, type='float', const='0')
    g.blk('to_cplx', 'blocks_float_to_complex', 200, y2, vlen='1')
    C = []
    for sg, yy, t in (('+', y2, '1'), ('-', y3, '0')):
        recepteur_ol(g, 200, yy + 90, freq=f'-(fc{sg}excursion+df_ol)', nom=f'lo_{t}')
        g.blk(f'mix_{t}', 'blocks_multiply_xx', 400, yy, type='complex', num_inputs='2', vlen='1')
        g.blk(f'correlateur_{t}', 'blocks_moving_average_xx', 560, yy, type='complex', length='sps', scale='1.0/sps',
              max_iter='4000', vlen='1')
        g.blk(f'module_{t}', 'blocks_complex_to_mag', 760, yy, vlen='1')
        C += [('to_cplx', f'mix_{t}', 0, 0), (f'lo_{t}', f'mix_{t}', 0, 1), (f'mix_{t}', f'correlateur_{t}'),
              (f'correlateur_{t}', f'module_{t}')]
    g.blk('difference', 'blocks_sub_xx', 960, y2, type='float', num_inputs='2', vlen='1')
    eye_sink(g, 'oeil', 'float', '2,0,1,2', 1100, y2 + 90, lim='1.5')
    g.blk('decalage_ech', 'blocks_skiphead', 1100, y2, type='float', num_items='0', vlen='1')
    g.blk('echantillonnage', 'blocks_keep_one_in_n', 1280, y2, type='float', n='sps', vlen='1')
    g.blk('decision', 'digital_binary_slicer_fb', 1450, y2)
    g.blk('dec_float', 'blocks_char_to_float', 1600, y2, scale='1', vlen='1')
    g.blk('paire', 'blocks_float_to_complex', 960, y3, vlen='1')
    g.blk('decalage_paire', 'blocks_skiphead', 1100, y3, type='complex', num_items='0', vlen='1')
    g.blk('echantillonnage_paire', 'blocks_keep_one_in_n', 1280, y3, type='complex', n='sps', vlen='1')
    const_sink(g, 'constellation_rx', '2,2,1,2', 1450, y3, lim='1.5')
    recepteur_ol(g, 200, y4 + 90)
    g.blk('mix_rx', 'blocks_multiply_xx', 400, y4, type='complex', num_inputs='2', vlen='1')
    g.blk('passe_bas', 'low_pass_filter', 560, y4, type='fir_filter_ccf', decim='1', interp='1', gain='1',
          samp_rate='samp_rate', cutoff_freq='3000', width='1000', win='window.WIN_HAMMING', beta='6.76')
    g.blk('discriminateur', 'analog_quadrature_demod_cf', 760, y4, gain='samp_rate/(2*math.pi*excursion)')
    g.blk('integrateur', 'blocks_moving_average_xx', 960, y4, type='float', length='sps', scale='1.0/sps', max_iter='4000', vlen='1')
    g.blk('decalage_d', 'blocks_skiphead', 1100, y4, type='float', num_items='0', vlen='1')
    g.blk('echantillonnage_d', 'blocks_keep_one_in_n', 1280, y4, type='float', n='sps', vlen='1')
    g.blk('decision_d', 'digital_binary_slicer_fb', 1450, y4)
    g.blk('dec_d_float', 'blocks_char_to_float', 1600, y4, scale='1', vlen='1')
    time_sink(g, 'scope_discriminateur', 'float', 2, ['bits NRZ émis', 'sortie du discriminateur'], '3,2,1,2', 960, y4 + 110, size='384')
    g.blk('tx_float', 'blocks_char_to_float', 200, y1 - 110, scale='1', vlen='1')
    time_sink(g, 'symboles', 'float', 3, ['bit émis', 'bit décidé (corrélateurs)', 'bit décidé (discriminateur)'], '3,0,1,2',
              1800, y3, size='64', ymin='-1', ymax='2')
    C += [('src', 'nrz'), ('nrz', 'mise_en_forme'), ('mise_en_forme', 'fsk_mod'), ('fsk_mod', 'throttle'),
          ('throttle', 'mix_tx', 0, 0), ('lo_tx', 'mix_tx', 0, 1), ('mix_tx', 'emis_reel'), ('emis_reel', 'canal', 0, 0),
          ('bruit', 'canal', 0, 1), ('emis_reel', 'scope_emis', 0, 0), ('mise_en_forme', 'scope_emis', 0, 1),
          ('emis_reel', 'spectre_emis'), ('canal', 'to_cplx', 0, 0), ('zero', 'to_cplx', 0, 1),
          ('module_1', 'difference', 0, 0), ('module_0', 'difference', 0, 1), ('difference', 'oeil'),
          ('difference', 'decalage_ech'), ('decalage_ech', 'echantillonnage'), ('echantillonnage', 'decision'),
          ('decision', 'dec_float'), ('module_1', 'paire', 0, 0), ('module_0', 'paire', 0, 1), ('paire', 'decalage_paire'),
          ('decalage_paire', 'echantillonnage_paire'), ('echantillonnage_paire', 'constellation_rx'),
          ('to_cplx', 'mix_rx', 0, 0), ('lo_rx', 'mix_rx', 0, 1), ('mix_rx', 'passe_bas'), ('passe_bas', 'discriminateur'),
          ('discriminateur', 'integrateur'), ('integrateur', 'decalage_d'), ('decalage_d', 'echantillonnage_d'),
          ('echantillonnage_d', 'decision_d'), ('decision_d', 'dec_d_float'),
          ('mise_en_forme', 'scope_discriminateur', 0, 0), ('discriminateur', 'scope_discriminateur', 0, 1),
          ('src', 'tx_float'), ('tx_float', 'symboles', 0, 0), ('dec_float', 'symboles', 0, 1), ('dec_d_float', 'symboles', 0, 2)]
    for c in C: g.c(*c)
    g.save()
    return g.fid


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for k, m in MODS.items():
        print(fsk(m) if k == 'fsk2' else lineaire(m))
