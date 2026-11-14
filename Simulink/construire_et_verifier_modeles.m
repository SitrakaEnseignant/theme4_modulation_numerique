%% CONSTRUIRE_ET_VERIFIER_MODELES - Bloc 2, Thème 4 : modulations numériques
% Construit, avec des blocs de la bibliothèque Simulink, les modèles mod_ook, mod_ask4, mod_bpsk,
% mod_qpsk, mod_psk8, mod_qam16 et mod_fsk2 (simulation pure, sans matériel). Chaque modèle reprend
% bloc pour bloc le graphe de flux GNU Radio de même nom (GNU_Radio/*_tx_rx_sim.grc) : mêmes noms de
% blocs, même ordre, mêmes paramètres, y compris l'oscillateur local du récepteur décalé en phase et
% en fréquence, la boucle de Costas (BPSK, QPSK, 8-PSK) et le récepteur de constellation (16-QAM).
% Le script exécute chaque modèle sans bruit puis pour plusieurs Eb/N0 et vérifie : absence d'erreur
% sans bruit, constellation, spectre, diagramme de l'oeil et taux d'erreur comparés à la théorie.
% Les mêmes paramètres figurent dans En_ligne/modulations.py.
clear; clc;
ici = fileparts(mfilename('fullpath')); cd(ici);
if exist(fullfile(ici, 'journal_matlab.txt'), 'file'), delete(fullfile(ici, 'journal_matlab.txt')); end
diary(fullfile(ici, 'journal_matlab.txt'));
fprintf('MATLAB %s\n', version);
P = struct('fs', 32000, 'sps', 32, 'Rs', 1000, 'fc', 4000, 'alpha', 0.35, 'ntaps', 353, 'df', 1000, ...
    'phi_ol', 0.7, 'df_ol', 1.0, 'w_costas', 2 * pi / 12800, 'w_recepteur', 2 * pi / 400, 'conv', 300);
mods = definir_modulations();
echecs = {}; R = struct();
for i = 1:numel(mods)
    m = mods(i);
    fprintf('\n===== mod_%s : %s =====\n', m.nom, m.titre);
    setappdata(0, 'liste_blocs', {});
    if strcmp(m.nom, 'fsk2'), mdl = construire_fsk(m, P); else, mdl = construire_lineaire(m, P); end
    fid = fopen(fullfile(ici, [mdl '_blocs.json']), 'w'); fprintf(fid, '%s', jsonencode(getappdata(0, 'liste_blocs'), 'PrettyPrint', true)); fclose(fid);
    save_system(mdl, fullfile(ici, [mdl '.slx']));
    try, print(['-s' mdl], '-dpng', '-r150', fullfile(ici, [mdl '_schema.png'])); catch err, fprintf('Schema non exporte : %s\n', err.message); end
    [r, e] = verifier(mdl, m, P, ici);
    R.(m.nom) = r; echecs = [echecs, e]; %#ok<AGROW>
    close_system(mdl, 0);
end
figure_taux(mods, R, ici);
fid = fopen(fullfile(ici, 'mesures_simulink.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
fprintf('\nControles en echec : %d\n', numel(echecs));
for k = 1:numel(echecs), fprintf('  ECHEC %s\n', echecs{k}); end
diary off;
if ~isempty(echecs), error('Des controles ont echoue (voir journal_matlab.txt).'); end

%% ------------------------------------------------------------------ modulations
function mods = definir_modulations()
g = @(k) bitxor(k, bitshift(k, -1));
psk = @(M, rot) exp(1j * (rot + 2 * pi * (0:M - 1) / M));
niv = [-1, -1/3, 1/3, 1]; q16 = zeros(1, 16); l16 = zeros(1, 16); n = 0;
for a = 0:3, for b = 0:3, n = n + 1; q16(n) = niv(a + 1) + 1j * niv(b + 1); l16(n) = bitor(bitshift(g(a), 2), g(b)); end, end
D = {
 'ook',   'OOK (tout ou rien)', [0 1],          [0 1],  'rect', 'enveloppe', 1, [6 8 10 12 14],    'TEB'
 'ask4',  'ASK a 4 niveaux',    [0 1/3 2/3 1],  g(0:3), 'rect', 'enveloppe', 1, [10 12 14 16 18],  'TES'
 'bpsk',  'BPSK',               [1 -1],         [0 1],  'rect', 'costas',    2, [0 2 4 6 8],       'TES'
 'qpsk',  'QPSK',               psk(4, pi/4),   g(0:3), 'rect', 'costas',    4, [0 2 4 6 8],       'TES'
 'psk8',  '8-PSK',              psk(8, pi/8),   g(0:7), 'rrc',  'costas',    8, [4 6 8 10 12],     'TES'
 'qam16', '16-QAM',             q16,            l16,    'rrc',  'recepteur', 4, [4 6 8 10 12],     'TES'
 'fsk2',  'FSK binaire',        [-1 1],         [0 1],  'rect', 'energie',   1, [4 6 8 10 12],     'TEB'};
for k = 1:size(D, 1)
    p = D{k, 3}; am = mean(abs(p)); if strcmp(D{k, 1}, 'fsk2'), am = 1; end
    mods(k) = struct('nom', D{k, 1}, 'titre', D{k, 2}, 'pts', p, 'lab', D{k, 4}, 'forme', D{k, 5}, 'detection', D{k, 6}, ...
        'ordre', D{k, 7}, 'ebn0', D{k, 8}, 'compare', D{k, 9}, 'M', numel(p), 'bits', log2(numel(p)), 'amp_moy', am, ...
        'ideal', p / am); %#ok<AGROW>
end
end

%% ------------------------------------------------------------------ construction
function b = ajouter(mdl, lib, nom, col, lig, varargin)
x = 30 + 130 * col; y = 30 + 90 * lig;
b = [mdl '/' nom];
L = getappdata(0, 'liste_blocs'); L{end + 1} = struct('bloc', nom, 'bibliotheque', strrep(lib, 'simulink/', ''), 'parametres', {varargin}); setappdata(0, 'liste_blocs', L);
add_block(lib, b, 'Position', [x y x + 75 y + 45], varargin{:});
end

function relier(mdl, a, b)
add_line(mdl, a, b, 'autorouting', 'on');
end

function fonction(mdl, nom, col, lig, code)
ajouter(mdl, 'simulink/User-Defined Functions/MATLAB Function', nom, col, lig);
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', [mdl '/' nom]); ch.Script = code;
end

function mdl = nouveau(nom, P, m)
mdl = ['mod_' nom];
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(1 / P.fs, 17), ...
    'StopTime', '1', 'SolverMode', 'SingleTasking');
hws = get_param(mdl, 'ModelWorkspace');
v = {'fs', P.fs; 'sps', P.sps; 'Ts', 1 / P.Rs; 'fc', P.fc; 'sigma2', 0; 'graine', 7; 'M', m.M; 'phi_ol', P.phi_ol; ...
     'df_ol', P.df_ol; 'ptsI', real(m.pts); 'ptsQ', imag(m.pts); 'points', m.ideal; 'ordre', m.ordre};
for k = 1:size(v, 1), hws.assignin(v{k, 1}, v{k, 2}); end
end

function emetteur_porteuse(mdl, entree, c0)
% lo_tx = exp(j 2 pi fc t), mix_tx, emis_reel = partie réelle, canal = emis_reel + bruit
B = 'simulink/';
ajouter(mdl, [B 'Sources/Sine Wave'], 'lo_tx cos', c0, 1.6, 'Amplitude', '1', 'Frequency', '2*pi*fc', 'Phase', 'pi/2', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Sources/Sine Wave'], 'lo_tx sin', c0, 2.4, 'Amplitude', '1', 'Frequency', '2*pi*fc', 'Phase', '0', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'lo_tx', c0 + 1, 2, 'Input', 'Real and imag');
ajouter(mdl, [B 'Math Operations/Product'], 'mix_tx', c0 + 2, 1, 'Inputs', '2');
ajouter(mdl, [B 'Math Operations/Complex to Real-Imag'], 'emis_reel', c0 + 3, 1, 'Output', 'Real');
ajouter(mdl, [B 'Sources/Random Number'], 'bruit', c0 + 3, 2.2, 'Mean', '0', 'Variance', 'sigma2', 'Seed', '1', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Sum'], 'canal', c0 + 4, 1.6, 'Inputs', '++');
ajouter(mdl, [B 'Sinks/To Workspace'], 's', c0 + 4, 0.4, 'VariableName', 's', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/Scope'], 'scope_emis', c0 + 5, 0.4, 'NumInputPorts', '1');
relier(mdl, 'lo_tx cos/1', 'lo_tx/1'); relier(mdl, 'lo_tx sin/1', 'lo_tx/2');
relier(mdl, entree, 'mix_tx/1'); relier(mdl, 'lo_tx/1', 'mix_tx/2'); relier(mdl, 'mix_tx/1', 'emis_reel/1');
relier(mdl, 'emis_reel/1', 'canal/1'); relier(mdl, 'bruit/1', 'canal/2'); relier(mdl, 'emis_reel/1', 's/1');
relier(mdl, 'emis_reel/1', 'scope_emis/1');
end

function oscillateur_rx(mdl, nom, freq, col, lig)
% nom = 2 exp(-j (2 pi freq t + phi_ol)) : oscillateur local du récepteur, décalé en phase et en fréquence
B = 'simulink/';
ajouter(mdl, [B 'Sources/Sine Wave'], [nom ' cos'], col, lig - 0.4, 'Amplitude', '2', 'Frequency', ['2*pi*(' freq ')'], 'Phase', 'phi_ol+pi/2', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Sources/Sine Wave'], [nom ' -sin'], col, lig + 0.4, 'Amplitude', '-2', 'Frequency', ['2*pi*(' freq ')'], 'Phase', 'phi_ol', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], nom, col + 1, lig, 'Input', 'Real and imag');
relier(mdl, [nom ' cos/1'], [nom '/1']); relier(mdl, [nom ' -sin/1'], [nom '/2']);
end

function mdl = construire_lineaire(m, P)
mdl = nouveau(m.nom, P, m); hws = get_param(mdl, 'ModelWorkspace');
rrc = strcmp(m.forme, 'rrc'); B = 'simulink/';
if rrc
    h = rrc_gr(P.sps, P.fs, P.Rs, P.alpha, P.ntaps); hr = rrc_gr(1, P.fs, P.Rs, P.alpha, P.ntaps);
    hws.assignin('h_tx', h); hws.assignin('h_rx', hr); gmf = max(conv(h, hr)); off = 0;
else
    hws.assignin('h_rx', ones(1, P.sps) / P.sps); gmf = 1; off = (P.sps - 1) / P.fs;
end
hws.assignin('gain_ech', 1 / (gmf * m.amp_moy)); hws.assignin('decalage', off);
hws.assignin('w_boucle', P.w_costas); if strcmp(m.detection, 'recepteur'), hws.assignin('w_boucle', P.w_recepteur); end
% émetteur
ajouter(mdl, [B 'Sources/Uniform Random Number'], 'src', 0, 1, 'Minimum', '0', 'Maximum', 'M', 'Seed', 'graine', 'SampleTime', 'Ts');
ajouter(mdl, [B 'Math Operations/Rounding Function'], 'indice', 1, 1, 'Operator', 'floor');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'table_I', 2, 0.6, 'NumberOfTableDimensions', '1', 'Table', 'ptsI');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'table_Q', 2, 1.4, 'NumberOfTableDimensions', '1', 'Table', 'ptsQ');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'symboles_tx', 3, 1, 'Input', 'Real and imag');
ajouter(mdl, [B 'Sinks/To Workspace'], 'sym', 2, 2.4, 'VariableName', 'sym', 'SaveFormat', 'Array');
relier(mdl, 'src/1', 'indice/1'); relier(mdl, 'indice/1', 'table_I/1'); relier(mdl, 'indice/1', 'table_Q/1'); relier(mdl, 'indice/1', 'sym/1');
relier(mdl, 'table_I/1', 'symboles_tx/1'); relier(mdl, 'table_Q/1', 'symboles_tx/2');
if rrc
    ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'maintien', 4, 1, 'OutPortSampleTime', '1/fs');
    ajouter(mdl, [B 'Sources/Pulse Generator'], 'impulsions', 4, 2.2, 'PulseType', 'Sample based', 'Period', 'sps', ...
        'PulseWidth', '1', 'PhaseDelay', '0', 'Amplitude', '1', 'SampleTime', '1/fs');
    ajouter(mdl, [B 'Math Operations/Product'], 'train', 5, 1, 'Inputs', '2');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], 'mise_en_forme', 6, 1, 'Numerator', 'h_tx', 'Denominator', '1');
    relier(mdl, 'symboles_tx/1', 'maintien/1'); relier(mdl, 'maintien/1', 'train/1'); relier(mdl, 'impulsions/1', 'train/2');
    relier(mdl, 'train/1', 'mise_en_forme/1'); c0 = 7;
else
    ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'mise_en_forme', 4, 1, 'OutPortSampleTime', '1/fs');
    relier(mdl, 'symboles_tx/1', 'mise_en_forme/1'); c0 = 5;
end
ajouter(mdl, [B 'Math Operations/Complex to Real-Imag'], 'bande_base_I', c0, 0, 'Output', 'Real');
ajouter(mdl, [B 'Sinks/To Workspace'], 'bb', c0 + 1, 0, 'VariableName', 'bb', 'SaveFormat', 'Array');
relier(mdl, 'mise_en_forme/1', 'bande_base_I/1'); relier(mdl, 'bande_base_I/1', 'bb/1');
emetteur_porteuse(mdl, 'mise_en_forme/1', c0);
% récepteur : mélange avec l'oscillateur local décalé, filtre adapté, normalisation
oscillateur_rx(mdl, 'lo_rx', 'fc+df_ol', c0, 4);
ajouter(mdl, [B 'Math Operations/Product'], 'mix_rx', c0 + 2, 3.4, 'Inputs', '2');
ajouter(mdl, [B 'Discrete/Discrete Filter'], 'filtre_adapte', c0 + 3, 3.4, 'Numerator', 'h_rx', 'Denominator', '1');
ajouter(mdl, [B 'Math Operations/Gain'], 'normalisation', c0 + 4, 3.4, 'Gain', 'gain_ech');
relier(mdl, 'canal/1', 'mix_rx/1'); relier(mdl, 'lo_rx/1', 'mix_rx/2'); relier(mdl, 'mix_rx/1', 'filtre_adapte/1');
relier(mdl, 'filtre_adapte/1', 'normalisation/1');
ech = '[Ts, decalage]'; c = c0 + 5;
ajouter(mdl, [B 'Sources/Constant'], 'points', c + 1, 6, 'Value', 'points');
switch m.detection
    case 'enveloppe'
        ajouter(mdl, [B 'Math Operations/Complex to Magnitude-Angle'], 'module', c, 3.4, 'Output', 'Magnitude');
        ajouter(mdl, [B 'Sinks/To Workspace'], 'y', c + 1, 2.4, 'VariableName', 'y', 'SaveFormat', 'Array');
        ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage', c + 1, 3.4, 'SampleTime', ech);
        ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage_c', c, 5, 'SampleTime', ech);
        ajouter(mdl, [B 'Sinks/To Workspace'], 'ech', c + 1, 5, 'VariableName', 'ech', 'SaveFormat', 'Array');
        relier(mdl, 'normalisation/1', 'module/1'); relier(mdl, 'module/1', 'y/1'); relier(mdl, 'module/1', 'echantillonnage/1');
        relier(mdl, 'normalisation/1', 'echantillonnage_c/1'); relier(mdl, 'echantillonnage_c/1', 'ech/1');
        if strcmp(m.nom, 'ook')
            ajouter(mdl, [B 'Logic and Bit Operations/Compare To Constant'], 'decision', c + 2, 3.4, 'relop', '>', 'const', '1.0');
            relier(mdl, 'echantillonnage/1', 'decision/1');
        else
            fonction(mdl, 'decision', c + 2, 3.4, code_decision());
            relier(mdl, 'echantillonnage/1', 'decision/1'); relier(mdl, 'points/1', 'decision/2');
        end
    case 'costas'
        ajouter(mdl, [B 'Sources/Constant'], 'ordre', c, 4.6, 'Value', 'ordre');
        ajouter(mdl, [B 'Sources/Constant'], 'w_boucle', c, 5.4, 'Value', 'w_boucle');
        fonction(mdl, 'costas', c + 1, 3.4, code_costas());
        ajouter(mdl, [B 'Sinks/To Workspace'], 'y', c + 2, 2.4, 'VariableName', 'y', 'SaveFormat', 'Array');
        ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage', c + 2, 3.4, 'SampleTime', ech);
        ajouter(mdl, [B 'Sinks/To Workspace'], 'ech', c + 3, 2.4, 'VariableName', 'ech', 'SaveFormat', 'Array');
        ajouter(mdl, [B 'Sinks/Terminator'], 'phase_costas', c + 2, 4.6);
        fonction(mdl, 'decision', c + 3, 3.4, code_decision());
        relier(mdl, 'normalisation/1', 'costas/1'); relier(mdl, 'ordre/1', 'costas/2'); relier(mdl, 'w_boucle/1', 'costas/3');
        relier(mdl, 'costas/1', 'y/1'); relier(mdl, 'costas/1', 'echantillonnage/1'); relier(mdl, 'costas/2', 'phase_costas/1');
        relier(mdl, 'echantillonnage/1', 'ech/1'); relier(mdl, 'echantillonnage/1', 'decision/1'); relier(mdl, 'points/1', 'decision/2');
    case 'recepteur'
        ajouter(mdl, [B 'Sinks/To Workspace'], 'y', c, 2.4, 'VariableName', 'y', 'SaveFormat', 'Array');
        ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage', c, 3.4, 'SampleTime', ech);
        ajouter(mdl, [B 'Sources/Constant'], 'w_boucle', c, 5.4, 'Value', 'w_boucle');
        fonction(mdl, 'recepteur', c + 2, 3.4, code_recepteur());
        ajouter(mdl, [B 'Sinks/To Workspace'], 'ech', c + 3, 2.4, 'VariableName', 'ech', 'SaveFormat', 'Array');
        ajouter(mdl, [B 'Sinks/To Workspace'], 'phase', c + 3, 4.4, 'VariableName', 'phase', 'SaveFormat', 'Array');
        relier(mdl, 'normalisation/1', 'y/1'); relier(mdl, 'normalisation/1', 'echantillonnage/1');
        relier(mdl, 'echantillonnage/1', 'recepteur/1'); relier(mdl, 'points/1', 'recepteur/2'); relier(mdl, 'w_boucle/1', 'recepteur/3');
        relier(mdl, 'recepteur/2', 'ech/1'); relier(mdl, 'recepteur/3', 'phase/1');
end
ajouter(mdl, [B 'Sinks/To Workspace'], 'dec', c + 4, 3.4, 'VariableName', 'dec', 'SaveFormat', 'Array');
relier(mdl, iff(strcmp(m.detection, 'recepteur'), 'recepteur/1', 'decision/1'), 'dec/1');
end

function mdl = construire_fsk(m, P)
mdl = nouveau(m.nom, P, m); hws = get_param(mdl, 'ModelWorkspace'); B = 'simulink/';
hws.assignin('excursion', P.df); hws.assignin('h_rx', ones(1, P.sps) / P.sps); hws.assignin('decalage', (P.sps - 1) / P.fs);
hws.assignin('h_pb', passe_bas_gr(P.fs, 3000, 1000));
ajouter(mdl, [B 'Sources/Uniform Random Number'], 'src', 0, 1, 'Minimum', '0', 'Maximum', 'M', 'Seed', 'graine', 'SampleTime', 'Ts');
ajouter(mdl, [B 'Math Operations/Rounding Function'], 'indice', 1, 1, 'Operator', 'floor');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'nrz', 2, 1, 'NumberOfTableDimensions', '1', 'Table', 'ptsI');
ajouter(mdl, [B 'Sinks/To Workspace'], 'sym', 2, 2.4, 'VariableName', 'sym', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'mise_en_forme', 3, 1, 'OutPortSampleTime', '1/fs');
ajouter(mdl, [B 'Sinks/To Workspace'], 'bb', 4, 0, 'VariableName', 'bb', 'SaveFormat', 'Array');
% fsk_mod (équivalent du bloc Frequency Mod) : phase += 2 pi excursion / fs * x, puis exp(j phase)
ajouter(mdl, [B 'Math Operations/Gain'], 'sensibilite', 4, 1, 'Gain', '2*pi*excursion/fs');
ajouter(mdl, [B 'Discrete/Discrete-Time Integrator'], 'accumulateur', 5, 1, 'IntegratorMethod', 'Accumulation: Forward Euler', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Trigonometric Function'], 'cos', 6, 0.6, 'Operator', 'cos');
ajouter(mdl, [B 'Math Operations/Trigonometric Function'], 'sin', 6, 1.4, 'Operator', 'sin');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'fsk_mod', 7, 1, 'Input', 'Real and imag');
relier(mdl, 'src/1', 'indice/1'); relier(mdl, 'indice/1', 'nrz/1'); relier(mdl, 'indice/1', 'sym/1');
relier(mdl, 'nrz/1', 'mise_en_forme/1'); relier(mdl, 'mise_en_forme/1', 'bb/1'); relier(mdl, 'mise_en_forme/1', 'sensibilite/1');
relier(mdl, 'sensibilite/1', 'accumulateur/1'); relier(mdl, 'accumulateur/1', 'cos/1'); relier(mdl, 'accumulateur/1', 'sin/1');
relier(mdl, 'cos/1', 'fsk_mod/1'); relier(mdl, 'sin/1', 'fsk_mod/2');
emetteur_porteuse(mdl, 'fsk_mod/1', 8);
% deux corrélateurs non cohérents (bit 1 : fc + excursion, bit 0 : fc - excursion)
noms = {'1', '+', 3.6; '0', '-', 5.6};
for t = 1:2
    n = noms{t, 1}; l = noms{t, 3};
    oscillateur_rx(mdl, ['lo_' n], ['fc' noms{t, 2} 'excursion+df_ol'], 8, l);
    ajouter(mdl, [B 'Math Operations/Product'], ['mix_' n], 10, l, 'Inputs', '2');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], ['correlateur_' n], 11, l, 'Numerator', 'h_rx', 'Denominator', '1');
    ajouter(mdl, [B 'Math Operations/Complex to Magnitude-Angle'], ['module_' n], 12, l, 'Output', 'Magnitude');
    relier(mdl, 'canal/1', ['mix_' n '/1']); relier(mdl, ['lo_' n '/1'], ['mix_' n '/2']);
    relier(mdl, ['mix_' n '/1'], ['correlateur_' n '/1']); relier(mdl, ['correlateur_' n '/1'], ['module_' n '/1']);
end
ajouter(mdl, [B 'Math Operations/Sum'], 'difference', 13, 4.6, 'Inputs', '+-');
ajouter(mdl, [B 'Sinks/To Workspace'], 'y', 14, 3.8, 'VariableName', 'y', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage', 14, 4.6, 'SampleTime', '[Ts, decalage]');
ajouter(mdl, [B 'Logic and Bit Operations/Compare To Zero'], 'decision', 15, 4.6, 'relop', '>');
ajouter(mdl, [B 'Sinks/To Workspace'], 'dec', 16, 4.6, 'VariableName', 'dec', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'paire', 13, 6.4, 'Input', 'Real and imag');
ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage_paire', 14, 6.4, 'SampleTime', '[Ts, decalage]');
ajouter(mdl, [B 'Sinks/To Workspace'], 'ech', 15, 6.4, 'VariableName', 'ech', 'SaveFormat', 'Array');
relier(mdl, 'module_1/1', 'difference/1'); relier(mdl, 'module_0/1', 'difference/2'); relier(mdl, 'difference/1', 'y/1');
relier(mdl, 'difference/1', 'echantillonnage/1'); relier(mdl, 'echantillonnage/1', 'decision/1'); relier(mdl, 'decision/1', 'dec/1');
relier(mdl, 'module_1/1', 'paire/1'); relier(mdl, 'module_0/1', 'paire/2'); relier(mdl, 'paire/1', 'echantillonnage_paire/1');
relier(mdl, 'echantillonnage_paire/1', 'ech/1');
% discriminateur de fréquence (équivalent du bloc Quadrature Demod), comme dans la leçon Learn SDR sur la FSK
oscillateur_rx(mdl, 'lo_rx', 'fc+df_ol', 8, 8);
ajouter(mdl, [B 'Math Operations/Product'], 'mix_rx', 10, 8, 'Inputs', '2');
ajouter(mdl, [B 'Discrete/Discrete Filter'], 'passe_bas', 11, 8, 'Numerator', 'h_pb', 'Denominator', '1');
ajouter(mdl, [B 'Discrete/Unit Delay'], 'retard', 12, 8.8, 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Math Function'], 'conj', 13, 8.8, 'Operator', 'conj');
ajouter(mdl, [B 'Math Operations/Product'], 'produit', 14, 8, 'Inputs', '2');
ajouter(mdl, [B 'Math Operations/Complex to Magnitude-Angle'], 'argument', 15, 8, 'Output', 'Angle');
ajouter(mdl, [B 'Math Operations/Gain'], 'discriminateur', 16, 8, 'Gain', 'fs/(2*pi*excursion)');
ajouter(mdl, [B 'Discrete/Discrete Filter'], 'integrateur', 17, 8, 'Numerator', 'h_rx', 'Denominator', '1');
ajouter(mdl, [B 'Discrete/Zero-Order Hold'], 'echantillonnage_d', 18, 8, 'SampleTime', '[Ts, decalage]');
ajouter(mdl, [B 'Logic and Bit Operations/Compare To Zero'], 'decision_d', 19, 8, 'relop', '>');
ajouter(mdl, [B 'Sinks/To Workspace'], 'dec_d', 20, 8, 'VariableName', 'dec_d', 'SaveFormat', 'Array');
relier(mdl, 'canal/1', 'mix_rx/1'); relier(mdl, 'lo_rx/1', 'mix_rx/2'); relier(mdl, 'mix_rx/1', 'passe_bas/1');
relier(mdl, 'passe_bas/1', 'produit/1'); relier(mdl, 'passe_bas/1', 'retard/1'); relier(mdl, 'retard/1', 'conj/1');
relier(mdl, 'conj/1', 'produit/2'); relier(mdl, 'produit/1', 'argument/1'); relier(mdl, 'argument/1', 'discriminateur/1');
relier(mdl, 'discriminateur/1', 'integrateur/1'); relier(mdl, 'integrateur/1', 'echantillonnage_d/1');
relier(mdl, 'echantillonnage_d/1', 'decision_d/1'); relier(mdl, 'decision_d/1', 'dec_d/1');
end

function r = iff(c, a, b), if c, r = a; else, r = b; end, end

function c = code_decision()
c = sprintf(['function indice = decision(x, points)\n' ...
    '%% decision au plus proche voisin, comme le bloc Constellation Decoder de GNU Radio\n' ...
    '[~, k] = min(abs(x - points));\nindice = k - 1;\nend\n']);
end

function c = code_costas()
c = sprintf(['function [y, phase] = costas(x, ordre, w)\n' ...
    '%% boucle de Costas, meme algorithme que le bloc Costas Loop de GNU Radio 3.10\n' ...
    'persistent ph fr\nif isempty(ph), ph = 0; fr = 0; end\n' ...
    'z = sqrt(2) / 2; den = 1 + 2 * z * w + w * w; a = 4 * z * w / den; b = 4 * w * w / den;\n' ...
    'y = x * exp(-1j * ph); r = real(y); q = imag(y); K = sqrt(2) - 1;\n' ...
    'if ordre == 2\n    e = r * q;\nelseif ordre == 4\n    e = sg(r) * q - sg(q) * r;\n' ...
    'elseif abs(r) >= abs(q)\n    e = sg(r) * q - K * sg(q) * r;\nelse\n    e = K * sg(r) * q - sg(q) * r;\nend\n' ...
    'e = max(-1, min(1, e)); fr = fr + b * e; ph = ph + fr + a * e;\n' ...
    'ph = mod(ph + pi, 2 * pi) - pi; fr = max(-1, min(1, fr)); phase = ph;\nend\n' ...
    'function s = sg(v)\nif v > 0, s = 1; else, s = -1; end\nend\n']);
end

function c = code_recepteur()
c = sprintf(['function [indice, y, phase] = recepteur(x, points, w)\n' ...
    '%% recepteur de constellation : boucle de phase pilotee par les decisions, comme le bloc\n' ...
    '%% Constellation Receiver de GNU Radio (erreur = angle entre l''echantillon et le point decide)\n' ...
    'persistent ph fr\nif isempty(ph), ph = 0; fr = 0; end\n' ...
    'z = sqrt(2) / 2; den = 1 + 2 * z * w + w * w; a = 4 * z * w / den; b = 4 * w * w / den;\n' ...
    'y = x * exp(-1j * ph); [~, k] = min(abs(y - points)); e = angle(y * conj(points(k)));\n' ...
    'fr = fr + b * e; ph = ph + fr + a * e; ph = mod(ph + pi, 2 * pi) - pi; fr = max(-0.1, min(0.1, fr));\n' ...
    'indice = k - 1; phase = ph;\nend\n']);
end

function h = passe_bas_gr(fs, fcoup, largeur)
% même filtre que firdes.low_pass(1, fs, fcoup, largeur, fenêtre de Hamming) de GNU Radio
n = floor(53 * fs / (22 * largeur)); n = bitor(n, 1); M = (n - 1) / 2; k = -M:M;
w = 0.54 - 0.46 * cos(2 * pi * (0:n - 1) / (n - 1)); wc = 2 * pi * fcoup / fs;
h = zeros(1, n); h(k ~= 0) = sin(wc * k(k ~= 0)) ./ (pi * k(k ~= 0)); h(k == 0) = wc / pi;
h = h .* w; h = h / sum(h);
end

function h = rrc_gr(gain, fs, rs, alpha, ntaps)
% même formule que firdes.root_raised_cosine de GNU Radio (somme des coefficients = gain)
spb = fs / rs; h = zeros(1, ntaps); c = floor(ntaps / 2);
for i = 0:ntaps - 1
    x = i - c; x1 = pi * x / spb; x2 = 4 * alpha * x / spb; x3 = x2 * x2 - 1;
    if abs(x3) >= 1e-6
        if i ~= c, num = cos((1 + alpha) * x1) + sin((1 - alpha) * x1) / (4 * alpha * x / spb);
        else, num = cos((1 + alpha) * x1) + (1 - alpha) * pi / (4 * alpha); end
        den = x3 * pi;
    else
        a3 = (1 - alpha) * x1; a2 = (1 + alpha) * x1;
        num = sin(a2) * (1 + alpha) * pi - cos(a3) * ((1 - alpha) * pi * spb) / (4 * alpha * x) + sin(a3) * spb * spb / (4 * alpha * x * x);
        den = -32 * pi * alpha * alpha * x / spb;
    end
    h(i + 1) = 4 * alpha * num / den;
end
h = h * gain / sum(h);
end

%% ------------------------------------------------------------------ exécution et mesures
function so = simuler(mdl, nsym, sigma2, P)
hws = get_param(mdl, 'ModelWorkspace'); hws.assignin('sigma2', sigma2);
so = sim(mdl, 'ReturnWorkspaceOutputs', 'on', 'StopTime', num2str((nsym - 1) / P.Rs + (P.sps - 1) / P.fs, 17));
end

function v = lire(so, nom), v = so.get(nom); v = v(:).'; end

function P = perms_ambiguite(m)
P = struct('k', {}, 'perm', {}, 'rot', {});
for k = 0:m.ordre - 1
    r = exp(2j * pi * k / m.ordre); [~, perm] = min(abs((m.pts(:) * r) - m.pts(:).'), [], 2); perm = perm(:).' - 1;
    if numel(unique(perm)) == m.M, P(end + 1) = struct('k', k, 'perm', perm, 'rot', r); end %#ok<AGROW>
end
end

function A = aligner(m, sym, dec, P)
% retard (en symboles) et rotation d'ambiguïté de la boucle, levée sur 500 symboles connus après l'accrochage
best = -1; c = P.conv;
for p = perms_ambiguite(m)
    for L = 0:59
        if c + L + 500 > numel(dec), continue; end
        acc = mean(p.perm(sym(c + 1:c + 500) + 1) == dec(c + L + 1:c + L + 500));
        if acc > best, best = acc; A = struct('L', L, 'k', p.k, 'perm', p.perm, 'rot', p.rot); end
    end
end
end

function [es, ns, eb, nb, gl] = compter(m, sym, dec, A, P)
% erreurs comptées après l'accrochage ; avec une boucle, les 100 premiers symboles de chaque bloc de 1000
% sont connus (pilotes) et lèvent l'ambiguïté du bloc après un éventuel glissement de cycle
n = min(numel(sym), numel(dec) - A.L) - 5; aa = sym(P.conv + 1:n); ba = dec(P.conv + 1 + A.L:n + A.L);
if m.ordre > 1, pp = perms_ambiguite(m); else, pp = struct('perm', A.perm); end
es = 0; ns = 0; eb = 0; gl = 0; prec = 0; BLOC = 1000; PIL = 100;
for i = 1:BLOC:numel(aa) - BLOC + 1
    a = aa(i:i + BLOC - 1); b = ba(i:i + BLOC - 1); sc = zeros(1, numel(pp));
    for j = 1:numel(pp), sc(j) = mean(pp(j).perm(a(1:PIL) + 1) == b(1:PIL)); end
    [~, j] = max(sc); if prec > 0 && j ~= prec, gl = gl + 1; end, prec = j;
    inv = zeros(1, m.M); inv(pp(j).perm + 1) = 0:m.M - 1; bb = inv(b(PIL + 1:end) + 1); a = a(PIL + 1:end);
    es = es + sum(a ~= bb); ns = ns + numel(a); x = bitxor(m.lab(a + 1), m.lab(bb + 1));
    for v = x, eb = eb + sum(dec2bin(v) == '1'); end
end
nb = ns * m.bits;
end

function p = theorie(m, g)
switch m.nom
    case 'ook', p = ber_ook_enveloppe(g);
    case 'ask4', p = ser_ask_enveloppe(abs(m.pts), g);
    case {'bpsk', 'qpsk', 'psk8'}, M = m.M; gs = g * log2(M); a = sin(pi / M) ^ 2;
        p = integral(@(t) exp(-gs * a ./ sin(t) .^ 2), 0, (M - 1) * pi / M) / pi;
    case 'qam16', gs = 4 * g; q = 2 * (1 - 1 / 4) * Qf(sqrt(3 * gs / 15)); p = 1 - (1 - q) ^ 2;
    case 'fsk2', p = 0.5 * exp(-g / 2);
end
end

function q = Qf(x), q = 0.5 * erfc(x / sqrt(2)); end

function q = marcum(a, b)
% Q1(a, b) de Marcum, par intégration de la loi de Rice
q = integral(@(x) x .* exp(-(x - a) .^ 2 / 2) .* besseli(0, a * x, 1), b, b + a + 40);
end

function p = ber_ook_enveloppe(g)
p = 0.5 * (exp(-g / 2) + 1 - marcum(2 * sqrt(g), sqrt(g)));
end

function p = ser_ask_enveloppe(niv, g)
% M-ASK, détection d'enveloppe, seuils au milieu des niveaux
M = numel(niv); gs = g * log2(M); s = sqrt(mean(niv .^ 2) / (2 * gs)); t = (niv(1:end - 1) + niv(2:end)) / 2; p = 0;
for i = 1:M
    if i > 1, bas = marcum(niv(i) / s, t(i - 1) / s); else, bas = 1; end
    if i < M, haut = marcum(niv(i) / s, t(i) / s); else, haut = 0; end
    p = p + 1 - (bas - haut);
end
p = p / M;
end

function [f, Pn] = densite(s, fs)
N = 4096; w = 0.5 - 0.5 * cos(2 * pi * (0:N - 1) / N); P = zeros(1, N / 2 + 1); k = 0;
for i = 1:N / 2:numel(s) - N
    X = fft(s(i:i + N - 1) .* w); P = P + abs(X(1:N / 2 + 1)) .^ 2; k = k + 1;
end
P = P / k; Pn = P / max(P); f = (0:N / 2) * fs / N;
end

function [r, echecs] = verifier(mdl, m, P, ici)
echecs = {}; r = struct('titre', m.titre); boucle = any(strcmp(m.detection, {'costas', 'recepteur'}));
so = simuler(mdl, 4000, 0, P);
s = lire(so, 's'); Ps = mean(s(20 * P.sps:end) .^ 2); r.puissance_signal = Ps;
sym = round(lire(so, 'sym')); dec = round(double(lire(so, 'dec'))); A = aligner(m, sym, dec, P);
[es, ns] = compter(m, sym, dec, A, P);
r.sans_bruit = struct('retard_symboles', A.L, 'rotation_ambiguite_deg', 360 * A.k / m.ordre, 'erreurs_symboles', es, 'symboles', ns);
if strcmp(m.nom, 'fsk2')
    dd = round(double(lire(so, 'dec_d'))); Ad = aligner(m, sym, dd, P); r.sans_bruit.erreurs_discriminateur = compter(m, sym, dd, Ad, P);
    if r.sans_bruit.erreurs_discriminateur > 0, echecs{end + 1} = sprintf('%s : erreurs sans bruit au discriminateur', m.nom); end
else
    ech = lire(so, 'ech'); n2 = min(numel(ech) - A.L, numel(sym)) - 5;
    ref = m.ideal(sym(P.conv + 1:n2) + 1) * A.rot; mes = ech(P.conv + 1 + A.L:n2 + A.L);
    if strcmp(m.detection, 'enveloppe'), ref = abs(ref); mes = abs(mes); end
    g = real(mes * ref') / real(ref * ref'); evm = sqrt(mean(abs(mes - ref) .^ 2) / mean(abs(ref) .^ 2));
    r.sans_bruit.gain = g; r.sans_bruit.evm = evm;
    if evm > 0.05 || abs(g - 1) > 0.05, echecs{end + 1} = sprintf('%s : EVM %.3f, gain %.3f', m.nom, evm, g); end
end
fprintf('Sans bruit : retard %d symboles, rotation %g degres, %d erreurs sur %d symboles\n', A.L, 360 * A.k / m.ordre, es, ns);
if es > 0, echecs{end + 1} = sprintf('%s : %d erreurs sans bruit', m.nom, es); end
[f, Pn] = densite(s, P.fs); bb = lire(so, 'bb');
if strcmp(m.forme, 'rrc')
    B = (1 + P.alpha) * P.Rs / 2; frac = sum(Pn(f > P.fc - B & f < P.fc + B)) / sum(Pn);
    r.spectre = struct('part_puissance_bande', frac); fprintf('Spectre : %.4f de la puissance entre %g et %g Hz\n', frac, P.fc - B, P.fc + B);
    if frac < 0.98, echecs{end + 1} = sprintf('%s : %.3f de la puissance dans la bande', m.nom, frac); end
elseif strcmp(m.nom, 'fsk2')
    [~, i1] = max(Pn .* (f < P.fc)); [~, i2] = max(Pn .* (f > P.fc));
    r.spectre = struct('raies_Hz', [f(i1), f(i2)]); fprintf('Spectre : raies a %.1f et %.1f Hz\n', f(i1), f(i2));
    if abs(f(i1) - (P.fc - P.df)) > 20 || abs(f(i2) - (P.fc + P.df)) > 20, echecs{end + 1} = sprintf('%s : raies %g %g', m.nom, f(i1), f(i2)); end
else
    lv = 10 * log10([max(Pn(abs(f - (P.fc - P.Rs)) < 20)), max(Pn(abs(f - (P.fc + P.Rs)) < 20))]);
    r.spectre = struct('niveau_premiers_zeros_dB', lv); fprintf('Spectre : premiers zeros a %.1f et %.1f dB\n', lv);
    if max(lv) > -20, echecs{end + 1} = sprintf('%s : zeros du spectre a %.1f dB', m.nom, max(lv)); end
end
tol = 1.5; if boucle, tol = 2; end
taux = struct('EbN0_dB', {}, 'variance_bruit', {}, 'TES', {}, 'TEB', {}, 'erreurs', {}, 'mesure', {}, 'theorie', {}, 'conforme', {}, 'TEB_discriminateur', {}, 'glissements_de_cycle', {}, 'hors_accrochage', {});
for e = m.ebn0
    g = 10 ^ (e / 10); s2 = Ps * P.sps / (2 * m.bits * g);
    so2 = simuler(mdl, 16000, s2, P); sym2 = round(lire(so2, 'sym')); dec2 = round(double(lire(so2, 'dec')));
    A2 = aligner(m, sym2, dec2, P); [es, ns, eb, nb, gl] = compter(m, sym2, dec2, A2, P);
    if strcmp(m.compare, 'TEB'), mes = eb / nb; ne = eb; else, mes = es / ns; ne = es; end
    th = theorie(m, g); ok = ne < 50 || (mes / th > 1 / tol && mes / th < tol); td = NaN;
    % boucle qui glisse : point rapporte hors accrochage ; aucun glissement admis lorsque la theorie est sous 1 %
    hc = boucle && gl > 0; if hc, ok = th > 0.01; end
    if strcmp(m.nom, 'fsk2')
        dd = round(double(lire(so2, 'dec_d'))); Ad = aligner(m, sym2, dd, P); [~, ~, ebd, nbd] = compter(m, sym2, dd, Ad, P); td = ebd / nbd;
    end
    taux(end + 1) = struct('EbN0_dB', e, 'variance_bruit', s2, 'TES', es / ns, 'TEB', eb / nb, 'erreurs', ne, ...
        'mesure', mes, 'theorie', th, 'conforme', ok, 'TEB_discriminateur', td, 'glissements_de_cycle', gl, 'hors_accrochage', hc); %#ok<AGROW>
    fprintf('Eb/N0 = %4.1f dB : %s mesure %.3e (%d erreurs), theorie %.3e, %d glissements  %s\n', e, m.compare, mes, ne, th, gl, string(ok));
    if ~ok, echecs{end + 1} = sprintf('%s : Eb/N0 %g dB, mesure %.3g, theorie %.3g', m.nom, e, mes, th); end %#ok<AGROW>
    if e == m.ebn0(3), xmid = lire(so2, 'ech'); end
end
r.taux = taux;
% figure : signal émis, spectre, constellation, oeil
so3 = simuler(mdl, 1500, Ps * P.sps / (2 * m.bits * 100), P); y = lire(so3, 'y');
hws = get_param(mdl, 'ModelWorkspace'); ph0 = mod(round(P.fs * hws.getVariable('decalage')), P.sps);
switch m.detection
    case 'costas'
        A3 = aligner(m, round(lire(so3, 'sym')), round(double(lire(so3, 'dec'))), P); oe = real(y * conj(A3.rot)); to = 'oeil de la voie I apres la boucle de Costas';
    case 'recepteur'
        ph = lire(so3, 'phase'); j = floor(((0:numel(y) - 1) - ph0) / P.sps); best = inf;
        for q = -20:20
            z = y .* exp(-1j * ph(min(max(j + q + 1, 1), numel(ph))));
            sm = z(ph0 + 1 + P.sps * (P.conv:P.conv + 700));
            for k = 0:3
                ro = exp(2j * pi * k / 4); d = mean(min(abs(sm(:) * conj(ro) - m.ideal(:).'), [], 2));
                if d < best, best = d; oe = real(z * conj(ro)); end
            end
        end
        to = 'oeil de la voie I corrigee par la phase de la boucle';
    otherwise
        oe = real(y); to = 'oeil du module';
        if strcmp(m.nom, 'fsk2'), to = 'oeil de |y_1| - |y_0|'; end
end
fh = figure('Visible', 'off', 'Position', [0 0 1300 900]);
subplot(2, 2, 1); i0 = 40 * P.sps; t = (0:8 * P.sps - 1) / P.fs * 1e3;
plot(t, s(i0 + 1:i0 + 8 * P.sps)); hold on; plot(t, bb(i0 + 1:i0 + 8 * P.sps), '--', 'LineWidth', 1.4); grid on;
xlabel('temps (ms)'); title([m.titre ' : signal emis et voie I en bande de base']);
subplot(2, 2, 2); plot(f / 1e3, 10 * log10(Pn + 1e-16)); xlim([0 10]); ylim([-80 3]); grid on;
xlabel('frequence (kHz)'); ylabel('dB'); title('spectre du signal emis');
subplot(2, 2, 3); x0 = lire(so, 'ech'); x0 = x0(P.conv + 1:end); xm = xmid(P.conv + 1:end);
plot(real(xm), imag(xm), '.', 'MarkerSize', 3); hold on; plot(real(x0), imag(x0), 'r+', 'MarkerSize', 14, 'LineWidth', 2); axis equal; grid on;
switch m.detection
    case 'enveloppe', title(sprintf('constellation avant detection d''enveloppe, Eb/N0 = %g dB', m.ebn0(3)));
    case 'energie', title(sprintf('points de decision (|y_1|, |y_0|), Eb/N0 = %g dB', m.ebn0(3)));
    otherwise, title(sprintf('constellation apres la boucle, Eb/N0 = %g dB', m.ebn0(3)));
end
subplot(2, 2, 4); hold on;
for j = P.conv + 60:P.conv + 260
    a = j * P.sps + ph0 - P.sps / 2; plot((0:2 * P.sps) / P.sps - 0.5, oe(a + 1:a + 2 * P.sps + 1), 'Color', [0 0.45 0.74]);
end
grid on; xlabel('temps (symboles)'); title([to ', Eb/N0 = 20 dB']);
sgtitle(['mod\_' m.nom '  -  MATLAB R' version('-release')]);
exportgraphics(fh, fullfile(ici, ['mod_' m.nom '_resultat.png']), 'Resolution', 130); close(fh);
end

function figure_taux(mods, R, ici)
fh = figure('Visible', 'off', 'Position', [0 0 1000 700]); co = lines(numel(mods));
for i = 1:numel(mods)
    m = mods(i); t = R.(m.nom).taux; e = [t.EbN0_dB]; y = [t.mesure]; y(y == 0) = NaN;
    ee = linspace(min(e) - 1, max(e) + 1, 50); th = arrayfun(@(x) theorie(m, 10 ^ (x / 10)), ee);
    semilogy(ee, th, '-', 'Color', co(i, :)); hold on;
    hc = logical([t.hors_accrochage]);
    semilogy(e(~hc), y(~hc), 'o', 'Color', co(i, :), 'DisplayName', [m.titre ' (' m.compare ')']);
    if any(hc), semilogy(e(hc), y(hc), 'x', 'Color', co(i, :), 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off'); end
end
semilogy(NaN, NaN, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'boucle hors accrochage (glissements de cycle)');
ylim([1e-5 0.5]); grid on; xlabel('Eb/N0 (dB)'); ylabel('taux d''erreur');
title(['Taux d''erreur mesures dans Simulink (points) et theorie (traits) - MATLAB R' version('-release')]);
h = [flipud(findobj(gca, 'Marker', 'o')); findobj(gca, 'Marker', 'x', 'HandleVisibility', 'on')]; legend(h, 'Location', 'southwest');
exportgraphics(fh, fullfile(ici, 'taux_erreur_simulink.png'), 'Resolution', 130); close(fh);
end
