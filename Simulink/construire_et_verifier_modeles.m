%% CONSTRUIRE_ET_VERIFIER_MODELES - Thème 4 : modulations numériques
% Construit, avec les seuls blocs de la bibliothèque de base de Simulink, les modèles
% mod_ook, mod_ask4, mod_bpsk, mod_qpsk, mod_psk8, mod_qam16 et mod_fsk2 (simulation pure,
% sans matériel), les exécute sans bruit puis pour plusieurs rapports Eb/N0, et vérifie :
% absence d'erreur sans bruit, constellation, spectre, diagramme de l'oeil et taux d'erreur
% comparés à la théorie. Les mêmes paramètres figurent dans En_ligne/modulations.py.
clear; clc;
ici = fileparts(mfilename('fullpath')); cd(ici);
if exist(fullfile(ici, 'journal_matlab.txt'), 'file'), delete(fullfile(ici, 'journal_matlab.txt')); end
diary(fullfile(ici, 'journal_matlab.txt'));
fprintf('MATLAB %s\n', version);
P = struct('fs', 32000, 'sps', 32, 'Rs', 1000, 'fc', 4000, 'alpha', 0.35, 'ntaps', 353, 'df', 500);
mods = definir_modulations();
echecs = {}; R = struct();
for i = 1:numel(mods)
    m = mods(i);
    fprintf('\n===== mod_%s : %s =====\n', m.nom, m.titre);
    if strcmp(m.nom, 'fsk2'), mdl = construire_fsk(m, P); else, mdl = construire_lineaire(m, P); end
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
 'ook',   'OOK (tout ou rien)', [0 1],               [0 1],     'rect', 'enveloppe', [6 8 10 12 14], 'TEB'
 'ask4',  'ASK a 4 niveaux',    [0 1/3 2/3 1],       g(0:3),    'rect', 'coherente', [8 10 12 14 16], 'TES'
 'bpsk',  'BPSK',               [1 -1],              [0 1],     'rrc',  'coherente', [0 2 4 6 8],    'TES'
 'qpsk',  'QPSK',               psk(4, pi/4),        g(0:3),    'rrc',  'coherente', [0 2 4 6 8],    'TES'
 'psk8',  '8-PSK',              psk(8, 0),           g(0:7),    'rrc',  'coherente', [4 6 8 10 12],  'TES'
 'qam16', '16-QAM',             q16,                 l16,       'rrc',  'coherente', [4 6 8 10 12],  'TES'
 'fsk2',  'FSK binaire',        [-1 1],              [0 1],     'rect', 'energie',   [4 6 8 10 12],  'TEB'};
for k = 1:size(D, 1)
    mods(k) = struct('nom', D{k, 1}, 'titre', D{k, 2}, 'pts', D{k, 3}, 'lab', D{k, 4}, 'forme', D{k, 5}, ...
        'detection', D{k, 6}, 'ebn0', D{k, 7}, 'compare', D{k, 8}, 'M', numel(D{k, 3}), 'bits', log2(numel(D{k, 3}))); %#ok<AGROW>
end
end

%% ------------------------------------------------------------------ construction
function b = ajouter(mdl, lib, nom, col, lig, varargin)
x = 40 + 150 * col; y = 40 + 95 * lig;
b = [mdl '/' nom];
add_block(lib, b, 'Position', [x y x + 70 y + 50], varargin{:});
end

function relier(mdl, a, b)
add_line(mdl, a, b, 'autorouting', 'on');
end

function mdl = nouveau(nom, P, m)
mdl = ['mod_' nom];
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(1 / P.fs, 17), ...
    'StopTime', '1', 'SolverMode', 'SingleTasking');
hws = get_param(mdl, 'ModelWorkspace');
hws.assignin('fs', P.fs); hws.assignin('sps', P.sps); hws.assignin('Ts', 1 / P.Rs); hws.assignin('fc', P.fc);
hws.assignin('sigma2', 0); hws.assignin('graine', 7);
hws.assignin('ptsI', real(m.pts)); hws.assignin('ptsQ', imag(m.pts)); hws.assignin('M', m.M);
end

function mdl = construire_lineaire(m, P)
mdl = nouveau(m.nom, P, m); hws = get_param(mdl, 'ModelWorkspace');
rrc = strcmp(m.forme, 'rrc');
if rrc
    h = rrc_gr(P.sps, P.fs, P.Rs, P.alpha, P.ntaps); hws.assignin('h_tx', h);
    hws.assignin('h_rx', h / sum(h));
else
    hws.assignin('h_rx', ones(1, P.sps) / P.sps);
end
B = 'simulink/';
ajouter(mdl, [B 'Sources/Uniform Random Number'], 'Source de symboles', 0, 1, 'Minimum', '0', 'Maximum', 'M', 'Seed', 'graine', 'SampleTime', 'Ts');
ajouter(mdl, [B 'Math Operations/Rounding Function'], 'Indice', 1, 1, 'Operator', 'floor');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'Table I', 2, 0, 'NumberOfTableDimensions', '1', 'Table', 'ptsI');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'Table Q', 2, 2, 'NumberOfTableDimensions', '1', 'Table', 'ptsQ');
ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'Maintien I', 3, 0, 'OutPortSampleTime', '1/fs');
ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'Maintien Q', 3, 2, 'OutPortSampleTime', '1/fs');
ajouter(mdl, [B 'Sinks/To Workspace'], 'sym', 2, 3.2, 'VariableName', 'sym', 'SaveFormat', 'Array');
relier(mdl, 'Source de symboles/1', 'Indice/1'); relier(mdl, 'Indice/1', 'Table I/1'); relier(mdl, 'Indice/1', 'Table Q/1');
relier(mdl, 'Indice/1', 'sym/1');
relier(mdl, 'Table I/1', 'Maintien I/1'); relier(mdl, 'Table Q/1', 'Maintien Q/1');
if rrc
    ajouter(mdl, [B 'Sources/Pulse Generator'], 'Impulsions', 3, 1, 'PulseType', 'Sample based', 'Period', 'sps', ...
        'PulseWidth', '1', 'PhaseDelay', '0', 'Amplitude', '1', 'SampleTime', '1/fs');
    ajouter(mdl, [B 'Math Operations/Product'], 'Train I', 4, 0, 'Inputs', '2');
    ajouter(mdl, [B 'Math Operations/Product'], 'Train Q', 4, 2, 'Inputs', '2');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], 'Mise en forme I', 5, 0, 'Numerator', 'h_tx', 'Denominator', '1');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], 'Mise en forme Q', 5, 2, 'Numerator', 'h_tx', 'Denominator', '1');
    relier(mdl, 'Maintien I/1', 'Train I/1'); relier(mdl, 'Impulsions/1', 'Train I/2');
    relier(mdl, 'Maintien Q/1', 'Train Q/1'); relier(mdl, 'Impulsions/1', 'Train Q/2');
    relier(mdl, 'Train I/1', 'Mise en forme I/1'); relier(mdl, 'Train Q/1', 'Mise en forme Q/1');
    sI = 'Mise en forme I/1'; sQ = 'Mise en forme Q/1'; c0 = 6;
else
    sI = 'Maintien I/1'; sQ = 'Maintien Q/1'; c0 = 4;
end
ajouter(mdl, [B 'Sources/Sine Wave'], 'Porteuse cos', c0, 1, 'Amplitude', '1', 'Frequency', '2*pi*fc', 'Phase', 'pi/2', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Sources/Sine Wave'], 'Porteuse -sin', c0, 3, 'Amplitude', '-1', 'Frequency', '2*pi*fc', 'Phase', '0', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Product'], 'Modulateur I', c0 + 1, 0, 'Inputs', '2');
ajouter(mdl, [B 'Math Operations/Product'], 'Modulateur Q', c0 + 1, 2, 'Inputs', '2');
ajouter(mdl, [B 'Math Operations/Sum'], 'Signal emis', c0 + 2, 1, 'Inputs', '++');
ajouter(mdl, [B 'Sources/Random Number'], 'Bruit', c0 + 2, 3, 'Mean', '0', 'Variance', 'sigma2', 'Seed', '1', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Sum'], 'Canal', c0 + 3, 2, 'Inputs', '++');
ajouter(mdl, [B 'Sinks/To Workspace'], 's', c0 + 3, 0, 'VariableName', 's', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'bb', c0 + 1, -1, 'VariableName', 'bb', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Signal Routing/Mux'], 'Mux', c0 + 4, 0, 'Inputs', '2');
ajouter(mdl, [B 'Sinks/Scope'], 'Scope', c0 + 5, 0, 'NumInputPorts', '1');
relier(mdl, sI, 'Modulateur I/1'); relier(mdl, 'Porteuse cos/1', 'Modulateur I/2');
relier(mdl, sQ, 'Modulateur Q/1'); relier(mdl, 'Porteuse -sin/1', 'Modulateur Q/2');
relier(mdl, 'Modulateur I/1', 'Signal emis/1'); relier(mdl, 'Modulateur Q/1', 'Signal emis/2');
relier(mdl, 'Signal emis/1', 'Canal/1'); relier(mdl, 'Bruit/1', 'Canal/2'); relier(mdl, 'Signal emis/1', 's/1');
relier(mdl, sI, 'bb/1'); relier(mdl, 'Signal emis/1', 'Mux/1'); relier(mdl, sI, 'Mux/2'); relier(mdl, 'Mux/1', 'Scope/1');
% récepteur en quadrature, filtre adapté
r0 = 5;
ajouter(mdl, [B 'Sources/Sine Wave'], 'OL cos', c0, r0, 'Amplitude', '2', 'Frequency', '2*pi*fc', 'Phase', 'pi/2', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Sources/Sine Wave'], 'OL -sin', c0, r0 + 2, 'Amplitude', '-2', 'Frequency', '2*pi*fc', 'Phase', '0', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Product'], 'Melangeur I', c0 + 1, r0, 'Inputs', '2');
ajouter(mdl, [B 'Math Operations/Product'], 'Melangeur Q', c0 + 1, r0 + 2, 'Inputs', '2');
ajouter(mdl, [B 'Discrete/Discrete Filter'], 'Filtre adapte I', c0 + 2, r0, 'Numerator', 'h_rx', 'Denominator', '1');
ajouter(mdl, [B 'Discrete/Discrete Filter'], 'Filtre adapte Q', c0 + 2, r0 + 2, 'Numerator', 'h_rx', 'Denominator', '1');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'Re-Im vers complexe', c0 + 3, r0 + 1, 'Input', 'Real and imag');
ajouter(mdl, [B 'Sinks/To Workspace'], 'z', c0 + 4, r0 + 1, 'VariableName', 'z', 'SaveFormat', 'Array');
relier(mdl, 'Canal/1', 'Melangeur I/1'); relier(mdl, 'OL cos/1', 'Melangeur I/2');
relier(mdl, 'Canal/1', 'Melangeur Q/1'); relier(mdl, 'OL -sin/1', 'Melangeur Q/2');
relier(mdl, 'Melangeur I/1', 'Filtre adapte I/1'); relier(mdl, 'Melangeur Q/1', 'Filtre adapte Q/1');
relier(mdl, 'Filtre adapte I/1', 'Re-Im vers complexe/1'); relier(mdl, 'Filtre adapte Q/1', 'Re-Im vers complexe/2');
relier(mdl, 'Re-Im vers complexe/1', 'z/1');
end

function mdl = construire_fsk(m, P)
mdl = nouveau(m.nom, P, m); hws = get_param(mdl, 'ModelWorkspace');
hws.assignin('df', P.df); hws.assignin('phi', 0.7); hws.assignin('h_rx', ones(1, P.sps) / P.sps);
B = 'simulink/';
ajouter(mdl, [B 'Sources/Uniform Random Number'], 'Source de bits', 0, 1, 'Minimum', '0', 'Maximum', 'M', 'Seed', 'graine', 'SampleTime', 'Ts');
ajouter(mdl, [B 'Math Operations/Rounding Function'], 'Indice', 1, 1, 'Operator', 'floor');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'NRZ', 2, 1, 'NumberOfTableDimensions', '1', 'Table', 'ptsI');
ajouter(mdl, [B 'Signal Attributes/Rate Transition'], 'Maintien', 3, 1, 'OutPortSampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Gain'], '2 pi df', 4, 1, 'Gain', '2*pi*df');
ajouter(mdl, [B 'Discrete/Discrete-Time Integrator'], 'Integrateur', 5, 1, 'IntegratorMethod', 'Integration: Forward Euler', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Sources/Clock'], 'Horloge', 4, 2.4);
ajouter(mdl, [B 'Math Operations/Gain'], '2 pi fc', 5, 2.4, 'Gain', '2*pi*fc');
ajouter(mdl, [B 'Sources/Constant'], 'phi', 5, 3.4, 'Value', 'phi');
ajouter(mdl, [B 'Math Operations/Sum'], 'Phase', 6, 2, 'Inputs', '+++');
ajouter(mdl, [B 'Math Operations/Trigonometric Function'], 'cos', 7, 2, 'Operator', 'cos');
ajouter(mdl, [B 'Sources/Random Number'], 'Bruit', 7, 3.4, 'Mean', '0', 'Variance', 'sigma2', 'Seed', '1', 'SampleTime', '1/fs');
ajouter(mdl, [B 'Math Operations/Sum'], 'Canal', 8, 2.6, 'Inputs', '++');
ajouter(mdl, [B 'Sinks/To Workspace'], 's', 8, 1.2, 'VariableName', 's', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'bb', 4, 0, 'VariableName', 'bb', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'sym', 2, 2.4, 'VariableName', 'sym', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Signal Routing/Mux'], 'Mux', 9, 1.2, 'Inputs', '2');
ajouter(mdl, [B 'Sinks/Scope'], 'Scope', 10, 1.2, 'NumInputPorts', '1');
relier(mdl, 'Source de bits/1', 'Indice/1'); relier(mdl, 'Indice/1', 'NRZ/1'); relier(mdl, 'Indice/1', 'sym/1');
relier(mdl, 'NRZ/1', 'Maintien/1'); relier(mdl, 'Maintien/1', '2 pi df/1'); relier(mdl, '2 pi df/1', 'Integrateur/1');
relier(mdl, 'Maintien/1', 'bb/1');
relier(mdl, 'Horloge/1', '2 pi fc/1'); relier(mdl, 'Integrateur/1', 'Phase/1'); relier(mdl, '2 pi fc/1', 'Phase/2');
relier(mdl, 'phi/1', 'Phase/3'); relier(mdl, 'Phase/1', 'cos/1'); relier(mdl, 'cos/1', 'Canal/1'); relier(mdl, 'Bruit/1', 'Canal/2');
relier(mdl, 'cos/1', 's/1'); relier(mdl, 'cos/1', 'Mux/1'); relier(mdl, 'Maintien/1', 'Mux/2'); relier(mdl, 'Mux/1', 'Scope/1');
% deux corrélateurs non cohérents, à fc + df (bit 1) et fc - df (bit 0)
noms = {'1', '+'; '0', '-'};
for t = 1:2
    n = noms{t, 1}; sg = noms{t, 2}; r0 = 5 + 3.2 * (t - 1);
    ajouter(mdl, [B 'Sources/Sine Wave'], ['OL cos ' n], 7, r0, 'Amplitude', '2', 'Frequency', ['2*pi*(fc' sg 'df)'], 'Phase', 'pi/2', 'SampleTime', '1/fs');
    ajouter(mdl, [B 'Sources/Sine Wave'], ['OL -sin ' n], 7, r0 + 2, 'Amplitude', '-2', 'Frequency', ['2*pi*(fc' sg 'df)'], 'Phase', '0', 'SampleTime', '1/fs');
    ajouter(mdl, [B 'Math Operations/Product'], ['Melangeur I ' n], 8, r0, 'Inputs', '2');
    ajouter(mdl, [B 'Math Operations/Product'], ['Melangeur Q ' n], 8, r0 + 2, 'Inputs', '2');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], ['Correlateur I ' n], 9, r0, 'Numerator', 'h_rx', 'Denominator', '1');
    ajouter(mdl, [B 'Discrete/Discrete Filter'], ['Correlateur Q ' n], 9, r0 + 2, 'Numerator', 'h_rx', 'Denominator', '1');
    ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], ['Re-Im ' n], 10, r0 + 1, 'Input', 'Real and imag');
    ajouter(mdl, [B 'Math Operations/Complex to Magnitude-Angle'], ['Module ' n], 11, r0 + 1, 'Output', 'Magnitude');
    ajouter(mdl, [B 'Sinks/To Workspace'], ['y' n], 12, r0 + 2, 'VariableName', ['y' n], 'SaveFormat', 'Array');
    relier(mdl, 'Canal/1', ['Melangeur I ' n '/1']); relier(mdl, ['OL cos ' n '/1'], ['Melangeur I ' n '/2']);
    relier(mdl, 'Canal/1', ['Melangeur Q ' n '/1']); relier(mdl, ['OL -sin ' n '/1'], ['Melangeur Q ' n '/2']);
    relier(mdl, ['Melangeur I ' n '/1'], ['Correlateur I ' n '/1']); relier(mdl, ['Melangeur Q ' n '/1'], ['Correlateur Q ' n '/1']);
    relier(mdl, ['Correlateur I ' n '/1'], ['Re-Im ' n '/1']); relier(mdl, ['Correlateur Q ' n '/1'], ['Re-Im ' n '/2']);
    relier(mdl, ['Re-Im ' n '/1'], ['Module ' n '/1']); relier(mdl, ['Module ' n '/1'], ['y' n '/1']);
end
ajouter(mdl, [B 'Math Operations/Sum'], 'Difference', 12, 6.6, 'Inputs', '+-');
ajouter(mdl, [B 'Sinks/To Workspace'], 'z', 13, 6.6, 'VariableName', 'z', 'SaveFormat', 'Array');
relier(mdl, 'Module 1/1', 'Difference/1'); relier(mdl, 'Module 0/1', 'Difference/2'); relier(mdl, 'Difference/1', 'z/1');
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

function [x, n] = echantillons(so, m, P, cal)
z = so.get('z'); z = z(:).'; sym = so.get('sym'); sym = round(sym(:).');
x = z(cal.off + 1:P.sps:end); x = x(1:min(numel(x), numel(sym) + 40));
n = numel(sym);
end

function [d] = decider(x, m, cal)
if strcmp(m.detection, 'enveloppe'), d = double(abs(x) > cal.gain / 2);
elseif strcmp(m.detection, 'energie'), d = double(real(x) > 0);
else
    [~, i] = min(abs(x(:) / cal.gain - m.pts(:).'), [], 2); d = i(:).' - 1;
end
end

function [es, ns, eb, nb] = compter(m, sym, d, L)
n = min(numel(sym), numel(d) - L) - 20;
a = sym(21:n); b = d(21 + L:n + L);
es = sum(a ~= b); ns = numel(a);
x = bitxor(m.lab(a + 1), m.lab(b + 1)); eb = 0;
for v = x, eb = eb + sum(dec2bin(v) == '1'); end
nb = ns * m.bits;
end

function cal = calibrer(so, m, P)
% retard du récepteur : phase d'échantillonnage et retard en symboles qui alignent le mieux
% les échantillons du filtre adapté sur les symboles émis (exécution sans bruit)
z = so.get('z'); z = z(:).'; sym = round(so.get('sym')); sym = sym(:).';
ref = m.pts(sym + 1); best = inf;
for off = 0:P.sps - 1
    x = z(off + 1:P.sps:end);
    for L = 0:20
        n = min(numel(ref), numel(x) - L) - 20; a = ref(21:n); b = x(21 + L:n + L);
        if strcmp(m.detection, 'energie'), b = real(b); end
        g = real(b * a') / real(a * a'); e = mean(abs(b - g * a) .^ 2) / mean(abs(g * a) .^ 2);
        if strcmp(m.detection, 'enveloppe'), g = real(abs(b) * a') / real(a * a'); e = mean(abs(abs(b) - g * a) .^ 2) / mean(abs(g * a) .^ 2); end
        if e < best, best = e; cal = struct('off', off, 'L', L, 'gain', g, 'evm', sqrt(e)); end
    end
end
end

function p = theorie(m, g)
switch m.nom
    case 'ook',  p = ber_ook_enveloppe(g);
    case 'ask4', p = 2 * 3 / 4 * Qf(0.5 * sqrt(2 * 2 * g / mean((0:3) .^ 2)));
    case {'bpsk', 'qpsk', 'psk8'}, M = m.M; gs = g * log2(M); a = sin(pi / M) ^ 2;
        p = integral(@(t) exp(-gs * a ./ sin(t) .^ 2), 0, (M - 1) * pi / M) / pi;
    case 'qam16', gs = 4 * g; q = 2 * (1 - 1 / 4) * Qf(sqrt(3 * gs / 15)); p = 1 - (1 - q) ^ 2;
    case 'fsk2', p = 0.5 * exp(-g / 2);
end
end

function q = Qf(x), q = 0.5 * erfc(x / sqrt(2)); end

function p = ber_ook_enveloppe(g)
% OOK, détection d'enveloppe, seuil A/2 : p = 1/2 [exp(-g/2) + 1 - Q1(2 sqrt(g), sqrt(g))]
a = 2 * sqrt(g); b = sqrt(g);
Q1 = integral(@(x) x .* exp(-(x .^ 2 + a ^ 2) / 2 + a * x) .* besseli(0, a * x, 1), b, b + a + 40);
p = 0.5 * (exp(-g / 2) + 1 - Q1);
end

function [f, Pn] = densite(s, fs)
N = 4096; w = 0.5 - 0.5 * cos(2 * pi * (0:N - 1) / N); P = zeros(1, N / 2 + 1); k = 0;
for i = 1:N / 2:numel(s) - N
    X = fft(s(i:i + N - 1) .* w); P = P + abs(X(1:N / 2 + 1)) .^ 2; k = k + 1;
end
P = P / k; Pn = P / max(P); f = (0:N / 2) * fs / N;
end

function [r, echecs] = verifier(mdl, m, P, ici)
echecs = {}; r = struct('titre', m.titre);
so = simuler(mdl, 4000, 0, P);
s = so.get('s'); s = s(:).'; Ps = mean(s(20 * P.sps:end) .^ 2); r.puissance_signal = Ps;
cal = calibrer(so, m, P); sym = round(so.get('sym')); sym = sym(:).';
[x, ~] = echantillons(so, m, P, cal); d = decider(x, m, cal); [es, ns] = compter(m, sym, d, cal.L);
r.sans_bruit = struct('phase_echantillonnage', cal.off, 'retard_symboles', cal.L, 'gain', cal.gain, 'evm', cal.evm, ...
    'erreurs_symboles', es, 'symboles', ns);
fprintf('Sans bruit : phase %d, retard %d symboles, gain %.4f, EVM %.4f, %d erreurs sur %d symboles\n', cal.off, cal.L, cal.gain, cal.evm, es, ns);
if es > 0, echecs{end + 1} = sprintf('%s : %d erreurs sans bruit', m.nom, es); end
if cal.evm > 0.05, echecs{end + 1} = sprintf('%s : EVM %.3f', m.nom, cal.evm); end
[f, Pn] = densite(s, P.fs); x0 = x(21 + cal.L:end - 20); bb = so.get('bb'); bb = bb(:).';
if strcmp(m.forme, 'rrc')
    B = (1 + P.alpha) * P.Rs / 2; frac = sum(Pn(f > P.fc - B & f < P.fc + B)) / sum(Pn);
    r.spectre = struct('part_puissance_bande', frac); fprintf('Spectre : %.4f de la puissance entre %g et %g Hz\n', frac, P.fc - B, P.fc + B);
    if frac < 0.98, echecs{end + 1} = sprintf('%s : %.3f de la puissance dans la bande', m.nom, frac); end
elseif strcmp(m.nom, 'fsk2')
    [~, i1] = max(Pn .* (f < P.fc)); [~, i2] = max(Pn .* (f > P.fc));
    r.spectre = struct('raies_Hz', [f(i1), f(i2)]); fprintf('Spectre : raies a %.1f et %.1f Hz\n', f(i1), f(i2));
    if abs(f(i1) - (P.fc - P.df)) > 20 || abs(f(i2) - (P.fc + P.df)) > 20, echecs{end + 1} = sprintf('%s : raies %g %g', m.nom, f(i1), f(i2)); end
else
    lv = [max(Pn(abs(f - (P.fc - P.Rs)) < 20)), max(Pn(abs(f - (P.fc + P.Rs)) < 20))]; lv = 10 * log10(lv);
    r.spectre = struct('niveau_premiers_zeros_dB', lv); fprintf('Spectre : premiers zeros a %.1f et %.1f dB\n', lv);
    if max(lv) > -20, echecs{end + 1} = sprintf('%s : zeros du spectre a %.1f dB', m.nom, max(lv)); end
end
taux = struct('EbN0_dB', {}, 'variance_bruit', {}, 'TES', {}, 'TEB', {}, 'erreurs', {}, 'mesure', {}, 'theorie', {}, 'conforme', {});
xmid = [];
for e = m.ebn0
    g = 10 ^ (e / 10); s2 = Ps * P.sps / (2 * m.bits * g);
    so = simuler(mdl, 16000, s2, P); sym = round(so.get('sym')); sym = sym(:).';
    [x, ~] = echantillons(so, m, P, cal); d = decider(x, m, cal); [es, ns, eb, nb] = compter(m, sym, d, cal.L);
    if strcmp(m.compare, 'TEB'), mes = eb / nb; ne = eb; else, mes = es / ns; ne = es; end
    th = theorie(m, g); ok = ne < 50 || (mes / th > 1 / 1.5 && mes / th < 1.5);
    taux(end + 1) = struct('EbN0_dB', e, 'variance_bruit', s2, 'TES', es / ns, 'TEB', eb / nb, 'erreurs', ne, ...
        'mesure', mes, 'theorie', th, 'conforme', ok); %#ok<AGROW>
    fprintf('Eb/N0 = %4.1f dB : %s mesure %.3e (%d erreurs), theorie %.3e  %s\n', e, m.compare, mes, ne, th, string(ok));
    if ~ok, echecs{end + 1} = sprintf('%s : Eb/N0 %g dB, mesure %.3g, theorie %.3g', m.nom, e, mes, th); end %#ok<AGROW>
    if e == m.ebn0(3), xmid = x / cal.gain; end
end
r.taux = taux;
so = simuler(mdl, 1500, Ps * P.sps / (2 * m.bits * 100), P); z = so.get('z'); z = z(:).' / cal.gain;
% figure : signal émis, spectre, constellation, oeil
fh = figure('Visible', 'off', 'Position', [0 0 1300 900]);
subplot(2, 2, 1); i0 = 40 * P.sps; t = (0:8 * P.sps - 1) / P.fs * 1e3;
plot(t, s(i0 + 1:i0 + 8 * P.sps)); hold on; plot(t, bb(i0 + 1:i0 + 8 * P.sps), '--', 'LineWidth', 1.4); grid on;
xlabel('temps (ms)'); title([m.titre ' : signal emis et voie I en bande de base']);
subplot(2, 2, 2); plot(f / 1e3, 10 * log10(Pn + 1e-16)); xlim([0 10]); ylim([-80 3]); grid on;
xlabel('frequence (kHz)'); ylabel('dB'); title('spectre du signal emis');
subplot(2, 2, 3);
if strcmp(m.nom, 'fsk2')
    y1 = so.get('y1'); y0 = so.get('y0'); plot(y1(cal.off + 1:P.sps:end) / cal.gain, y0(cal.off + 1:P.sps:end) / cal.gain, '.', 'MarkerSize', 4);
    xlabel('|y_1|'); ylabel('|y_0|'); title('points de decision (|y_1|, |y_0|), Eb/N0 = 20 dB');
else
    plot(real(xmid), imag(xmid), '.', 'MarkerSize', 3); hold on; plot(real(x0 / cal.gain), imag(x0 / cal.gain), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
    axis equal; title(sprintf('constellation, Eb/N0 = %g dB (points) et sans bruit (+)', m.ebn0(3)));
end
grid on;
subplot(2, 2, 4); hold on; zz = real(z); if strcmp(m.detection, 'enveloppe'), zz = abs(z); end
for j = 60:260
    a = j * P.sps + cal.off - P.sps / 2; plot((0:2 * P.sps) / P.sps - 0.5, zz(a + 1:a + 2 * P.sps + 1), 'Color', [0 0.45 0.74]);
end
grid on; xlabel('temps (symboles)'); title('diagramme de l''oeil, voie I, Eb/N0 = 20 dB');
sgtitle(['mod\_' m.nom '  -  MATLAB R' version('-release')]);
exportgraphics(fh, fullfile(ici, ['mod_' m.nom '_resultat.png']), 'Resolution', 130); close(fh);
end

function figure_taux(mods, R, ici)
fh = figure('Visible', 'off', 'Position', [0 0 1000 700]); co = lines(numel(mods));
for i = 1:numel(mods)
    m = mods(i); t = R.(m.nom).taux; e = [t.EbN0_dB]; y = [t.mesure]; y(y == 0) = NaN;
    ee = linspace(min(e) - 1, max(e) + 1, 50); th = arrayfun(@(x) theorie(m, 10 ^ (x / 10)), ee);
    semilogy(ee, th, '-', 'Color', co(i, :)); hold on;
    semilogy(e, y, 'o', 'Color', co(i, :), 'DisplayName', [m.titre ' (' m.compare ')']);
end
ylim([1e-5 0.5]); grid on; xlabel('Eb/N0 (dB)'); ylabel('taux d''erreur');
title(['Taux d''erreur mesures dans Simulink (points) et theorie (traits) - MATLAB R' version('-release')]);
h = findobj(gca, 'Marker', 'o'); legend(flipud(h), 'Location', 'southwest');
exportgraphics(fh, fullfile(ici, 'taux_erreur_simulink.png'), 'Resolution', 130); close(fh);
end
