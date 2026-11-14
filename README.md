# theme4_modulation_numerique

Bloc 2, Thème 4 : modulations numériques OOK, ASK à 4 niveaux, BPSK, QPSK, 8-PSK, 16-QAM et FSK binaire.
Chaque modulation dispose d'un modèle MATLAB Simulink et d'un graphe de flux GNU Radio 3.10 construits bloc pour bloc
de la même façon (mêmes noms de blocs, mêmes paramètres), en simulation pure : aucun matériel de radio logicielle
n'est utilisé, le canal est un bruit blanc gaussien.

Le récepteur n'est pas idéal : son oscillateur local est décalé de 0,7 rad en phase et de 1 Hz en fréquence.
La synchronisation de porteuse est donc assurée, comme dans un vrai récepteur, par une boucle de Costas
(BPSK, QPSK, 8-PSK) ou par un récepteur de constellation piloté par les décisions (16-QAM).

| Paramètre | Valeur |
|---|---|
| Cadence d'échantillonnage | 32 kHz |
| Rapidité de modulation | 1 000 symboles par seconde (32 échantillons par symbole) |
| Porteuse | 4 kHz |
| Mise en forme | rectangulaire (OOK, ASK, BPSK, QPSK, FSK), racine de cosinus surélevé de retombée 0,35 (8-PSK, 16-QAM) |
| OOK, ASK | détection d'enveloppe (module), seuils au milieu des niveaux |
| BPSK, QPSK, 8-PSK | boucle de Costas d'ordre 2, 4, 8 ; bande normalisée 2π/12800 |
| 16-QAM | récepteur de constellation (boucle pilotée par les décisions) ; bande normalisée 2π/400 |
| FSK | tons à 3 et 5 kHz (indice 2) ; deux corrélateurs non cohérents et un discriminateur de fréquence |
| Ambiguïté de phase | levée par 100 symboles connus (pilotes) au début de chaque bloc de 1 000 symboles |

## Contenu

- `Simulink/construire_et_verifier_modeles.m` : construit les sept modèles, les exécute et vérifie les résultats.
- `Simulink/modeles/` : modèles `mod_*.slx` et leurs schémas `mod_*_schema.png`, produits par l'exécution sur GitHub Actions.
- `Simulink/resultats/` : figures `mod_*_resultat.png` (signal émis, spectre, constellation, diagramme de l'œil),
  `taux_erreur_simulink.png`, `mesures_simulink.json`, `journal_matlab.txt`.
- `GNU_Radio/*_tx_rx_sim.grc` : graphes de flux à ouvrir dans GNU Radio Companion 3.10 ; générés par `outils/generer_grc.py`.
- `GNU_Radio/resultats/` : figures `*_tx_rx_sim.png`, `taux_erreur_gnuradio.png`, `mesures.json`, `journal_gnuradio.txt`.
- `En_ligne/executer_grc_ci.py` : exécution des graphes sans écran, mesures et figures.
- `En_ligne/modulations.py`, `En_ligne/theorie.py` : paramètres communs et formules théoriques.
- `.github/workflows/simulations.yml` : exécution automatique ; les résultats sont recopiés dans `main` et dans la branche `resultats`.

## Contrôles automatiques

Pour chaque modulation et dans chaque outil : aucune erreur sans bruit, constellation et gain corrects (EVM),
spectre conforme (zéros à fc ± k·Rs pour l'impulsion rectangulaire, 98 % de la puissance dans la bande pour le
cosinus surélevé, raies à fc ± 1 kHz pour la FSK), puis taux d'erreur mesuré comparé à la théorie pour cinq Eb/N0.
Aux faibles Eb/N0, une boucle peut glisser de 2π/ordre (glissement de cycle) : ces points sont tracés par une croix
et notés « hors accrochage » ; lorsque le taux d'erreur théorique est inférieur à 1 %, un glissement isolé est toléré et signalé, plusieurs glissements font échouer le contrôle.

## Références

Les graphes de flux suivent les leçons de modulation numérique de Learn SDR (dépôt `gallicchio/learnSDR`,
leçons OOK, ASK, BPSK, QPSK et FSK) : détection d'enveloppe pour OOK et ASK, impulsions rectangulaires et filtre à
moyenne glissante pour BPSK et QPSK, boucle de Costas, FSK d'indice 2 démodulée par discriminateur.
