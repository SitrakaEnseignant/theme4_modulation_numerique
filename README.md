# theme4_modulation_numerique

Thème 4 : modulations numériques OOK, ASK à 4 niveaux, BPSK, QPSK, 8-PSK, 16-QAM et FSK binaire.
Chaque modulation dispose d'un modèle MATLAB Simulink (blocs de la bibliothèque de base) et d'un graphe de flux
GNU Radio 3.10, en simulation pure : aucun matériel de radio logicielle n'est utilisé, le canal est un bruit blanc gaussien.

Pour chaque simulateur, l'exécution réelle sur GitHub Actions produit le signal émis, son spectre, la constellation,
le diagramme de l'œil et le taux d'erreur en fonction de Eb/N0, comparé à la théorie.

| Paramètre | Valeur |
|---|---|
| Cadence d'échantillonnage | 32 kHz |
| Rapidité de modulation | 1 000 symboles par seconde (32 échantillons par symbole) |
| Porteuse | 4 kHz |
| Mise en forme | rectangulaire (OOK, ASK, FSK), racine de cosinus surélevé, retombée 0,35 (PSK, QAM) |
| FSK | tons à 3,5 et 4,5 kHz (indice 1), détection non cohérente |

## Contenu

- `Simulink/construire_et_verifier_modeles.m` : construit les sept modèles `mod_*.slx`, les exécute et vérifie les résultats.
- `GNU_Radio/*_tx_rx_sim.grc` : graphes de flux à ouvrir dans GNU Radio Companion 3.10 ; générés par `outils/generer_grc.py`.
- `En_ligne/executer_grc_ci.py` : exécution des graphes sans écran, mesures et figures (utilisé par GitHub Actions).
- `En_ligne/modulations.py`, `En_ligne/theorie.py` : paramètres communs et formules théoriques.
- `.github/workflows/simulations.yml` : exécution automatique ; résultats dans la branche `resultats`.
