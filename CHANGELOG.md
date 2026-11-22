## 0.7.2 (corrective release preparation)
- Corrects discrete resonance wavenumbers by including the screening term `mu**2`: `k_n**2 = lambda_n + mu**2`.
- Preserves complex prescribed Dirichlet data in 1D–3D boundary evaluation, right-hand-side assembly, and field expansion.
- Adds regression tests for shifted resonances and complex Dirichlet data.
- Adds reproducible benchmark scripts and separates modal setup/solve from sparse assembly/factorization/solve timings.
- Adds independent sparse-reference checks for Green columns, frequency sweeps and discrete resonances in 1D–3D, plus 2D Dirichlet obstacles on two geometries and two grids.
- Local release validation: 166 tests passed; see `RELEASE_NOTES_v0.7.2.md`.

## 0.7.1 (release-preparation update)
- Adds the MIT license file, citation metadata (`CITATION.cff`), project URLs, package author metadata, and a security policy.
- Hardens the local browser interface by checking the Host authority and, when present, the Origin authority against the loopback server and its exact port; rejects invalid/negative request lengths and unsupported transfer encoding.
- Adds regression tests for Host-port and cross-origin request rejection.

## 0.7.0
- **Interface graphique interactive** (`python -m gueye.gui`, `gueye-gui`, ou `gy.launch_gui()`) : page HTML locale servie par la bibliotheque standard (aucune dependance en plus), 1D/2D/3D.
  Champ en temps reel (curseur k, CL par face D/N/R/A/P, source placee a la souris, obstacle circulaire 2D, coupes 3D), balayage en frequence avec resonances exactes,
  animation d'onde (Newmark), exemples predefinis, theme clair/sombre. Garde-fous de gueye conserves (resonance, source sur un bord de Dirichlet, PML refusee en temps).
  Serveur lie a 127.0.0.1, refuse les en-tetes Host non locaux, plafonne la taille des grilles.

## 0.6.0
- **Solveur modal vectorisé** (`kernel.ModeStack`) : tous les modes d'un axe en un balayage O(N) vectorisé, formule fermée conservée ; identique à la LU à 1e-14 ; 3D ×1,8–×5, 2D sans gain (mesuré). `method="lu"` garde l'ancienne boucle.
- **CL dépendant du temps** : `bc_time={face: s(t)}` dans Wave2D/Wave3D (g(x,t)=g0(x)s(t)) ; face inconnue -> ValueError.
- **Temporel ouvert** : `recipes.time_response` (synthèse de Fourier à fréquence complexe, toutes CL dont PML/transparente, 1D/2D/3D).
- Intègre les corrections de la 0.5.1 (points strictement dans [0, L]).
- Tests : +7 (149 au total). Limites restantes documentées dans le README (c1 en 2D/3D, Newmark+PML pas à pas, obstacle 3D).

## 0.5.1
- **Coordonnees physiques strictes** : `Solver.node` et `WaveSolver.node` rejettent maintenant toute coordonnee hors de `[0, L]`, y compris si elle est a moins d'une demi-maille du bord. Les points valides dans la boite restent associes au noeud le plus proche.
- Tests de regression ajoutes pour les points `-0.01` et `L + 0.01` sur un maillage de pas `0.05`, pour les problemes stationnaires et temporels.
- Le reste de la version 0.5.0 est conserve.

# Changelog

## 0.5.0
Corrections issues de l'audit de la 0.4.0 :
- **Singularités** : `Solver`/`problem` ont `check=True` par défaut (SingularOperator si déterminant relatif < 1e-10, avertissement < 1e-6) ;
  `with_k` hérite ; `sweep` renvoie NaN aux fréquences singulières. Test de régression sur la résonance 1D D-D (k = 3,1407852607) et Neumann pur.
- **Bords de Dirichlet** : source sur noeud Dirichlet et point hors domaine -> ValueError (fin du `np.clip` silencieux) ; récepteur Dirichlet -> 0 (`value`)
  ou valeur imposée (`at`) ; récepteurs du temporel (Wave2D/3D) lus sur le champ complet (plus de lecture par indice négatif) ; mêmes règles dans `wave`.
- **`resonances`** : k complexes renvoyés (avec avertissement) quand le spectre n'est pas réel positif (Robin complexe, alpha < 0).
- **Entrées** : validation de `size`, `cells`, `bc`, PML (m, R0) -> ValueError explicites.
- **Performance** : factorisation LU tridiagonale (gttrf) mise en cache par noyau (réutilisée à chaque pas et par mode).
- **Tests** : `tests/test_core.py` (validations essentielles) et `tests/test_api.py` étendu.
Nouveauté :
- **PML sur toutes les faces** (2D : gauche/droite en plus de bas/haut ; 3D : x, y en plus de z) par couches explicites sur les axes diagonalisés ;
  validée à 1e-15 contre l'assemblage explicite, et contre un domaine 3x plus grand (2,7e-5). `recipes.open_field(all_open=True)`, exemple `ex5`.

## 0.4.0
API unifiée, recettes, installation, auto-test.
