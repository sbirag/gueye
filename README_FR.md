# gueye 0.7.1 — Laplace / Helmholtz discrétisé, 1D-3D, toutes conditions aux limites

Résout  **−Δu + μ²u − k²u = f**  sur une boîte 1D, 2D ou 3D (différences finies d'ordre 2), avec des conditions
**Dirichlet / Neumann / Robin / transparente / PML** indépendantes sur chaque face, et l'équation d'onde en temps.
Méthode : inversion semi-analytique (formule fermée tridiagonale, axes transverses diagonalisés) — pas de théorème à
connaître pour s'en servir.

## Installation
    pip install -e .            # NumPy et SciPy suffisent ; matplotlib pour les exemples
    python -m gueye             # auto-test (< 1 s) : doit finir par « TOUT EST OK »

## Interface graphique
    python -m gueye.gui        # ouvre le navigateur ; ou  gueye-gui ;  ou  import gueye; gueye.launch_gui()
Curseur k en temps reel, conditions aux limites face par face, clic pour placer la source (Maj+clic : recepteur), obstacle circulaire 2D,
coupes 3D, balayage en frequence avec resonances exactes, animation d'onde. Le serveur n'ecoute que 127.0.0.1.

## Cinq minutes
    import numpy as np, gueye as gy

    # boîte 3 x 2, grille 300 x 200, k = 9 ; axe x : (gauche Dirichlet, droite Robin alpha = 3) ; axe y : (bas Neumann, haut Dirichlet)
    pb = gy.problem(size=(3.0, 2.0), cells=(300, 200), k=9.0, bc=[("D", ("R", 3.0)), ("N", "D")])

    u = pb.solve(lambda x, y: np.exp(-((x - 1)**2 + (y - 1)**2) / 0.01))     # champ complet, tableau (y, x)
    g = pb.green((0.9, 0.7))                                                 # réponse à une source ponctuelle
    pb.resonances(6)                                                         # k_n exacts du schéma
    G = pb.sweep(np.linspace(1, 12, 500), src=(0.9, 0.7), rec=[(2.0, 1.2)])  # réponse en fréquence
    pb.at([(2.0, 1.2)], f)                                                   # valeur en quelques points, sans le champ entier

    W = gy.wave(size=(3.0, 2.0), cells=(120, 80), bc=[("D", ("R", 1.5)), ("N", "D")], c=1.0, eta=0.02)
    t, U, R = W.march(0.02, 400, source=W.point_source((0.9, 0.7), signal), receivers=[(1.0, 1.0)])

## Conditions aux limites (une paire (bas, haut) par axe x, y, z ; normale sortante)
| écriture | sens |
|---|---|
| `"D"`, `("D", g)` | u = g (g : réel ou fonction des coordonnées de la face) |
| `"N"`, `("N", g)` | ∂u/∂n = g |
| `("R", alpha)`, `("R", alpha, g)` | ∂u/∂n + alpha·u = g ; alpha complexe permis (`1j*k` : absorbant en incidence normale) |
| `"A"` | transparente exacte du schéma (demi-espace discret), **dernier axe** seulement |
| `"P"`, `("P", m, R0)` | PML de m mailles, réflexion visée R0 (défaut 16 mailles, 1e-6), **sur tous les axes** (domaine ouvert sur toutes les faces) |
`k` complexe = absorption ; `mu` = masse (milieu écranté). Les 9 combinaisons D/N/R par axe sont couvertes.

## Recettes (`gueye.recipes`)
`resonance_table`, `frequency_response` (pics vs résonances), `open_field` (milieu ouvert), `dirichlet_obstacle`
(obstacle conducteur 2D, capacité embarquée ; aussi en domaine ouvert sur 4 faces), `grid_dispersion` (mailles par longueur d'onde). Exemples complets avec figures : `examples/`.

## Garde-fous (comportement garanti, testé)
- **Résonance / singularité** : `problem(..., check=True)` (défaut) lève `SingularOperator` si le déterminant relatif de l'opérateur
  est < 1e-10 (le résultat n'aurait plus de sens) et avertit au-delà de 1e-6 (perte de chiffres). `check=False` désactive tout
  contrôle, à vos risques ; `check=1e-12` fixe un seuil. `sweep` met **NaN** aux fréquences singulières et regroupe les avertissements.
- **Points** : un point hors du domaine (au-delà d'une demi-maille) → `ValueError`. Une **source** sur un nœud de Dirichlet → `ValueError`
  (placer la source à une maille du bord). Un **récepteur** sur une face de Dirichlet → 0 pour `value` (Green homogène) et valeur
  imposée pour `at`. Les récepteurs du temporel lisent le champ complet (bords inclus).
- **Entrées** : `size` > 0 et fini, `cells` entiers ≥ 2, une paire (bas, haut) par axe, PML avec m ≥ 2 et 0 < R0 < 1 → sinon `ValueError`.
- **`resonances(n)`** : tableau **réel** si le spectre est réel positif ; sinon (Robin complexe ou α < 0, opérateur non auto-adjoint)
  des k **complexes** avec un avertissement — jamais de fausses résonances réelles.

## Quand l'utiliser (honnêtement)
- **Oui** : conditions mixtes par face ; balayage en fréquence (une fréquence = un élément de l'inverse) ; fonctions de Green ;
  obstacles/conducteurs par capacité ; PML sur tous les axes sans refactoriser.
- **Pas forcément** : Poisson tout-Dirichlet sur un cube où l'on veut le champ complet — la DST (FFT) est 3 à 20 fois plus rapide
  (mesuré, 3D, N = 24…96). **`at` n'est pas une accélération universelle** : pour M = 10⁶ en 1D, 3 points par la formule ≈ 1,9 s contre
  ≈ 0,1 s pour la résolution bande du champ complet ; son intérêt est la mémoire quand le champ entier est trop gros.
- Domaines non rectangulaires : l'obstacle passe par la capacité embarquée, pas par un maillage adapté.
- Règle de maille : au moins ~10 mailles par longueur d'onde ; l'erreur de dispersion de grille vaut ≈ (kh)²/24 (`grid_dispersion`).
- Précision : les écarts ≤ 1e-12 annoncés valent pour les cas bien conditionnés ; près d'une résonance (déterminant relatif ε) on
  perd ≈ log10(1/ε) chiffres (cas le moins précis du banc de noyau : 2,6e-7 pour d = 2 + 1e-9, Neumann pur).

## Temporel : CL dépendant du temps et domaine ouvert (0.6.0)
    W = gy.wave(size=(1, 1), cells=(80, 80), bc=[("D", "D"), ("N", "N")], bc_time={"left": lambda t: np.sin(6*t)})   # g(x,t) = g0(x) s(t) par face
    from gueye.recipes import time_response       # domaine OUVERT en temps (PML/transparente) par synthèse de Fourier
    U = time_response(pb, src, [rec], signal, t, wmax=14.0, nw=140, shift=0.5)       # (nrec, len(t))
Validé : fermé vs Newmark (2e-4), ouvert PML 2D vs boîte 4x plus grande avant retour des échos (1,5e-3 relatif).

## Limites connues
- Transparente exacte (`"A"`) : dernier axe seulement. La PML est admise sur tous les axes ; sur les axes diagonalisés (x, y) elle est
  rendue explicite (couches de nœuds) et sa diagonalisation dépend de k : `with_k`/`sweep` refont une diagonalisation complexe par
  fréquence (O((Nx+2m)³) en 2D) — nettement plus lent qu'un balayage sans PML sur x.
- Pas de PML dans `wave(...).march` (Newmark pas à pas) : l'ouvert en temps passe par `recipes.time_response` (une résolution par fréquence ; plus coûteux qu'un pas de temps, pas de PML « ADE » en pas à pas).
- Termes d'ordre 1 (c1) : 1D seulement ; en 2D/3D, coefficients anisotropes cx, cy, cz admis. Obstacle embarqué : 2D.
- Performance : le solveur modal vectorisé (0.6.0) donne ×1,8 à ×5 sur le 3D, mais **pas de gain en 2D** (mesuré) ; la DST reste préférable pour un cube tout-Dirichlet.
- Qualité de la PML : l'élimination algébrique est exacte (écart ~1e-14 à l'assemblage explicite), mais la réflexion résiduelle dépend
  de m, R0, de l'angle et de la discrétisation (10⁻⁵ à 10⁻² dans les tests) : ce n'est pas un milieu ouvert exact.
- `Problem1D` accepte un terme du premier ordre (`c1`) ; 2D/3D : opérateur de Laplace seul (c1 = 0), coefficients isotropes pour la PML.
- `dirichlet_obstacle` : ~6-10 s pour ~1000-1400 nœuds d'obstacle sur 200×150 à 240×180 (coût dominé par le bloc de Green).
- Temporel : une factorisation LU tridiagonale est calculée une fois par mode et réutilisée à chaque pas (dt constant).

## Validation (dossier `validation/`, journal `validation/validation_all.txt`)
Noyau vs inversion dense (89 cas, ≤ 1e-12 hors cas quasi singuliers) ; 9 CL 1D, 9 combinaisons 2D, 8 combinaisons 3D vs
assemblage creux (≤ 1e-13) ; ordre de convergence 2,00 avec CL inhomogènes ; Newmark vs LU creux (≤ 1e-14) ; PML vs
transparente exacte et vs grand domaine, PML sur toutes les faces vs assemblage explicite (2D, 3D) et vs domaine 3x plus grand ;
`tests/` : `pytest` (test_core.py = validations essentielles rapides, test_api.py = API, garde-fous, régressions de l'audit).


## Project and release information (English)

- **Release:** 0.7.1
- **Author:** Serigne Bira Gueye, Department of Physics, Faculty of Sciences and Techniques, Cheikh Anta Diop University, Dakar, Senegal.
- **Contact:** bira.gueye@ucad.edu.sn
- **License:** MIT (see `LICENSE`).
- **Repository:** https://github.com/sbirag/gueye
- **Citation metadata:** see `CITATION.cff`.

The optional browser interface is a local exploratory and visualization tool; it is not a substitute for the documented numerical validation. The server binds to loopback and is not intended for network exposure.

### AI-assisted development disclosure

ChatGPT (OpenAI) and Claude (Anthropic) were used as aids for drafting, code review and verification, troubleshooting, and technical advice during software development and manuscript preparation. The maintainer remains responsible for the software, the execution and interpretation of tests, and all scientific claims.
