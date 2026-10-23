# gueye 0.7.1 release notes

This release candidate is based on the submitted `gueye 0.7.0` source archive and includes packaging/documentation updates for a reproducible public release:

- Adds MIT license files (`LICENSE` and `Licence.txt`).
- Adds `CITATION.cff`, project URLs, author metadata, and a security policy.
- Provides English `README.md` and retains the French documentation in `README_FR.md`.
- Adds a GitHub Actions test workflow for Python 3.9–3.13.
- Tightens the local GUI server's Host/Origin checks and request-length handling.
- Adds regression tests for cross-origin and mismatched-port request rejection.

Validation in the preparation environment (Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0): 158 tests passed; `python -m gueye` ended with `TOUT EST OK`; wheel and source distribution built successfully. These checks do not constitute a cross-platform CI run or prove that the remote GitHub repository has been updated.
