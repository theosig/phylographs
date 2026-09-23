# PhyloGraphs

Code accompanying:

> Signor, T., Neto, M., Jofré, P., Vitali, S., Hua, X., Tissera, P. B., Aguilera-Gómez, C., Das, P., Tapia-Contreras, B., Yates, R. M., Martí, L., & Sánchez-Pi, N. (2026).
> **Comparing phylogenetic trees of stars with spectral graph distances.**
> *Astronomy & Astrophysics*, accepted. <!-- TODO: add DOI/arXiv link once assigned -->

We compare stellar phylogenetic trees, built with Neighbor-Joining from chemical abundances, using classical tree-distance measures and a set of graph-spectral distances (adjacency, Laplacian, and normalized Laplacian spectral distance) that don't require the trees to share leaves. This lets trees built from different, disjoint stellar samples be compared directly — something the classical measures can't do.

## Repository contents

- `phylograph.py` — core module: pairwise chemical distances, Neighbor-Joining tree construction, tree-to-graph conversion, the three spectral distance measures (ASD, LSD, NLSD), and the `RunConfig`/benchmark-running machinery used throughout.
- `abundances.py` — converts [X/Fe] abundance columns to [X/H] (Appendix B robustness check on the choice of reference element).
- `run_histories.ipynb` — data cleaning and population-matching for the four benchmark pairs (MW–GSE, MW–Sgr, GSE–Sgr, Nissen older/younger); reproduces Table 2 and Figs. 4–5, and the Appendix B sensitivity analysis (tau/R(w) calibration, standardization, leave-one-element-out) — Table B.1, Table B.2, and Figs. B.1–B.2.

## Dependencies



## Reproducing the results

1. Obtain the input catalogs (see **Data** below) and place them under `data/` following the paths referenced in `run_histories.ipynb`.
2. Run `run_histories.ipynb` top to bottom. It generates the trees, computes ASD/LSD/NLSD distances, and reproduces Table 2, Figs. 4–5, Table B.1, Table B.2, and Figs. B.1–B.2.

## Data

This repository does not redistribute the underlying survey catalogs. The abundance samples used are:

- **Solar twins:** Nissen et al. (2020, A&A 640, A81).
- **Milky Way & Gaia–Sausage–Enceladus:** chemo-kinematical catalog of Re Fiorentin et al. (2024, ApJ 977, 278).
- **Sagittarius dSph:** candidate members from Vitali et al. (2022, MNRAS 517, 6121), cross-matched with APOGEE DR17 (Abdurro'uf et al. 2022, ApJS 259, 35).

See Sect. 3 of the paper for the exact selection cuts (`X_FE_FLAG`, `STARFLAG`, metallicity matching) applied to each sample.


## Citation

If you use this code, please cite the paper above.
<!-- TODO: add a BibTeX entry once the DOI is assigned -->

## Contact

Theosamuele Signor — theo.signor@gmail.com
