# HumanCD4CoDEGNet — transcriptional buffering of essential genes

**Built with Claude: Life Sciences** (Researcher Track), 7–13 July 2026.
Finding: the yeast transcriptional-buffering law generalises to primary human CD4+ T cells —
but is masked, and sign-inverted, by an expression–power confound.

Method lineage: **YeastCoDEGNet** (Nasar, **Rehman**, Ott & Alam, *NAR* 2026, gkaf1410).
Data: Marson-lab genome-scale CRISPRi Perturb-seq (Zhu, Dann et al. 2025).

## Result in one line

| view | essential vs non | reading |
|---|---|---|
| raw responsiveness | Cliff's δ = **+0.26**, p=8×10⁻⁷ | essential look MORE responsive (confound) |
| baseMean-matched | δ = **−0.17**, p=8×10⁻⁶ (9/10 deciles) | essential are BUFFERED |
| CEGv2 essentials | δ = **−0.19**, p=1×10⁻¹⁴ | second list, same |
| shet constraint | ρ = **−0.09**, p=3×10⁻²⁰ | continuous, same |

Unanimous across 2 essential-gene sets × pooled + 3 stimulation states (12 strata). Strongest at Stim8hr.

## Run it

```bash
cd ~/CoDEG_Tcell
python3.12 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
jupyter lab notebooks/HumanCD4CoDEGNet_buffering.ipynb   # reproduces in ~2 min from one public file
```

The notebook streams a single public S3 file by byte-range (never downloads the 16.8 GB whole),
fetches gene lists from raw GitHub, and regenerates every number and figure. No credentials.

## Layout

```
notebooks/
  HumanCD4CoDEGNet_buffering.ipynb   THE submission notebook (executed, figures inline, ~2 min)
analyze_all.py        per-condition × multi-axis test table  -> artifacts/results_table.csv
buffering_test.py     the core headline test (Hart + shet)
robustness_matched.py non-parametric baseMean-matched validation
figures.py            fig1 sign-flip, fig2 robustness forest, fig3 shet dose-response
run_day1.py           independent binary-matrix build (cross-checks n_regulators to 3 decimals)
artifacts/
  figures/*.png                      the three figures
  results_table.csv                  all 12 strata
  per_gene_full.csv                  per-gene responsiveness + essentiality + baseMean + shet
  buffering_test_results.json, robustness_matched.json, day1_summary.json
SUBMISSION.md         abstract, figure captions, methods, reproducibility, limitations
```

## What we verified (not assumed)

- Dataset shapes/layers/flags read from the file, not the docs.
- QC filter → **17,260** usable perturbation-conditions (predicted from flags, confirmed).
- Density **0.74%** matches the authors' own `n_downstream` counts (independent cross-check).
- My independent responsiveness build agrees with the authors' `n_regulators` to 3 decimals.
- Read the preprint + analysis repo: the authors compute responsiveness and its power-correction
  (Supp Fig 6) but **never join it to essentiality** (their Hart list is loaded then commented out).
  That join — plus the yeast comparator — is the contribution.

## Four traps handled

1. Yeast's FC>2 cutoff does not transfer (median |log2FC|=0.095) → threshold on FDR alone.
2. `log_fc` is log2, not natural log.
3. `DE_stats.suppl_table.csv` is a truncated copy of `.obs` missing 4 QC flags → read flags from h5ad.
4. CRISPRi on-target self-DE masked (82.8% of rows).

## Prior work

Nasar MI, Rehman SSU, Ott S, Alam MT. *Uncovering coordinated pathway interactions through gene
co-differential expression in yeast.* Nucleic Acids Research 54(1), 2026, gkaf1410.
Yeast network: NDEx `794dd2cc-bcae-11f0-a218-005056ae3c32`.
