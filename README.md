# VEGF-pathway resistance architecture: reproducibility code and source data

This archive contains the analysis code, frozen figure source data, selected result tables, sample-level manifests, environment records, and reference figures for the manuscript *Recurrent and context-dependent transcriptional architectures of resistance to VEGF-pathway inhibition across experimental and human contexts*.

The analysis preserves study-specific contrasts and compartments. Discovery contains five datasets and eight resistance contexts. Primary experimental validation contains four contrasts from two datasets; the three GSE64052 contrasts are nested within one dataset. GSE81465 is supportive-only because replicate independence is unresolved. Human analyses are kept separate from experimental validation.

## Package contents

- `scripts/`: ordered R analysis and figure-rendering scripts, plus the archive verifier.
- `config/`: final dataset and contrast evidence roles.
- `metadata/`: dataset sources and sample-level inclusion records.
- `data/figure_source_data/`: frozen inputs used to render the publication figures.
- `results/`: selected manuscript-facing result tables.
- `reference_figures/`: reference PNG, PDF, and SVG files.
- `environment/`: software versions and the recorded R session.
- `checksums/`: SHA-256 manifest for archive integrity.

## Verify the archive

From the unpacked archive directory, run:

```bash
python3 scripts/verify_release.py
```

A successful check ends with `RELEASE_INTEGRITY=PASS`.

## Re-render the figures

Figure rendering is self-contained because all required source-data tables are included:

```bash
Rscript scripts/11_render_publication_figures.R
```

New files are written to `outputs/figures/`. To render selected figures only, provide a comma-separated list:

```bash
VEGF_FIGURES=Figure1,Figure4 Rscript scripts/11_render_publication_figures.R
```

The renderer uses the included frozen TSV files and does not rerun differential expression or enrichment analysis.

## Re-run analyses from expression data

The expression datasets are publicly available from NCBI GEO and are not duplicated in this archive. Dataset links and required input types are listed in `metadata/public_data_sources.tsv`; the expected directory layout is documented in `DATA_ACCESS.md`.

After reconstructing that layout, point the scripts to its root:

```bash
export VEGF_ANALYSIS_ROOT=/absolute/path/to/analysis_directory
Rscript scripts/01_discovery_preprocess.R
Rscript scripts/02_discovery_differential_expression.R
Rscript scripts/03_discovery_gsea.R
Rscript scripts/04_cross_context_architecture.R
```

Continue in the order given in `run_manifest.tsv`. The GSE249415 sensitivity script instead takes its input through `VEGF_GSE249415_INPUT` and can write to a user-selected `VEGF_OUTPUT_ROOT`.

## Reproducibility boundary

The figure-rendering path is directly runnable from this archive. Full expression-level reconstruction additionally requires downloading public GEO files and arranging them according to `DATA_ACCESS.md`. Randomized analyses use explicit seeds. Package versions are recorded in `environment/software_versions.tsv`; the complete recorded R session is in `environment/R_sessionInfo.txt`.

