# VEGF-pathway resistance architecture: reproducibility code and source data

This archive contains the analysis code, frozen figure source data, selected result tables, sample-level manifests, environment records, and reference figures for the manuscript *Recurrent and context-dependent transcriptional architectures of resistance to VEGF-pathway inhibition across experimental and human contexts*.

The analysis preserves study-specific contrasts and compartments. Discovery contains five datasets and eight resistance contexts. Final external validation contains 12 contexts from eight independent GEO datasets; contexts within GSE64052, GSE84048, and GSE78698 remain nested within their GEO series. Dataset-level recurrence is the primary summary and context-level results are secondary. GSE81465 and GSE45161 are supportive-only. Human analyses are kept separate from experimental validation.

## Package contents

- `scripts/`: ordered R/Python analysis, validation-summary, and figure-rendering scripts, plus the archive verifier. The `legacy_phase9A/` subdirectory preserves the superseded two-dataset validation implementation for provenance and is not part of the current execution path.
- `config/`: final dataset and contrast evidence roles.
- `metadata/`: dataset sources and sample-level inclusion records. `all_sample_inclusion_QC.tsv` is the authoritative cross-stage sample/QC record; `external_validation_sample_manifest.tsv` is its formal eight-dataset validation subset; and `legacy_phase9A_validation_sample_manifest.tsv` is retained only as an archival record of the superseded Phase 9A implementation.
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
Rscript scripts/14_render_phase10B_validation_figures.R
```

New files are written to `outputs/figures/`. The first command renders the eight figures unchanged by the expanded validation cohort; the second renders Figure 1, Figure 4, and Supplementary Figure S3. To render selected unchanged figures only, provide a comma-separated list:

```bash
VEGF_FIGURES=Figure2,Figure5 Rscript scripts/11_render_publication_figures.R
```

The renderers use included frozen TSV files and do not rerun differential expression or enrichment analysis.

## Rebuild expanded-validation summaries

The dataset-, program-, evidence-tier-, leave-one-dataset-out, and leave-one-context-out summaries can be regenerated from the bundled 72-row formal-validation result block:

```bash
python3 scripts/13_summarize_expanded_validation.py
```

A successful run ends with `EXPANDED_VALIDATION_SUMMARY=PASS`.

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

Continue in the order given in `run_manifest.tsv`, skipping entries marked `archived_provenance`. Candidate-dataset scripts use the external input layout documented in `DATA_ACCESS.md`. `VEGF_OUTPUT_ROOT` can redirect generated results without changing the archive.

## Reproducibility boundary

The figure-rendering and expanded-validation summary paths are directly runnable from this archive. Full expression-level reconstruction additionally requires downloading public GEO files and arranging them according to `DATA_ACCESS.md`. Dataset-specific ortholog handling, empirical-Bayes settings, and analysis seeds are documented in `metadata/dataset_methods_inventory.tsv` and the corresponding scripts. Package versions are recorded in `environment/software_versions.tsv`; the complete recorded R session is in `environment/R_sessionInfo.txt`. The expanded validation cohort was locked after discovery-program freezing and was not represented as prospectively prespecified.
