# Public data access and expected layout

All expression datasets used by the code are publicly available through NCBI Gene Expression Omnibus. Accession pages and analysis roles are listed in `metadata/public_data_sources.tsv`. Raw and deposited processed expression files are not redistributed in this archive.

## Analysis root

Set `VEGF_ANALYSIS_ROOT` to a directory containing the following structure:

```text
analysis_directory/
├── data/
│   ├── metadata/sample_manifest.tsv
│   ├── processed/
│   ├── raw/
│   └── raw_round1_signoff/GSE180687/extracted/
├── results/
├── reports/
├── round3a_validation/
│   └── data/raw/
└── round3b_human_translation/
    └── data/raw/
```

The included scripts expect these principal inputs:

| Analysis | Required public inputs |
|---|---|
| Discovery | GEO Series Matrix or deposited normalized matrices for GSE76068, GSE73571, GSE64472, and GSE26644; GSE180687 CEL files; applicable GEO platform annotations; `metadata/discovery_sample_manifest.tsv` |
| Experimental validation | GSE81465 IDAT files and GPL10558 manifest; GSE64052 Series Matrix and GPL570 annotation; GSE86525 Series Matrix and GPL16699 platform SOFT; optional qualitative GSE45161 Series Matrix and GPL9324 platform SOFT |
| Human longitudinal | GSE79671 deposited count matrix; `metadata/GSE79671_sample_manifest.tsv` |
| Human baseline association | GSE37138 Series Matrix; `metadata/GSE37138_sample_manifest.tsv` |
| Additional sensitivity | GSE249415 deposited `GSE249415_raw_HS.txt.gz` matrix |

Copy the included sample manifests into the matching metadata directories before expression-level reconstruction. The two evidence-role files in `config/` are the authoritative role definitions.

## Notes

- GSE76068 technical arrays were resolved to matched biological units before discovery analysis.
- GSE64472 human-array data were excluded after raw-data integrity review; the mouse-stromal data were retained.
- GSE81465 array-level directions are supportive-only because biological versus technical replicate identity could not be resolved.
- GSE45161 is qualitative directional support and does not enter the primary-validation denominator.
- Public raw files may be large; GSE37138 is especially large and should be downloaded directly from GEO rather than redistributed.

