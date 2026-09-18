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
│   │   └── candidate_validation/
│   │       ├── phase8c1/
│   │       ├── phase8c2/
│   │       └── platforms/
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
| Final external validation | GSE64052, GSE86525, GSE249415, GSE121153, GSE84048, GSE207976, GSE328515, and GSE78698 deposited files and required platform annotations; exact context definitions are in `config/external_validation_cohort.tsv` |
| Candidate-dataset sensitivity | GSE59476, GSE66346, GSE221557, GSE132568, and GSE80778 deposited files; final dispositions are in `config/candidate_dataset_disposition.tsv` |
| Supportive experimental evidence | GSE81465 IDAT files and GPL10558 manifest; optional GSE45161 Series Matrix and GPL9324 platform SOFT |
| Human longitudinal | GSE79671 deposited count matrix; `metadata/GSE79671_sample_manifest.tsv` |
| Human baseline association | GSE37138 Series Matrix; `metadata/GSE37138_sample_manifest.tsv` |
| GSE249415 entry point | Set `VEGF_GSE249415_INPUT` to the deposited `GSE249415_raw_HS.txt.gz` matrix |

Copy the included sample manifests into the matching metadata directories before expression-level reconstruction. The two evidence-role files in `config/` are the authoritative role definitions.

## Notes

- GSE76068 technical arrays were resolved to matched biological units before discovery analysis.
- GSE64472 human-array data were excluded after raw-data integrity review; the mouse-stromal data were retained.
- GSE81465 array-level directions are supportive-only because biological versus technical replicate identity could not be resolved.
- GSE45161 is qualitative directional support and does not enter the formal external-validation denominator.
- `VEGF_CANDIDATE_INPUT_ROOT` may override the default `data/raw/candidate_validation` location for scripts 11–12.
- `VEGF_R_LIBRARY` may point to a project-specific R library; it is optional.
- Public raw files may be large; GSE37138 is especially large and should be downloaded directly from GEO rather than redistributed.
