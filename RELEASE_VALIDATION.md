# Release validation

- Archive structure and checksum coverage: PASS
- R script syntax: PASS
- Python script syntax: PASS
- Personal absolute paths: none detected
- Assistant or automated-author attribution: none detected
- Credential-like strings: none detected
- Temporary/cache files: none retained
- Unresolved submission placeholders: none retained
- Expanded formal-validation block: 8 datasets, 12 contexts, and 72 dataset-context-program rows
- Expanded-validation summary regeneration: exact for all 5 bundled summary tables
- Figure source-data coverage: complete for all main and supplementary figures in this archive
- PNG reproduction: exact for the 8 figures carried forward without layout adjustment; Figure 1, Figure S1, and Figure S4 were successfully rendered and checked for content consistency, but the final submission exports include layout-only crop/margin adjustments and are therefore not asserted to be pixel-identical to direct renderer output
- PDF and SVG rendering: completed without errors; the same layout-only export qualification applies to Figure 1, Figure S1, and Figure S4
- Evidence-role configuration: included
- Evidence-tier taxonomy: two comparator tiers only (Tier A / Tier B); dataset-specific qualifications retained separately; supportive-only evidence remains outside the formal denominator
- Authoritative cross-stage sample/QC manifest: included (`metadata/all_sample_inclusion_QC.tsv`)
- Formal external-validation sample/QC subset: 93 records across 8 datasets (`metadata/external_validation_sample_manifest.tsv`)
- Superseded Phase 9A validation manifest: explicitly retained as archival provenance
- Dataset-specific ortholog handling, empirical-Bayes settings, Hallmark scope, and seeds: aligned with the executed scripts and methods inventory
- GSE45161 supportive inventory: 4 qualitative contrasts; excluded from the formal validation denominator
- Robustness table: publication-facing four-column table aligned with the workbook; release-record availability documented separately
- External dataset accession links: included
- Phase 10C manuscript package modified during Phase 10C.1 preparation: no
- Expression-level analyses or enrichment rerun during Phase 10C.1 preparation: no
- Publication figures modified during Phase 10C.1 preparation: layout-only corrections to Figure 1, Figure S1, and Figure S4; no numerical or scientific content changed

The publication-figure renderers and expanded-validation summary builder are directly runnable from the archive. Direct rendering reproduces the underlying figure content; the final submission versions of Figure 1, Figure S1, and Figure S4 additionally contain layout-only crop/margin adjustments. Expression-level analyses require the public GEO inputs described in `DATA_ACCESS.md`. The archived Phase 9A scripts are retained for provenance and are not part of the current Phase 10C.1 execution path.
