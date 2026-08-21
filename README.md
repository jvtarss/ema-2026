# Eucalyptus MicroRNA Archive (EMA)

A curated database of microRNAs in *Eucalyptus grandis*, integrating three independent public small RNA sequencing studies (vegetative tissue, somatic embryogenesis, tension wood formation) into a single, locus-resolved catalog with study-level evidence tracking.

**Current release:** v1.0 , with 99 curated miRNAs (31 known, 68 novel), 34 family-level groupings, 1,773 predicted target interactions.

## Citation

If you use EMA in your research, please cite:

> [FULL CITATION PENDING — add once published]

## Repository structure

```
ema-dashboard/
├── data/
│   ├── catalog.csv                          curated miRNA catalog (99 entries: accession, sequence, family, confidence tier, coordinates)
│   ├── evidence_table.csv                   per-study discovery evidence (canonical, same-mature-sequence, same-identifier, same-locus-overlap)
│   ├── expression_matrix.csv                raw counts + CPM, all 20 samples
│   ├── differential_expression_results.csv  DESeq2 results, QIN-2021 and TOLENTINO-2022 contrasts
│   ├── target_predictions.csv               1,773 psRNATarget miRNA-target interactions
│   ├── network_nodes.csv                    node list for the miRNA-target-PPI network
│   ├── network_edges_mirna_target.csv       miRNA-target edges (cleavage/translation)
│   └── appendix_a_sample_metadata.csv       full per-sample metadata
├── scripts/
│   ├── export_curated_files.sh              regenerates all files in data/ from the SQLite database
│   ├── 01_process_string_ppi.py             harmonizes and filters STRING PPI export
│   ├── 02_build_network.py                  merges miRNA-target and PPI edges, classifies node status
│   ├── 05_generate_node_labels.py           generates short display labels and common names for protein nodes
│   ├── 06_functional_categories.py          assigns broad functional categories from STRING GO enrichment
│   └── 07_group_mirnas_by_target_go.py      groups miRNAs by the dominant GO term of their target set
├── db/
│   └── ema_v1.db                            SQLite database (source of truth for data/)
└── webapp/                                  dashboard frontend/backend source
```

## Data provenance

Source studies (see manuscript Methods section 2.1 and Table 1 for full sample-level detail):

| Study | Focus | SRA BioProject / accessions |
|---|---|---|
| LIN-2018 | Vegetative tissue (leaf, stem) | SRR5433935, SRR5433936 |
| QIN-2021 | Somatic embryogenesis (stem vs. callus, two genotypes) | SRR11749594–SRR11749605 |
| TOLENTINO-2022 | Mechanically induced tension wood | SRR14138141–SRR14138146 |

## Reproducing the curated data files

All files in `data/` are generated directly from `db/ema_v1.db`:

```bash
bash scripts/export_curated_files.sh
```
## License

[PENDING — to be finalized once the target journal/venue is confirmed]

## Contact

For questions, bug reports, or to propose additional datasets for a future EMA release, contact the corresponding author: [EMAIL PENDING]
