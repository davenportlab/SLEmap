Git repo to accompany SLEmap manuscript

Authors
- Haerin Jang
- Catherine Sutherland
- Wanseon Lee

Corresponding Author
- Emma Davenport

Structure
- 1_YASCP: single cell data demultiplexing, Azimuth cell type annotation, doublet detection
- 2_singlecell_QC_annotation: single cell data QC, clustering, manual annotation
- 3_variant_QC: filter genotype file for eQTL mapping, ancestry assignment
- 4_TCRBCR: TCR BCR data processing
- 5_eQTL mapping: identify optimal number of expression PCs, map eQTLs by cell type, flanders pipeline
- 6_mashr: identify proportion of eQTLs shared between cell types
- 7_colocalization: colocalization analysis between eQTL and GWAS
- 8_ancestry: ancestry interaction eQTL mapping and eQTL mapping by ancestry
- 9_Onek1k: reannotate onek1k cell type, map eQTLs, colocalize with GWAS, and compare with SLEmap