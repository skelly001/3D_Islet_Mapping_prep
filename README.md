# 3D Islet Mapping — Spatial Proteomics

Pseudo-single-cell spatial proteomics of human pancreatic islet tissue. This repository contains the complete analysis pipeline and 2,215 protein abundance maps at cellular resolution.

**[Browse the interactive protein maps](https://skelly001.github.io/3D_Islet_Mapping_prep/)**

## Overview

This project maps protein expression in human pancreatic islets by combining:
- **CellPose 2.0** deep learning cell segmentation
- **FIJI/ImageJ** cell boundary refinement and pixel-to-cell mapping
- **Immunofluorescence** (INS, GCG, DAPI) for cell type classification (alpha, beta, acinar)
- **Nanopots mass spectrometry** for protein quantification
- **Azimuth** single-cell RNA-seq reference for cell-type deconvolution

## Pipeline

| Stage | Script | Description |
|-------|--------|-------------|
| 0 | `0a_1-3D_spatial_proteomics.R` | 3D spatial analysis with Giotto |
| 1 | `1b_1-ROI_mapping_and_cell_type_assignment.R` | Cell segmentation and type assignment |
| 2 | `2b_1-ROI_to_pixel.R` | ROI-to-pixel coordinate mapping |
| 3 | `3a_1-ROI_and_pixel_to_MS.R` | Mass spectrometry data integration |
| 4 | `4a_1-RNA-Seq.R` | RNA-seq deconvolution |
| 5 | `5a_1-cell_type_adjusted_protein_maps.R` | Final protein map generation |

## Output

2,215 protein abundance maps in `output/RD5-final_protein_maps/final_protein_maps/`, each showing cell-type-adjusted relative intensity (navy blue = low, yellow = high) across segmented cell polygons.

## Website

The companion website at [skelly001.github.io/3D_Islet_Mapping_prep](https://skelly001.github.io/3D_Islet_Mapping_prep/) features:
- Interactive slider comparisons between cell type maps and protein expression
- Searchable gallery of all 2,215 protein maps
- Pipeline visualization and methods summary
