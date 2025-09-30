# Gridded copepod abundance
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub top language](https://img.shields.io/github/languages/top/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods) ![GitHub contributors](https://img.shields.io/github/contributors/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods) ![GitHub last commit](https://img.shields.io/github/last-commit/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods)     


## Introduction

Copepods are found everywhere in the ocean. They eat diatoms and phytoplankton, and constitute an essential link in the food web.
Gridded maps are useful for data visualisation, but also to assess the change over long time periods.

### The method

To create the gridded maps we will rely on [`DINCAE`](https://github.com/gher-uliege/DINCAE.jl) (Data-Interpolating Convolutional Auto-Encoder), a neural network designed to reconstruct missing data in satellite images, but that can be adapted to grid sparse, in situ measurements. 

## Directory structure

```
{{directory_name}}/
├── analysis
├── data/
│   ├── derived_data/
│   └── raw_data/
├── docs/
├── product/
└── scripts/
```

* **analysis** - Markdown or Jupyter notebooks
* **data** - Raw and derived data
* **docs** - Rendered reports
* **product** - Output product files
* **scripts** - Reusable code

## Data series

{{data_series}}

```
{{data_wfs_request}}
```

## Data product

{{data_product_description}}

## More information:

### References

### Code and methodology

{{link_code}}

### Citation and download link

The input data should be cited as:
> Pierre Helaouët (Marine Biological Association of the United Kingdom) (2025): CPR Data request - DINCAE - 29/09/2025. The Archive for Marine Species and Habitats Data (DASSH). (Dataset). https://doi.org/10.17031/68da4a97650f1

This product should be cited as:

{{product_citation}}

Available to download in:

{{link_download}}

### Authors

Charles Troupin (ULiège), Pierre Helaouët (MBA)      
Copyright (c) 2025 ULiège (GHER)
