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
|   ├── aux_data/
│   ├── derived_data/
│   └── raw_data/
├── docs/
├── product/
└── scripts/
```

* **analysis** - Jupyter notebooks for the preparation of the input files and for the interpolation 
* **data** - Raw data (CPR observations), auxiliary data (bathymetry, coastline, etc) and derived data (data files to be ingested by `DINCAE`)
* **docs** - Rendered reports
* **product** - Output product files and figures
* **scripts** - Reusable code (main module containing functions used in the **analysis**)

## Data series

Raw data have to be downloaded from https://doi.mba.ac.uk/data/3548/

and then extracted in the corresponding directory (`./data/raw_data/`)
```bash
unzip CPR_DataRequest_DINCAE_29Sep25.zip
```

Three files are obtained:
1. CPR_Data_DINCAE_290925.docx: the documentation,
2. CPR_DINCAE_ControlMap_29092025.png: a map representing the selected samples,
3. CPR_DINCAE_Data_290925.csv: the selected samples.

## Data product

The data product consists of a netCDF file storing gridded fields for two variables: the _small copepods_ and the _large copepods_
.
## More information:

### References

### Code and methodology

`DINCAE` software is available at https://github.com/gher-uliege/DINCAE.jl

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
