# Gridded copepod abundance
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub top language](https://img.shields.io/github/languages/top/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods) ![GitHub contributors](https://img.shields.io/github/contributors/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods) ![GitHub last commit](https://img.shields.io/github/last-commit/gher-uliege/EMODnet-Biology-Interpolation-NN-Copepods)     


## Introduction

Copepods are found everywhere in the ocean. They eat diatoms and phytoplankton, and constitute an essential link in the food web.
Gridded maps are useful for data visualisation, but also to assess the change over long time periods.

Two taxonomic groups are considered in this work:
- the _small copepods_: they correspond to the mean values of the 46 small copepods taxa (<= 2mm) found in the selected area;
- the _large copepods_: they correspond to the mean values of the 110 large copepods taxa (> 2mm) found in the selected area.

```geojson
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "coordinates": [
          [
            [
              -30.0,
              67.0
            ],
            [
              -30.0,
              42.0
            ],
            [
              9.0,
             42.0
            ],
            [
              9.0,
              67.0
            ],
            [
              -30.0,
              67.0
            ]
          ]
        ],
        "type": "Polygon"
      }
    }
  ]
}
```

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
1. `CPR_Data_DINCAE_290925.docx`: the documentation,
2. `CPR_DINCAE_ControlMap_29092025.png`: a map representing the selected samples,
3. `CPR_DINCAE_Data_290925.csv`: the selected samples.

## Data product

The data product consists of a netCDF file storing gridded fields for two variables: the _small copepods_ and the _large copepods_.

## More information:

### References

### Code and methodology

`DINCAE` software is available at https://github.com/gher-uliege/DINCAE.jl.     
The code is described in the following papers:

> Barth, A., Alvera-Azcárate, A., Troupin, C., & Beckers, J.-M. (2022). DINCAE 2.0: multivariate convolutional neural network with error estimates to reconstruct sea surface temperature satellite and altimetry observations. _Geoscientific Model Development_, __15(5)__, 2183–2196. https://doi.org/10.5194/gmd-15-2183-2022
<details>

<summary>BibTeX entry</summary>

```bash
@Article{gmd-15-2183-2022,
AUTHOR = {Barth, A. and Alvera-Azc\'arate, A. and Troupin, C. and Beckers, J.-M.},
TITLE = {DINCAE 2.0: multivariate convolutional neural network with error estimates to reconstruct sea surface temperature satellite and altimetry observations},
JOURNAL = {Geoscientific Model Development},
VOLUME = {15},
YEAR = {2022},
NUMBER = {5},
PAGES = {2183--2196},
URL = {https://gmd.copernicus.org/articles/15/2183/2022/},
DOI = {10.5194/gmd-15-2183-2022}
}
``` 

</details>

<br>

> Barth, A., Alvera-Azcárate, A., Licer, M., & Beckers, J.-M. (2020). DINCAE 1.0: a convolutional neural network with error estimates to reconstruct sea surface temperature satellite observations. _Geoscientific Model Development_, __13(3)__, 1609–1622. https://doi.org/10.5194/gmd-13-1609-2020

<details>

<summary>BibTeX entry</summary>

```bash
@Article{gmd-13-1609-2020,
AUTHOR = {Barth, A. and Alvera-Azc\'arate, A. and Licer, M. and Beckers, J.-M.},
TITLE = {DINCAE 1.0: a convolutional neural network with error estimates to reconstruct sea surface temperature satellite observations},
JOURNAL = {Geoscientific Model Development},
VOLUME = {13},
YEAR = {2020},
NUMBER = {3},
PAGES = {1609--1622},
URL = {https://gmd.copernicus.org/articles/13/1609/2020/},
DOI = {10.5194/gmd-13-1609-2020}
}
``` 

</details>

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

## Procedure

1. Download the data file from https://doi.org/10.17031/68da4a97650f1.
   <details>   

   The zip file has to be extracted in the corresponding directory (`./data/raw_data/`). Three files are obtained: 
   1. `CPR_Data_DINCAE_290925.docx`: the documentation,
   2. `CPR_DINCAE_ControlMap_29092025.png`: a map representing the selected samples,
   3. `CPR_DINCAE_Data_290925.csv`: the selected samples. 
   </details>

2. Edit and run `read_copepods_data.ipynb` to convert the CSV to netCDF files.
   <details>

   This will produce one netCDF file per group: `Small_copepods_DINCAE_$(regionname).nc` and `Large_copepods_DINCAE_$(regionname).nc` (where `regionname` is defined in `param.jl`).

   </details>
  
3. Edit and run `prepare_validation_data.ipynb` to prepare the data files that will be used for the training and for the validation.
   
4. Run `prepare_environment_variables.ipynb` to prepare the environmental data.
   
   <details>
   
   Presently we work with the following variables:
   - bathymetry, 
   - distance to nearest coastline and 
   - sea surface temperature.
  Other variables can be considered (depending on the application).
    </details>
   
5. Run the script `run_DINCAE_random_copepods.jl` to generate analysis with a random set of _hyperparameters_. 
    <details>

    > [!NOTE]
    > Ideally, the code should run with GPU, otherwise it will be very slow. You will get the following message:

    ```julia
    ┌ Warning: No supported GPU found. We will use the CPU which is very slow. Please check https://developer.nvidia.com/cuda-gpus
    ```
    </details>

6. Run the script `find_best_param.jl` to find the combination of hyperparameters and environment variables that yield the best results.

7. Run the script `run_best_params.jl` to generate the _best_ analysis using the full dataset.

## Acknowledgements

Computational resources have been provided by the Consortium des Équipements de Calcul Intensif ([CÉCI](https://www.ceci-hpc.be/)), funded by the Fonds de la Recherche Scientifique de Belgique ([F.R.S.-FNRS](https://frs-fnrs.be)) under Grant No. 2.5020.11 and by the Walloon Region.