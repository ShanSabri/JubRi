# JubRi
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) 
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/ShanSabri/deconR/graphs/commit-activity) 
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.1-6666ff.svg)](https://cran.r-project.org/) 
[![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/master) 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-ff69b4.svg)](https://www.gnu.org/licenses/gpl-3.0)

JubRi is a R package that performs literature-based keyword ontology. JubRi utilizes a curated keyword database of roughly 18,000 genes queried using the [`PubMedScrapeR`](https://github.com/ShanSabri/PubMedScrapeR) R package. 

## Installation

Use [devtools](https://github.com/r-lib/devtools) to install JubRi directly from GitHub.

```R
# if(!require(devtools)) {install.packages(devtools)}
devtools::install_github("ShanSabri/JubRi")

```

## Usage

See [`example/example.R`](https://github.com/ShanSabri/JubRi/blob/master/example/example.R) on how to generate a metaplot of term enrichments given a query set, shown below: 
<img src="man/figures/metaplot.png" align="center" alt="" width=800/>

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[ GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)

