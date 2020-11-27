[![DOI](https://zenodo.org/badge/DOI/?.svg)](https://zenodo.org/record/?)
[![Python 3.7](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-369/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# Diverse dynamical response to orographic gravity wave drag hotspots - a zonal mean perspective
**P. Sacha, A. Kuchar, R. Eichinger, P. Pisoft, Ch. Jacobi and H. Rieder**

Code used to process and visualise the model and other data outputs in order to reproduce figures in the manuscript.
Model data are available [here](http://climate-modelling.canada.ca/climatemodeldata/cmam/output/CMAM/CMAM30-SD/index.shtml). All datasets already preprocessed can be found [here](https://data.mendeley.com/datasets/j3hj7f9t67/3).

Notebooks for each individual figure as well as for two data tables are in the [`code/` directory](code), while the figures themselves are in the [`plots/` directory](plots).

### Figures
|  #  | Figure                                                                                                                                                                                                    | Notebook                                                                              | Dependencies                                                                                                                                                             |
|:---:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|  1 | [EPFD composite average at lag=0 for the HI, EA and WA hotspot, respectively](plots/EPFD+EPfluxes_anomalies_all_20days_zm_wEPFDsignificancetropopause_DJFonly+alllayers.pdf)                                                                              | [GRL_reproduce_Fig1.ipynb](code/GRL_reproduce_Fig1.ipynb)                       |                                                                                                                                       |
|  2 | [Decomposition into leading zonal wavenumbers of Eliassen-Palm flux composite average at lag=0 for the HI, EA and WA hotspot, respectively](plots/EPFD+EPfluxes+wavenumbers-123_anomalies_all_20days_zm_wEPFDsignificancetropopause_DJFonly_pvalue0.05.pdf)                                                      | [GRL_reproduce_Fig1.ipynb](code/GRL_reproduce_Fig2.ipynb)                 |                                                                                                                           |
|  3 | [Zonal wind composite at lag=0 representing the HI, EA and WA hotspot, respectively](plots/ua_anomalies_all_20days_zm_wsignificance_DJFonly_Berlin_cm.pdf)                | [GRL_reproduce_Fig3.ipynb](code/GRL_reproduce_Fig3.ipynb)                 | [ssw_composite_cmam_optimized2-wss_hotspots.py](code/ssw_composite_cmam_optimized2-wss_hotspots.py), [montecarlo_composites_script.py](code/montecarlo_composites_script.py)                                                                                                                            |
|  4 | [Zonal mean FAWA composite averages at lag=0 representing the HI, EA and WA hotspot, respectively](code/lwatend_anomalies_all_20days_zm_wsignificance_DJFonly.pdf) | [GRL_reproduce_Fig4.ipynb](code/GRL_reproduce_Fig4.ipynb)                     |                                                                                                                            |

#### Supplementary figures
|  #  | Figure                                                                                                                                                                                                    | Notebook                                                                              | Dependencies                                                                                                                                                             |
|:---:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|  S1 | [Residual anomalies during boreal winter representing the HI, EA and WA hotspot, respectively](plots//TEM-res3-new_anomalies_all_20days_zm_wosignificance_DJFonly.pdf)                                               | [GRL_reproduce_FigS1.ipynb](code/GRL_reproduce_FigS1.ipynb)                     | |
| S2 | [{Surface downward eastward wind stress anomalies during boreal winter representing the HI, EA and WA hotspot, respectively](plots/tauu_anomalies_allwclim_20days_wsignificancefrom10000_PlateCarree_DJFonly.pdf)                                                                              | [GRL_reproduce_FigS2.ipynb](code/GRL_reproduce_FigS2.ipynb)                       |                                                                                                                                    |
|  S3 | [Zonal mean of zonal wind anomalies (shading; units: m/s) at lag=0 representing the Himalaya hotspot composites differentiated according to QBO phases](plots/ua_anomalies_all_20days_zm_wosignificance_DJFonly_QBO_Himalyasonly.pdf)                                                                              | [GRL_reproduce_FigS3.ipynb](code/GRL_reproduce_FigS3.ipynb)                       |             [ssw_composite_cmam_optimized2-wss_hotspots_QBO.py](code/ssw_composite_cmam_optimized2-wss_hotspots_QBO.py),           [
GRL_QBO_timeseries4composites_CMAM.ipynb ](code/
GRL_QBO_timeseries4composites_CMAM.ipynb )                                                                                                            |

### Required package installation
`pip install -r requirements.txt`

### References

Kuchar, A., Sacha, P., Eichinger, R., Jacobi, C., Pisoft, P., and Rieder, H. E.: On the intermittency of orographic gravity wave hotspots and its importance for middle atmosphere dynamics, Weather Clim. Dynam., 1, 481-495, [https://doi.org/10.5194/wcd-2020-21](https://doi.org/10.5194/wcd-1-481-2020), 2020.
