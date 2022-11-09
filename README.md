[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3688864.svg)](https://doi.org/10.5281/zenodo.3688864)

# trawllight
The trawllight R package contains functions to implement an algorithm that derives apparent optical properties from light (downwelling irradiance) measurements collected during NOAA Alaska Fisheries Science Center bottom-trawl surveys using trawl-mounted archival tags. The light data collection protocol, algorithm subroutines, and algorithm performance are described in Rohan et al. (2020) and Rohan et al. (2021).

<br><br>
![](./assets/transect_S_2011_2017_wide.png)

<i>Transect along 60°N on the eastern Bering Sea shelf in 2017 showing the downwelling diffuse attenuation coefficient <i>K<sub>d</sub></i>(<i>z, tag</i>), mixed layer depth (MLD), bottom layer depth (BLD), and depth where downwelling irradiance is equal to 10% of irradiance and just beneath the sea surface (Z<sub>10%</sub>, i.e., optical depth = 2.3).</i>


# Installation

The trawllight R package can be installed from GitHub as follows:

```
require(remotes)
remotes::install_github("afsc-gap-products/trawllight")
```

# Documentation

- [Decode and synchronize data](1_process_mk9.Rmd)
- [Run trawllight: QA/QC and derive apparent optical properties](2_run_trawllight.Rmd)
- [Generate data prodcts](3_make_data_product.Rmd)

# References

## Algorithm and validation
Rohan, S.K, Kotwicki, S., Britt, L.L., Laman, E.A., and Aydin, K. 2020. Deriving apparent optical properties from light measurements obtained using bottom-trawl-mounted archival tags. U.S. Dep. Commer., NOAA Tech. Memo. NMFS-AFSC-403, 91 p. [https://doi.org/10.25923/42yn-1q79](https://doi.org/10.25923/42yn-1q79)

Rohan, S.K., Kotwicki, S., Kearney, K.A., Schulien, J.A., Laman, E.A., Cokelet, E.D., Beauchamp, D.A., Britt, L.L., Aydin, K.Y., Zador, S.G., 2021. Using bottom trawls to monitor subsurface water clarity in marine ecosystems. Prog. Oceanogr. 194, 102554. [https://doi.org/10.1016/j.pocean.2021.102554](https://doi.org/10.1016/j.pocean.2021.102554)

## Applications
Spies, I., Tarpey, C., Kristiansen, T., Fisher, M., Rohan, S., and Hauser, L. 2022. Mechanisms for genomic differentiation in Pacific cod using Pool-Seq. Evolutionary Applications. [https://doi.org/10.1111/eva.13488](https://doi.org/10.1111/eva.13488)


# Legal disclaimer

This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
