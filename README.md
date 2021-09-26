[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3688864.svg)](https://doi.org/10.5281/zenodo.3688864)

# trawllight
The trawllight R package contains functions to implement an algorithm that derives apparent optical properties from light measurements collected during NOAA bottom-trawl surveys. A description of the light data collection protocol, algorithm subroutines, and algorithm performance is provided in Rohan et al. (2020) and Rohan et al. (2021). Algorithm subroutines are contained in individual functions in the trawllight package.

The most recent version of the trawllight package was built in R 4.1.1.

# Installation

trawllight can be installed by starting R and running the following code. Installation requires the devtools package.

```
require(devtools)
install_github("sean-rohan-NOAA/trawllight")
```

# References
Rohan, S.K, Kotwicki, S., Britt, L.L., Laman, E.A., and Aydin, K. 2020. Deriving apparent optical properties from light measurements obtained using bottom-trawl-mounted archival tags. U.S. Dep. Commer., NOAA Tech. Memo. NMFS-AFSC-403, 91 p. [https://doi.org/10.25923/42yn-1q79](https://doi.org/10.25923/42yn-1q79)

Rohan, S.K., Kotwicki, S., Kearney, K.A., Schulien, J.A., Laman, E.A., Cokelet, E.D., Beauchamp, D.A., Britt, L.L., Aydin, K.Y., Zador, S.G., 2021. Using bottom trawls to monitor subsurface water clarity in marine ecosystems. Prog. Oceanogr. 194, 102554. [https://doi.org/10.1016/j.pocean.2021.102554](https://doi.org/10.1016/j.pocean.2021.102554)

### Legal disclaimer

This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.