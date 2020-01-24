# trawllight
The trawllight R package contains functions to implement an algorithm that derives apparent optical properties from light measurements collected during NOAA bottom-trawl surveys. A description of the light data collection protocol, algorithm subroutines, and algorithm performance is presented in Rohan et al. (in prep). Algorithm subroutines are contained in individual functions in the trawllight package.

The trawllight pacakge was built using R 3.6.1.


# Installation

trawllight can be installed by starting R and running the following code. Installation requires the install_github() function from the devtools package.

```
require(devtools)
install_github("sean-rohan/trawllight")
```

# Vignette

A vignette demonstrating how to use the trawllight package is in development (SR Jan 24, 2020).

# References
Rohan, S., Kotwicki, S., Laman, E., Britt, L., Aydin, K. In review. Deriving apparent optical properties from light measurements obtained using bottom-trawl-mounted archival tags.