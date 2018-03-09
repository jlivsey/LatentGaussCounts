## LatentGaussCounts

Modeling Integer-valued Time Series via Latent Gaussian Process

### Installation

As the package is not yet on CRAN, it needs to be installed from GitHub:

    library(devtools)  # if you don't have devtools: install.packages("devtools")
    install_github("jlivsey/LatentGaussCounts")

### Utility of LatentGaussCounts

- Fit count time series model for most classic discrete distributions
- Currently we support the following discrete dstributions
	- Poisson
	- mixture of Poissons
- Allow autocorrelation to depend on an underlying Gaussian process.
- Currently we support the following underlying Gaussian processes
	- AR(1)
	- FARIMA(0,d,0)
- Estimate model parameters with Gaussian likelihood of a particle filtering MLE.


### Authors

[James Livsey](http://www.census.gov/research/researchers/profile.php?cv_profile=3922&cv_submenu=title) (United States Census Bureau) and
[Vladas Pipiras](http://pipiras.web.unc.edu/)(University of North Carolina, Chapel Hill) 


### License

TBD