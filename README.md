
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spheresmooth

This package offers a method of fitting a smooth path to a given set of
noisy spherical data observed at known time. It implements a piecewise
geodesic curve fitting method on the unit sphere based on a
velocity-based penalization scheme. The proposed approach is implemented
using the Riemannian block coordinate descent algorithm. To understand
the method and algorithm, you can refer to Bak, K. Y., Shin, J. K., &
Koo, J. Y. (2023)
[\<doi:10.1080/02664763.2022.2054962\>](https://www.tandfonline.com/doi/full/10.1080/02664763.2022.2054962)
for the case of order 1. Additionally, this package includes various
functions necessary for handling spherical data.

- Authorss:
  - Seyoung Lee, Sungshin Women’s University, <20210861@sungshin.ac.kr>
  - Kwan-Young Bak, professor at Sungshin Women’s University,
    <kybak@sungshin.ac.kr>,
    [ORCID:0000-0002-4541-160X](https://orcid.org/0000-0002-4541-160X%7D%7BORCID:0000-0002-4541-160X)

## Installation

You can install the development version of spheresmooth from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kybak90/spheresmooth")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(spheresmooth)
## basic example code
```
