# CRAN comments

This is a new version of an existing package. It updates the package from version 0.1.0 to 0.1.1.

## Changes in version 0.1.1

None of the changes in version 0.1.1 alter the functionality of the package. They primarily address deprecated behavior and improve metadata. Specifically:

* To prepare for the future deprecation of the `fortify` function, we have updated the world map visualizations by migrating to the `sf` package.
* All instances where `fortify` was used have been replaced with `sf` transformations, which also improve the quality of the visualizations.
* A package website and GitHub repository have been added to the `URL` field in the package description.
* The GitHub repository has also been added to the `BugReports` field in the package description.

## R CMD check results

0 errors | 0 warnings | 0 notes
