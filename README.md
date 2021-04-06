# UQ-Project2

## How to use this code

Set `/src` as the working directory within R studio. Then run each line of `GPExample.R` one at a time. (Or all at once if you like!)

## Required Packages

### RSpectra

This is used for getting the eigenvalues and eigenvectors for the Karhunen-Loeve expansion.

```R
install.packages("RSpectra")
```

## Recommended Packages

### fields

This package is not required but has some nice functions for plotting functions for 2D GP samples.

```R
install.packages("fields")
```