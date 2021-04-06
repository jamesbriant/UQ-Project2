# Karhunen-Loeve implementation

# FIRST SET THIS DIRECTORY AS THE WORKING DIRECTORY.

source("include/GPFunctions2D.R")


# Try I=J=35 first, this isn't ideal but it looks okay.
# Note a I=J=100 created a matrix that is 10000x10000, this takes up 800MB of RAM...
#   I recommend only doing this if you have 16GB RAM installed.
I <- 35 # number of points in the x-direction
J <- I  # number of points in the y-direction
N <- 4  # number of samples

# Build the grid
x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)

# Get the Karhunen-Loeve decomposition.
# THIS IS *VERY* SLOW, but only needs running once.
KLDecomp <- GetKLDecomposition(I=I, J=J)

# This generates new samples from the KL decomposition.
# This function can be run as many times as you like, it's very quick. 
# Increase N if you like! eg. N=200 is still very fast
samples <- GenerateSamples(N, KLDecomp)

# the package `fields` is required for this.
par(mfrow=c(2,2))
fields::image.plot(x, y, samples[[1]])
fields::image.plot(x, y, samples[[2]])
fields::image.plot(x, y, samples[[3]])
fields::image.plot(x, y, samples[[4]])
par(mfrow=c(1,1))

# Alternatively, if `fields` is not installed.
persp(x, y, samples[[1]])
