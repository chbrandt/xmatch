Maximum Likelihood Estimator

The MLE method, learnt from [LaMassa et al. 2016], takes into account
the distance between the primary source and the ancillary objects,
their astrometic errors and magnitude distribution. Ancillary objects
will be the sources found within the search radius in the second catalog.

In MLE it is defined the Reliability value, `R`, which is a way to
distinguish between true counterparts and chance detections. The greater
the `R`, greater are the chances to be a true counterpart. `R` is computed
as follows:
$$
    R = \frac{LR}{\Sum_i LR_i + (1 - Q)}
$$

`Q` is the ratio of primary sources that present ancillary objects within
the search radius [DOUBT-1] divided by the total number of (primary) sources.

`LR` is the likelihood ratio, it is assigned to each ancillary object
within the search radius and it is the probability that one of those
objects is the true counterpart divided by the probability that a background
source is there by chance.
$$
    LR = \frac{q(m) f(r)}{n(m)}
$$

$q(m)$ is the expected normalized magnitude distribution of counterparts
within the radius search. $q(m)$ is estimated by subtracting the histogram
of the sources inside the the radius from the histogram of background
sources; the histograms are, each, normalized by the respective search
areas. Background objects are taken from an annulus around the primary
source, where the internal radius should be larger than the search radius.

$f(r)$ is the probability distribution of the astrometric errors; it is
modeled as a two-dimensional Gaussian distribution where the primary
ancillary positional errors are added in quadrature.

$n(m)$ is the normalized magnitude distribution of background sources.


LaMassa et al. 2016: http://doi.org/10.3847/0004-637X/817/2/172
