# Great circle cross-matching algorithm

The fundamental characteristic of a great-circle (gc) cross-matching is the
restriction on the region to search for counterparts, and that limiting region
is defined by a *radius* parameter around each target source.

Whereas a *nn* algorithm will *always* find a counterpart for each and every
target, the matching pair may be separated by a unreasonable distance.
The *gc* algorithm on the other hand, by restricting the search region, may not
find a counterpart for each target; and depending on the size of the search
region, it may find multiple counterparts, to decide which one to take is a
second step.

Typically, in the simplest scenario, the *nearest neighbour* is taken among the
counterpart candidates found by the *great-circle* algorithm.

The great circle approach provides a better performant algorithm (O(N*logN))
and provide a natural mechanism to sub-sample the datasets.
Such sampling mechanism may then be used to dynamically define the best match.
