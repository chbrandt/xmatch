
[kdtree]: https://en.wikipedia.org/wiki/K-d_tree

# Great circle algorithm

The basic characteristic of a *great-circle* (hereafter `GC`) algorithm is the
restriction on the region to search for counterparts, such region is limited
by a *radius* parameter around each object in the target catalog.

Besides adding a fundamental premise to the process -- i.e, it is worthless to
search for non nearby objects --, `GC` performs computationaly better as the
search space around each target is considerably small.
From the computational point-of-view, such restriction allows us to make use of
optimal data-structures to correlate our catalogs.

The great-circle method may find multiple counterpart candidates inside the
search radius around the target(s), and so a flitering process has to eventually
be applied to decide which one of the candidates is to be considered the right match.
In the simplest scenario, the *nearest neighbour* will be considered the right match.
This approach is typically used when matching catalogs of same wavebandas,
similar positional errors or sensitivity.

An extension to the nearest-neighbour approach to select the best candidate when
great-circle matches multiple candidates will be discussed in the next section,
`Maximum Likelihood Estimator`, where other objects features (e.g, luminosity)
are considered on estimating the correct counterpart.

Here, below, I take the chance to review the (computational) data-structure working
under a search-by-position algorithm.


## Space-partitioning trees
