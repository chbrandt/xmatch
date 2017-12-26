
# Nearest Neighbor algorithm

Nearest-neighbour is the fundamental idea of cross-matching.
Effectively, an object in the sky seen by two different instruments, can only be said to be
be the same if at least their (registered) position in each catalog is spatialy close enough.

The first approach to the solution is to consider the least distant pair of objects in
between both catalogs to be the correct match -- i.e, the same object.
Although we may foresee flaws in the approach -- for instance, the nearest neighbour may
actually be too far to be the same object --, it is a reasonable idea to consider.

Here the nearest-neighbour algorithm (hereafter `NN`) is furhter discussed to highlight
the subtleties of the implementation.
For instance, the `NN` implemented performs a bi-directional matching so that the
resulting matches hold a one-to-one relation, which is the only possible physically
reasonable answer.

In what follows we will consider two artifical list of objects positions,
the tables **A**, **B**.

`Figure 1`_ illustrates the spatial distribution of objects from tables *A*,*B*.
The size of the circles express the *positional error radius*, ``20 pixels``.

![Figure 1: Mock spatial distribution of objects.](images/xmatch/TOY_artificial_obj_distro_2_explain_matchAB_labeled.png)


## The search for neighbors

The algorithm does a bi-directional search, first taking catalog *A* as the *target catalog*
and searching for the nearest neighbor(NN) in *B* (*B->A*), then for the NN of *B* in *A*
(*A->B*).
Because there can be multiple matches with the same NN (different sources in *A* can have the
same object in *B* as NN, when *B->A* for example) a cleaning step is necessary to remove this
multiplicity.
Such cleaning is done by keeping only the pair of objects (*i.e.*, match) that exist on both
matching lists.

Searching for the nearest neighbor to the objects in catalog *A* we find the objects
in *B* as listed on the  table `Matching (NN) B->A`_ below.
On the other hand, looking for the NN in the other direction, *B* as reference and
looking for neighbors in *A*, we have the matches shown in `Matching (NN) A->B`_.
As we might already been expecting, multiple objects were found from the second
("neighborhood") catalog to match.

Next step is to clean those multiplicities. Going through the *multiple matches*
given from both tables we should agree that the real *nearest neighbors* are:

* **B->A**: [A:3,B:1], [A:7,B:4], [A:4,B:2], [A:12,B:5]
* **A->B**: [B:2,A:4]

That is to say that all other multiple matches should be taken out.
Doing such, we end up with `Table 3`_, where separation (in pixels) are shown,
where indeed only those 4 matchings-pairs are kept for the final result.

Notice, though, that we are not getting into the question whether the matches
are (physically) right or wrong.
The current discussion regards only the a search for a unique set of *nearest neighbors*
pairs as a basic idea to establish before going further on qualifying those matches.

Table 2.1: Matching (NN) B->A

A_ID | x | y | NN_in_B
--- | --- | --- | ---
1 | 301 | 34 | 1
2 | 51 | 74 | 1
3 | 230 | 145 | 1
4 | 404 | 232 | 2
5 | 229 | 265 | 4
6 | 52 | 286 | 4
7 | 83 | 288 | 4
8 | 346 | 289 | 2
9 | 317 | 339 | 2
10 | 214 | 376 | 4
11 | 465 | 388 | 5
12 | 401 | 455 | 5


Table 2.2: Matching (NN) A->B

B_ID | x | y | NN_in_A
--- | --- | --- | ---
1 | 299 | 126 | 3
2 | 391 | 246 | 4
3 | 445 | 249 | 4
4 | 120 | 359 | 7
5 | 428 | 441 | 12
6 | 125 | 446 | 10


Table 3: NN matching-table

A_ID | x | y | B_ID | x | y | dist
--- | --- | --- | --- | --- | --- | ---
1 | 301 | 34 | nan | nan | nan | nan
2 | 51 | 74 | nan | nan | nan | nan
3 | 230 | 145 | 1 | 299 | 126 | 71
4 | 404 | 232 | 2 | 391 | 246 | 19
5 | 229 | 265 | nan | nan | nan | nan
6 | 52 | 286 | nan | nan | nan | nan
7 | 83 | 288 | 4 | 120 | 359 | 80
8 | 346 | 289 | nan | nan | nan | nan
9 | 317 | 339 | nan | nan | nan | nan
10 | 214 | 376 | nan | nan | nan | nan
11 | 465 | 388 | nan | nan | nan | nan
12 | 401 | 455 | 5 | 428 | 441 | 30



### Pseudocode

The (pseudocode) inputs are ``catalog-A`` and ``catalog-B``, tables where at least
positional columns ``RA`` (Right Ascension), ``DEC`` (Declination) and an identifier
column ``ID`` (Object ID in each catalog) are given.

```
SkyCoord := ~astropy.coordinates.SkyCoord

A <- read columns 'RA','DEC','OBJID' from catalog-A
B <- read columns 'RA','DEC','OBJID' from catalog-B

coords_A <- build SkyCoord objects from 'RA','DEC'
coords_B <- build SkyCoord objects from 'RA','DEC'

match_A <- find nearest-neighbor from A in B
match_B <- find nearest-neighbor from B in A

match_AB <- filter intersection matches between match_A & match_B

matched_catalog <- join A,B following match_AB
```

The output, ``matched_catalog``, provides all columns in ``A`` and ``B``
-- for instance, ``RA``, ``DEC``, ``OBJID`` -- plus the angular separation
-- ``distance_AB`` -- between the matches.
The number of rows is the same as the given target catalog (i.e, *catalog A*).
Using the *relational database* jargon, the ``matched_catalog`` is the result of
a *left join* between tables ``A``, ``B``.

|     | A   |     |     | B   |     |     |  AB
| --- | --- | --- | --- | --- | --- | --- | ---
| index | OBJID | RA | DEC | OBJID | RA | DEC | dist
| 1    |    |    |    |    |    |    |    |
| 2    |    |    |    |    |    |    |    |
| ...  |    |    |    |    |    |    |    |


## A reasonable separation distance

So far, the "tolerance radius" we previously talked about has not been considered.
This is actually what the *great-circle* algorithm -- next to be presented --
adds to the discussion.
Nevertheless we here consider it as a post-processing step to prepare our next discussion.

On not considering a reasonable separation distance, we might end up with matches between
("nearest-neighbors") objects that are not actually close to each other.
In our current example, it can be argued that among all the possible matches `Table 3`_
only *[A:4,B:2]* and *[A:12,B:5]* do represent real matches -- if we consider ``30 pixels``
to be sufficiently close.

Let us say the objects in our tables to have an "error radius" of ``20 pixels``.
Considering the posistional errors between "A" and "B" to be independent,
their combined errors is given by a sum in quadrature,

$$
   \epsilon &= \sqrt{err_A^2 + err_B^2} \\
            &= \sqrt{20^2 + 20^2} \\
            &= 28.2
$$
, which gives us a tolerance radius of $\delta_r = 28$ pixels.

Using this value as a threshold for the separation distance, our final table of matches
is left with the pair *[A:4,B:2]*.
The final matching table is given by `Table 4`_ below.

Table 4: Final matching-table

A_ID | x | y | B_ID | x | y | dist
--- | --- | --- | --- | --- | --- | ---
1 | 301 | 34 | nan | nan | nan | nan
2 | 51 | 74 | nan | nan | nan | nan
3 | 230 | 145 | nan | nan | nan | nan
4 | 404 | 232 | 2 | 391 | 246 | 19
5 | 229 | 265 | nan | nan | nan | nan
6 | 52 | 286 | nan | nan | nan | nan
7 | 83 | 288 | nan | nan | nan | nan
8 | 346 | 289 | nan | nan | nan | nan
9 | 317 | 339 | nan | nan | nan | nan
10 | 214 | 376 | nan | nan | nan | nan
11 | 465 | 388 | nan | nan | nan | nan
12 | 401 | 455 | nan | nan | nan | nan
