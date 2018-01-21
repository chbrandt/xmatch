# xmatch

Implementation of cross-matching algorithms for astronomical catalogs.

Cross-match is a daily procedure in astrophysics, and although it's been around
virtually since ever it is a non-trivial hot topic when multiwavelength catalogs
are in hands.

Cross-matching has on its basis the relative position of objects in the projected
sky: the same object seeing by two different instruments should have the same
(or very close) registered positions in the sky. This is an absolute true when
if/when the object emits light similarly in the frequencies both instruments
observe, but may get a bit non-linear when the object differently in different
wavebands or the instruments observe in quite different resolutions.

The `xmatch` tool is python interface to the cross-matching procedure, implmenting
the most commong algorithms in use:
* bidirectional nearest-neighbour (nn)
* great-circle (gc): nn with maximum distance threshold
* maximum likelihood estimator (mle): a feature column is used to evaluate multiple
counterpart candidates


## Install

```bash
# python setup.py install
```
or
```bash
# pip install .
```


## Using it

Consider we have two catalogs "A" and "B" with at least the columns `ra`, `dec`,
`id` (for *right ascension*, *declination*, *identifier*, resp.) that we want to
cross-match.

If all we want is to "find the nearest object in catalog 'B' (aka, counterpart)
for each object from 'A' (target), if it is within a limiting 'radius'; if it is
not, say none was found", then we run `xmatch` on its default:
```python
>>> radius = 5 # arcsec
>>> from xmatch import xmatch
>>> out_x = xmatch(catalog_A[['ra','dec','id']],
                   catalog_B[['ra','dec','id']],
                   radius=radius)
```
`out_x` provides a table containing all entries from `catalog_A` and the
corresponding counterpart -- *if any* -- from `catalog_B`.


## Examples

Check out `xmatch/docs/notebooks` for some examples.

-- /.\
