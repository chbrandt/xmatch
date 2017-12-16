# -*- coding=utf-8 -*-
import logging

import astropy

# def gc(A_coord, B_coord, radius, parallel=False, nprocs=None):
def gc(A_coord, B_coord, radius):
    '''
    Match catalogs by position within a radial distance

    Input:
     - radius : float or ~astropy.Quantity
            If a float is given it is assumed to be in 'arcsec'
    '''
    from astropy.units import Quantity
    if not isinstance(radius,Quantity):
    try:
        radius = float(radius)
        radius = Quantity(radius,'arcsec')
    except Exception as e:
        print("Error: Got exception {}\n".format(e))
        return None
    return _gc_serial(A_coord, B_coord, radius)


def _gc_serial(A_coord, B_coord, radius):

    from astropy.coordinates import search_around_sky

    def assert_input(A_coord, B_coord, radius):
        from astropy.coordinates import SkyCoord
        from astropy.units import Quantity
        assert isinstance(A_coord,SkyCoord), "Was expecting an ~astropy.coordinates.SkyCoord instance for 'A_coord'."
        assert isinstance(B_coord,SkyCoord), "Was expecting an ~astropy.coordinates.SkyCoord instance for 'B_coord'."
        assert isinstance(radius,Quantity), "Was expecting an ~astropy.units.Quantity instance for 'radius'"
    assert_input(A_coord, B_coord, radius)

    logging.info("Searching B_coord {1} objects, {0} neighbors.".format(len(B_coord),len(A_coord)))
    _prau = A_coord.ra.unit
    _pdecu = A_coord.dec.unit
    _nrau = B_coord.ra.unit
    _ndecu = B_coord.dec.unit
    logging.debug("Unit of coordinates being matched: ({0},{1}) and ({2},{3})".format(_prau,_pdecu,_nrau,_ndecu))

    match_A_gc_idx, match_B_gc_idx, match_gc_sep, _d3d = search_around_sky(A_coord, B_coord, radius)

    from .utils import stats
    _sts = stats.basic(match_gc_sep.value)
    logging.info("Basic stats of distances between matchings: {}".format(_sts))

    assert len(match_A_gc_idx) == len(match_B_gc_idx)

    return match_A_gc_idx, match_B_gc_idx, match_gc_sep
