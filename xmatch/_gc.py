# -*- coding=utf-8 -*-
import logging


def gc(A_coord, B_coord, radius):
    '''
    Match catalogs by position within a radial distance

    Input:
     - A_coord: astropy.coordinates.SkyCoord
     - B_coord: astropy.coordinates.SkyCoord
     - radius : astropy.coordinates.Angle

    Output:
     - matched_A_indexes: indexes of A entry that match the corresponding position in B
     - matched_B_indexes: indexes of B entry that match the corresponding position in A
     - separation_AB_val: separation between matched_{AB}_indexes, astropy.coordinates.Angle

    * All outputs (1D arrays) have the same length.
    '''
    from astropy.coordinates import Angle, SkyCoord
    assert isinstance(A_coord, SkyCoord), "Was expecting an ~astropy.coordinates.SkyCoord instance for 'A_coord'."
    assert isinstance(B_coord, SkyCoord), "Was expecting an ~astropy.coordinates.SkyCoord instance for 'B_coord'."
    assert isinstance(radius, Angle), "Was expecting an ~astropy.coordinates.Angle instance for 'radius'"

    return _gc_serial(A_coord, B_coord, radius)


def _gc_serial(A_coord, B_coord, radius):
    from astropy.coordinates import search_around_sky

    logging.info("Searching B_coord {1} objects, {0} neighbors.".format(len(B_coord),len(A_coord)))

    match_A_gc_idx, match_B_gc_idx, match_gc_sep, _d3d = search_around_sky(A_coord, B_coord, radius)
    assert len(match_A_gc_idx) == len(match_B_gc_idx)

    from .utils import stats
    _sts = stats(match_gc_sep.value)
    logging.info("Basic stats of distances between matchings: {}".format(_sts))

    return match_A_gc_idx, match_B_gc_idx, match_gc_sep
