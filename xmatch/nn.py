# -*- coding=utf-8 -*-
import logging

import astropy
import numpy

def nn(A_coord, B_coord, parallel=False, nprocs=None):
    """
    Nearest-Neighbor search

    Input:
     - A_coord : ~astropy.coordinates.SkyCoord
            reference catalog (catalog "A")
     - B_coord : ~astropy.coordinates.SkyCoord
            matching catalog (catalog "B")

    Output:
     - tuple with ~numpy.ndarray , ~astropy.units.Quantity
            array of respective (to 'A_coord') index entries in 'B_coord'
            , array of respective pair distances
    """
    if parallel and nprocs > 0:
        from .utils import parallel
        dview = parallel.setup(nprocs)
        if not dview:
            return None

        logging.debug("Running in parallel with {} processors.".format(len(dview)))
        match_A_nn_idx, match_A_nn_sep = _nn_parallel(A_coord, B_coord, dview=dview)
    else:
        logging.debug("Running in serial mode.")
        match_A_nn_idx, match_A_nn_sep = _nn_serial(A_coord, B_coord)

    return match_A_nn_idx, match_A_nn_sep


def _nn_serial(A_coord, B_coord):

    from astropy.coordinates import SkyCoord
    assert isinstance(A_coord,SkyCoord), "Was expecting a ~astropy.coordinates.SkyCoord instance."
    assert isinstance(B_coord,SkyCoord), "Was expecting a ~astropy.coordinates.SkyCoord instance."

    logging.info("Searching among {0} neighbors, {1} reference objects.".format(len(B_coord),len(A_coord)))
    _prau = A_coord.ra.unit
    _pdecu = A_coord.dec.unit
    _nrau = B_coord.ra.unit
    _ndecu = B_coord.dec.unit
    logging.debug("Unit of coordinates being matched: ({0},{1}) and ({2},{3})".format(_prau,_pdecu,_nrau,_ndecu))

    from astropy.coordinates import match_coordinates_sky
    match_A_nn_idx, match_A_nn_sep, _d3d = match_coordinates_sky(A_coord,B_coord)

    from .utils import stats
    _sts = stats.basic(match_A_nn_sep.value)
    logging.info("Basic stats of distances between matchings: {}".format(_sts))

    assert len(match_A_nn_idx) == len(A_coord)
    assert match_A_nn_idx.max() < len(B_coord)

    return (match_A_nn_idx, match_A_nn_sep)


def _nn_parallel(A_coord,B_coord,dview=None):

    assert dview, "A cluster clients hub, ie. 'dview', must be given."

    # Encapsulate some variables to send for processing
    def make_nn_search_parallel(foo,cat2):
        def pkg_nn_search(cat1,foo=foo,cat2=cat2):
            return foo(cat1,cat2)
        return pkg_nn_search
    # ---

    # Split array (of coordinates) in N pieces
    def split_array(A_coord,N):
        from numpy import arange,array_split
        index = arange(len(A_coord))
        A_pieces = [ A_coord[idx]   for idx in array_split( index,N ) ]
        return A_pieces

    # Join array/list of tuples in N pieces
    def join_array(A_outs):
        from numpy import append
        match_A_nn_idx = None
        match_A_nn_sep = None
        for each_out in A_outs:
            match_idx, match_sep = each_out
            if match_A_nn_idx is None:
                assert match_A_nn_sep is None
                match_A_nn_idx = match_idx
                match_A_nn_sep = match_sep
            else:
                match_A_nn_idx = append(match_A_nn_idx,match_idx)
                match_A_nn_sep = append(match_A_nn_sep,match_sep)
        return (match_A_nn_idx,match_A_nn_sep)
    # ---

    # A-B
    foo_match_coordinates = make_nn_search_parallel(nn_serial, B_coord)

    A_pieces = split_array(A_coord,N=len(dview))

    A_outs = dview.map_sync( foo_match_coordinates, A_pieces )

    # This is a hack to recompose the 'unit' parameter from 'match_sep' below;
    unit_sep = A_outs[-1][-1].unit
    #
    match_A_nn_idx,match_A_nn_sep = join_array(A_outs)
    #
    # Do the hack: ('match_sep' loose its unit during numpy.append in 'join')
    from astropy.coordinates import Angle
    match_A_nn_sep = Angle(match_A_nn_sep.value, unit=unit_sep)
    #

    return (match_A_nn_idx,match_A_nn_sep)
