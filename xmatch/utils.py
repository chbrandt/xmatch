# -*- coding:utf-8 -*-
import logging

def stats(arrei):
    '''
    Outputs in a dictionary:
    - min
    - max
    - mean
    - std
    - quantiles (1,2,3)
    '''
    if not len(arrei):
        logging.error("There is no data to comute stats; given array is empty.")
        return None

    from collections import OrderedDict
    sts = OrderedDict()

    sts['length'] = len(arrei)

    if all([ hasattr(arrei,attr) for attr in ['min','max','mean','std'] ]):
        sts['min'] = arrei.min()
        sts['max'] = arrei.max()
        sts['mean'] = arrei.mean()
        sts['std'] = arrei.std()
    else:
        from numpy import amin, amax, mean, std
        sts['min'] = amin(arrei)
        sts['max'] = amax(arrei)
        sts['mean'] = mean(arrei)
        sts['std'] = std(arrei)

    from numpy import percentile
    _q1 = percentile(arrei, 25)
    _q2 = percentile(arrei, 50)
    _q3 = percentile(arrei, 75)
    sts['25%'] = _q1
    sts['50%'] = _q2
    sts['75%'] = _q3

    return sts
