#!/usr/bin/env python
from __future__ import absolute_import, print_function
import sys

from xmatch import xmatch


def open_table(filename, sep=';', columns=None):
    import pandas
    df = pandas.read_csv(filename, sep=sep, comment='#')
    if columns is not None and len(columns) > 0:
        assert isinstance(columns, (list, tuple))
        assert all(col in df.columns for col in columns)
        df = df[columns]
    return df


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Cross-Match catalogs")

    parser.add_argument('catalog1', type=str,
                        help="Target catalog (csv filename)")
    parser.add_argument('--ra1', type=str, default='RA',
                        help="RA column name in 'cat1'")
    parser.add_argument('--dec1', type=str, default='DEC',
                        help="Dec column name in 'cat1'")
    parser.add_argument('--id1', type=str, default=None,
                        help="ID (ObjID) column name in 'cat1'")
    parser.add_argument('--sep1', type=str, default=';',
                        help="Column separator in 'cat1'")

    parser.add_argument('catalog2', default=None, type=str,
                        help="Ancillary catalog (csv filename)")
    parser.add_argument('--ra2', type=str, default='RA',
                        help="RA column name in 'cat2'")
    parser.add_argument('--dec2', type=str, default='DEC',
                        help="Dec column name in 'cat2'")
    parser.add_argument('--id2', type=str, default=None,
                        help="ID (ObjID) column name in 'cat2'")
    parser.add_argument('--sep2', type=str, default=';',
                        help="Column separator in 'cat2'")

    parser.add_argument('-r', '--radius', type=float,
                        default=None,
                        help=("Radius (in 'arcsec' units) to search for ancillary objects."
                              "If 'None', Nearest-Neighbour method is used."))

    parser.add_argument('--method', default='gc', type=str,
                        choices=['gc','mle','nn','filter'],
                        help="Method for x-matching the catalogs")
    parser.add_argument('--filter-method', default=None, type=str,
                        help="Method for filtering catalogs")
    parser.add_argument('--mle-column', default=None, type=str,
                        help="Column to use from catalog2 for the likelihood")

    args = parser.parse_args()

    radius = args.radius
    if args.method == 'filter' or args.method == 'gc':
        assert radius != None and radius > 0, "Give me a radius value to search for target siblings"

    # Target catalog (cat1)
    cols1 = [args.ra1, args.dec1]
    if args.id1 is not None:
        cols1.append(args.id1)
    cat1 = open_table(args.catalog1, sep=args.sep1)#, columns=cols1)
    if args.id1 is None:
        cols1.append('ID')
        cat1.index.name = 'ID'
        cat1.reset_index(inplace=True)

    # Ancillary catalog (cat2)
    if args.catalog2:
        cols2 = dict(ra=args.ra2, dec=args.dec2)
        if args.id2 is not None:
            cols2['id'] = args.id2
        cat2 = open_table(args.catalog2, sep=args.sep2)#, columns=cols2)
        if args.id2 is None:
            cols2['id'] = 'ID'
            cat2.index.name = cols2['id']
            cat2.reset_index(inplace=True)
    else:
        cat2 = cat1
        cols2 = cols1

    from astropy.coordinates import SkyCoord
    from astropy import units
    coords = SkyCoord(cat1[cols1['ra']], cat1[cols1['dec']], unit=(units.deg,units.hourangle))
    cat1[cols1['ra']] = coords.icrs.ra
    cat1[cols1['dec']] = coords.icrs.dec
    del coords

    coords = SkyCoord(cat2[cols2['ra']], cat2[cols2['dec']], unit=(units.deg,units.hourangle))
    cat2[cols2['ra']] = coords.icrs.ra
    cat2[cols2['dec']] = coords.icrs.dec
    del coords

df['snr'] = df['nufnu_3keV(erg.s-1.cm-2)'] / df['nufnu_error_3keV(erg.s-1.cm-2)']

cols = dict(ra='RA', dec='DEC', id='OBJID')

from astropy.coordinates import Angle
rad = Angle(5,'arcsec')

from xmatch import xmatch
%time xcat = xmatch(df, df, cols, cols, radius=rad, snr_column='snr')
    print('Target catalog:\n', cat1.describe())#include='all'))
    print('Ancillary catalog:\n', cat2.describe())#include='all'))
    sys.exit()
    cat = xmatch(cat1, cols1,
                 cat2, cols2,
                 radius=rs, method=method)
