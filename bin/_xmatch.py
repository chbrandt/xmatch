#!/usr/bin/env python
from __future__ import absolute_import, print_function

from xmatch import xmatch


def open_table(filename, sep=';', columns=None):
    import pandas
    df = pandas.read_csv(filename, sep=sep)
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
    parser.add_argument('--ra1', type=str, default='ra', required=True,
                        help="RA column name in 'cat1'")
    parser.add_argument('--dec1', type=str, default='dec', required=True,
                        help="Dec column name in 'cat1'")
    parser.add_argument('--id1', type=str, default=None,
                        help="ID (ObjID) column name in 'cat1'")
    parser.add_argument('--sep1', type=str, default=';',
                        help="Column separator in 'cat1'")

    parser.add_argument('catalog2', default=None, type=str,
                        help="Ancillary catalog (csv filename)")
    parser.add_argument('--ra2', type=str, default='ra', required=True,
                        help="RA column name in 'cat2'")
    parser.add_argument('--dec2', type=str, default='dec', required=True,
                        help="Dec column name in 'cat2'")
    parser.add_argument('--id2', type=str, default=None,
                        help="ID (ObjID) column name in 'cat2'")
    parser.add_argument('--sep2', type=str, default=';',
                        help="Column separator in 'cat2'")

    parser.add_argument('-r', '--radius', metavar='rs', type=float,
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

    rs = args.rs
    if rs is None:
        method = 'nn'
    else:
        method = 'gc'
        from astropy import units
        rs = rs * units('arcsec')

    # Target catalog (cat1)
    cols1 = [args.ra1, args.dec1]
    if args.id1 is not None:
        cols1.append(args.id1)
    cat1 = open_table(args.cat1, sep=args.sep1, columns=cols1)
    if args.id1 is None:
        cols1.append('ID')
        cat1.index.name = 'ID'
        cat1.reset_index(inplace=True)

    # Ancillary catalog (cat2)
    cols2 = [args.ra2, args.dec2]
    if args.id2 is not None:
        cols2.append(args.id2)
    cat2 = open_table(args.cat2, sep=args.sep2, columns=cols2)
    if args.id2 is None:
        cols2.append('ID')
        cat2.index.name = 'ID'
        cat2.reset_index(inplace=True)


    cat = xmatch(cat1, cols1,
                 cat2, cols2,
                 radius=rs, method=method)
