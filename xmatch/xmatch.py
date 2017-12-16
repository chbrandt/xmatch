#!/usr/bin/env python
from .xmatchi import xmatch

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Cross-Match catalogs")

    parser.add_argument('--catA', type=str, required=True,
                        help="Target catalog (csv filename)")
    parser.add_argument('--raA', type=str, default='ra',
                        help="RA column name in 'catA'")
    parser.add_argument('--decA', type=str, default='dec',
                        help="Dec column name in 'catA'")

    parser.add_argument('--catB', type=str, required=True,
                        help="Ancillary catalog (csv filename)")
    parser.add_argument('--raB', type=str, default='ra',
                        help="RA column name in 'catB'")
    parser.add_argument('--decB', type=str, default='dec',
                        help="Dec column name in 'catB'")

    parser.add_argument('-r', '--search-radius', metavar='rs', type=float,
                        default=None,
                        help=("Radius (in 'arcsec' units) to search for ancillary objects."
                              "If 'None', Nearest-Neighbour method is used."))

    args = parser.parse_args()

    rs = args.rs
    if rs is None:
        method = 'nn'
    else:
        method = 'gc'
        from astropy import units
        rs = rs * units('arcsec')

    cat = xmatch(catA, columns_A,
                 catB, columns_B,
                 radius=rs, method=method)
