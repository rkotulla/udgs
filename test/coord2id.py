#!/usr/bin/env python3

import os
import sys
import astropy.table

import ephem



def coord2id(ra, dec):

    e = ephem.Equatorial(ra, dec)
    print(e, e.ra, e.dec)
    return "xxx"


if __name__ == "__main__":

    cat_fn = sys.argv[1]

    cat = astropy.table.Table.read(cat_fn)

    for src in cat[:5]:
        id = coord2id(src['ALPHA_J2000'], src['DELTA_J2000'])
        print(id)