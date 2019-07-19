#!/usr/bin/env python3



import astropy.table
import astropy.io.fits as pyfits

import os
import sys


def read_definitions(definitions):

    # convert definitions into lookup index
    array_suffix = {}

    if (definitions is None):
        return array_suffix

    array_names = definitions.split("::")
    # print(array_names)
    for item_def in array_names:
        key = item_def.split(":")[0]
        values = item_def.split(":")[1].split(",")
        # print(key, values)
        array_suffix[key] = values

    return array_suffix



def fitsldac2vot(ldac_fn, vot_fn=None, array_suffix=None):

    hdu = pyfits.open(ldac_fn)
    ldac_data = hdu[2].data
    cat = astropy.table.Table(ldac_data)
    #print(cat.info())

    # print("===\n"*5)

    # now convert arrays to multiple columns
    for col in cat.colnames:

        # print(cat[col].shape, cat[col].ndim)

        if (cat[col].ndim == 2):
            for idx in range(cat[col].shape[1]):
                try:
                    suffix = array_suffix[col][idx]
                    new_colname = "%s_%s" % (col, suffix)
                except:
                    new_colname = "%s_%d" % (col, idx+1)
                # print(new_colname)

                cat[new_colname] = cat[col][:,idx]
            del(cat[col])

    # print(cat.colnames)
    if (vot_fn is not None):
        cat.write(vot_fn, format='votable', overwrite=True)

    return cat

if __name__ == "__main__":

    ldac_fn = sys.argv[1]
    vot_fn = sys.argv[2]

    definitions = "FLUX_RADIUS:50,80"
    array_suffix = read_definitions(definitions)

    cat = fitsldac2vot(ldac_fn=ldac_fn, vot_fn=vot_fn, array_suffix=array_suffix)
    print(cat)


