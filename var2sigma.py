#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy

if __name__ == "__main__":

    var_fn = sys.argv[1]

    if (var_fn.endswith("_var.fits")):
        sigma_fn = var_fn[:-9]+"_sigma.fits"
    else:
        sigma_fn = sys.argv[2]

    hdulist = pyfits.open(var_fn)
    variance =  hdulist[0].data

    sigma = numpy.sqrt(variance)

    hdulist[0].data = sigma

    hdulist.writeto(sigma_fn, clobber=True)
