#!/usr/bin/env python3

import os
import sys
import pyfits
import numpy


if __name__ == "__main__":

    ref_file = sys.argv[1]
    ref_hdu = pyfits.open(ref_file)

    fits_sum = numpy.zeros_like(ref_hdu[0].data)

    for fn in sys.argv[2:]:

        print(fn)
        hdu = pyfits.open(fn)
        model = hdu[2].data

        fits_section = hdu[2].header['FITSECT'][1:-1]

        x1 = int(fits_section.split(",")[0].split(":")[0])
        x2 = int(fits_section.split(",")[0].split(":")[1])
        y1 = int(fits_section.split(",")[1].split(":")[0])
        y2 = int(fits_section.split(",")[1].split(":")[1])

        bg = float(hdu[2].header['2_SKY'].split(" ")[0])

        fits_sum[y1-1:y2, x1-1:x2] += (model - bg)

    ref_hdu[0].data = fits_sum
    ref_hdu.writeto("galfit_model_sum.fits", clobber=True)