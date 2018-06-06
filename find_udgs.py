#!/usr/bin/env python3

import os
import sys
import numpy
import pyfits



if __name__ == "__main__":

    image_fn = sys.argv[1]
    catalog_fn = sys.argv[2]
    segmentation_fn = sys.argv[3]

    # load catalog
    catalog = numpy.loadtxt(catalog_fn)
    print(catalog.shape)

    # load images
    img_hdu = pyfits.open(image_fn)
    img = img_hdu[0].data
    print(img.shape)

    segmentation_hdu = pyfits.open(segmentation_fn)
    segmenation = segmentation_hdu[0].data
    print(segmenation.shape)

    # figure out which sources we dont want
    fwhm = catalog[:, 9]
    bad = fwhm < 6
    bad_sources = catalog[:,2][bad].astype(numpy.int)
    print(bad_sources)

    # find all bad pixels
    bad_pixels = numpy.zeros(segmenation.shape, dtype=numpy.bool)
    print(bad_pixels[:2,:3])
    for i, bad_source_id in enumerate(bad_sources):

        print("Working on source %d of %d" % (i+1, bad_sources.shape[0]))
        bad_pixels_in_this_source = (segmenation == bad_source_id)
        # pyfits.PrimaryHDU(data=bad_pixels_in_this_source.astype(numpy.int)).writeto("bad_mask_%d.fits" % (bad_source_id), clobber=True)

        bad_pixels |= bad_pixels_in_this_source
        pass

    pyfits.PrimaryHDU(data=bad_pixels.astype(numpy.int)).writeto("bad_mask_combined.fits", clobber=True)

    masked_image = img.copy()
    masked_image[bad_pixels] = numpy.NaN

    pyfits.PrimaryHDU(
        data=masked_image,
        header=img_hdu[0].header).writeto("image_masked.fits", clobber=True)

