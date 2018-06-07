#!/usr/bin/env python3

import os
import sys
import numpy
import pyfits
import scipy.ndimage.filters

sexcolumns = {
    'ra': 0,
    'dec': 1,
    'fwhm': 9,
    'x': 3,
    'y': 4,
    'id': 2,
}


if __name__ == "__main__":

    image_fn = sys.argv[1]
    catalog_fn = sys.argv[2]
    segmentation_fn = sys.argv[3]
    background_fn = sys.argv[4]
    variance_fn = sys.argv[5]

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

    background_hdu = pyfits.open(background_fn)
    background = background_hdu[0].data
    variance_hdu = pyfits.open(variance_fn)
    variance = variance_hdu[0].data

    # figure out which sources we dont want
    fwhm = catalog[:, sexcolumns['fwhm']]
    bad = fwhm < 6
    bad_sources = catalog[:,sexcolumns['id']][bad].astype(numpy.int)
    bad_catalog = catalog[bad]

    print(bad_sources)

    reuse = False
    boxsize=100
    if (not reuse):
        # find all bad pixels
        bad_pixels = numpy.zeros(segmenation.shape, dtype=numpy.bool)
        print(bad_pixels[:2,:3])
        for i, bad_source_id in enumerate(bad_sources):

            x = int(bad_catalog[i, sexcolumns['x']])
            y = int(bad_catalog[i, sexcolumns['y']])
            x1 = int(numpy.max([0, x-boxsize]))
            x2 = int(numpy.min([x+boxsize, img.shape[1]]))
            y1 = int(numpy.max([0, y-boxsize]))
            y2 = int(numpy.min([y+boxsize, img.shape[0]]))

            print("\rWorking on source %d of %d" % (i+1, bad_sources.shape[0]), end='', flush=True)
            bad_pixels_in_this_source = (segmenation[y1:y2, x1:x2] == bad_source_id)
            # pyfits.PrimaryHDU(data=bad_pixels_in_this_source.astype(numpy.int)).writeto("bad_mask_%d.fits" % (bad_source_id), clobber=True)

            bad_pixels[y1:y2, x1:x2] |= bad_pixels_in_this_source
            pass

        print(" ... all sources masked")
        pyfits.PrimaryHDU(data=bad_pixels.astype(numpy.int)).writeto("bad_mask_combined.fits", clobber=True)

    else:
        hdu = pyfits.open("bad_mask_combined.fits")
        bad_pixels = hdu[0].data.astype(numpy.bool)

    masked_image = img.copy()
    masked_image[bad_pixels] = numpy.NaN

    pyfits.PrimaryHDU(
        data=masked_image,
        header=img_hdu[0].header).writeto("image_masked.fits", clobber=True)


    # grow mask
    mask_grown = scipy.ndimage.filters.convolve(
                            input=bad_pixels.astype(numpy.int),
                            weights=numpy.ones((6,6)),
                            output=None,
                            mode='constant', cval=0.0).astype(numpy.bool)

    pyfits.PrimaryHDU(data=mask_grown.astype(numpy.int)).writeto("bad_mask_grown.fits", clobber=True)
    masked_image[mask_grown] = numpy.NaN
    pyfits.PrimaryHDU(
        data=masked_image,
        header=img_hdu[0].header).writeto("image_masked_grown.fits", clobber=True)

    masked_image[mask_grown] = background[mask_grown]
    pyfits.PrimaryHDU(
        data=masked_image,
        header=img_hdu[0].header).writeto("image_masked_backgroundfixed.fits", clobber=True)

    # noise_image = numpy.random.randn(img.shape[0], img.shape[1]) * 0.0166987
    noise_image = numpy.random.randn(img.shape[0], img.shape[1]) * numpy.sqrt(variance)
    masked_image[mask_grown] += noise_image[mask_grown]
    pyfits.PrimaryHDU(
        data=masked_image,
        header=img_hdu[0].header).writeto("image_masked_backgroundnoisy.fits", clobber=True)

