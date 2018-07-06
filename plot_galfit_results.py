#!/usr/bin/env python3


import os
import sys
import pyfits
import numpy
import matplotlib
from PIL import Image


def scale_and_normalize(data, scaling_mode, flux_max):

    if (scaling_mode == "asinh" or True):
        return numpy.arcsinh((data + 0.2*flux_max)/(1.2*flux_max))


def plot_galfit_result(fits_fn, plot_fn, scaling_mode='asinh'):

    hdulist = pyfits.open(fits_fn)

    img_input = hdulist[1].data
    img_fit = hdulist[2].data
    img_residuals = hdulist[3].data

    # print(img_input.shape, img_fit.shape, img_residuals.shape)

    flux_max = numpy.max(img_fit)
    # print(flux_max)

    output_array = numpy.zeros((img_input.shape[0], 3*img_input.shape[1]+20))

    sx = img_input.shape[1]
    output_array[:, 0:sx]          = scale_and_normalize(img_input, scaling_mode, flux_max)
    output_array[:, sx+10:2*sx+10] = scale_and_normalize(img_fit, scaling_mode, flux_max)
    output_array[:, 2*sx+20:]      = scale_and_normalize(img_residuals, scaling_mode, flux_max)
    output_array[output_array <= 0] = 0.
    output_array[output_array >= 1] = 1.0

    img = Image.fromarray(numpy.uint8(output_array *255))
    img.transpose(Image.FLIP_TOP_BOTTOM).save(plot_fn)
    print("plot saves as %s" % (plot_fn))

if __name__ == "__main__":

    for fn in sys.argv[1:]:
        print(fn)
        out_fn = fn[:-5]+".png"

        plot_galfit_result(fn, out_fn)