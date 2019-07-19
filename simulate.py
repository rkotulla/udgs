#!/usr/bin/env python3

import os, sys, pyfits, numpy
from astLib import astWCS
import scipy.ndimage
import tempfile


tmpdir = '/tmp'
executable = "/work/sersic/makesersic_v6"

def create_model(imgsize=1000,
                 cx=500, cy=500,
                 sersic_n=1.0,
                 inst_magnitude=-7.,
                 halflight_semimajor_axis=80,
                 axis_ratio=0.9,
                 position_angle=30,
                 background=0,
                 cleanup=True
                 ):

    gauss_sigma = 0.
    gauss_n_sigma = 0.
    noise_level = 0.
    subsampling = 1

    # create some unique temporary filename
    _id, filename = tempfile.mkstemp(suffix="",
                     prefix="sersic_",
                     dir=tmpdir,
                     text=True)
    print(filename)

    conf_filename = filename+".conf"
    fits_filename = filename+".fits"

    with open(conf_filename, "w") as fh:


        conf_text = """\
%(imgsize)f # image size
%(cx)f #  center-x
%(cy)f #  center_y
%(sersic_n)f # sersic n
%(app_magnitude)f 1 # apparent magnitude (1) or I_at_Reff (0) in counts per sq.pixel
%(halflight_semimajor_axis)f  # Half-light semi-major axis (pixels)
%(axis_ratio)f             # Axis ratio
%(position_angle)f           # Position angle (degrees East from North)
0               # B4 !experimental! - Deviation from elliptical shape, i.e. "boxiness" - try 0.02 for boxy and -0.02 for disky - NOTE: seems to be some kind of bug producing a few bright pixels...use with care
%(background)f               # Background level
%(gauss_sigma)f             # PSF Gaussian sigma in pixels (Note: divide Gaussian FWHM by 2.3548 to get Gaussian sigma)
%(gauss_n_sigma)f                # Number of sigmas out to where the PSF should be calculated
%(noise_level)f  234   # Noise level (Gaussian); if not zero, also give some positive integer as random seed
%(subsampling)d                # Increase the default subpixel sampling resolution by this factor (values >1 mean finer sampling) -- expect computation time to increase by at least factor^3 !
""" % dict(
            imgsize=int(imgsize), cx=float(cx), cy=float(cy),
            sersic_n=float(sersic_n),
            app_magnitude=float(inst_magnitude),
            halflight_semimajor_axis=float(halflight_semimajor_axis),
            axis_ratio=float(axis_ratio), position_angle=float(position_angle),
            background=float(background), gauss_sigma=float(gauss_sigma), gauss_n_sigma=float(gauss_n_sigma),
            noise_level=float(noise_level), subsampling=int(subsampling)
        )
        fh.write(conf_text)
        print(conf_text)

    #
    # Now run the tool to create the FITS file from the config file
    #
    cmd = "%s %s %s" % (executable, conf_filename, fits_filename)
    # TODO: Make better call than using os.system
    os.system(cmd)

    #
    # Now we have a FITS file with the image in it
    #
    hdulist = pyfits.open(fits_filename)
    model_data = hdulist[0].data

    if (cleanup):
        os.remove(conf_filename)
        os.remove(fits_filename)

    return model_data


    # #
    # # Compute the model file using Thorsten's tool
    # #
    # model_conf_fn = (basename+".model.conf") % (_fn+1)
    # model_fits_fn = (basename+".model.fits") % (_fn+1)
    # with open(model_conf_fn_in, "r") as mcf:
    #
    #     lines = mcf.readlines()
    #     lines = [x.strip() for x in lines]
    #
    #     line_numbers = [0,1,2,5]
    #     line_formats = ["%d", "%.1f", "%.1f", "%f"]
    #
    #     for i in range(len(line_numbers)): #[0,1,2,5]:
    #         line_no = line_numbers[i]
    #         line_fmt = line_formats[i]
    #         _size = lines[line_no].split()
    #         size_arcsec = float(_size[0]) * 0.25
    #         size_odiraw = size_arcsec / pixelscale
    #         if (size_odiraw > 1020):
    #             size_odiraw = 1020.
    #         _size[0] = line_fmt % (size_odiraw)
    #         lines[line_no] = " ".join(_size)
    #
    #     with open(model_conf_fn, "w") as f:
    #         f.write("\n".join(lines))
    #
    # cmd = "./makesersic_v6 %s %s" % (model_conf_fn, model_fits_fn)
    # print cmd
    # os.system(cmd)
    # # model_fn = "xxx.fits"
    #
    # model_hdu = pyfits.open(model_fits_fn)
    # model = model_hdu[0].data
    #
    # mx = (model.shape[1]-1)/2
    # my = (model.shape[0]-1)/2
    #
    # zp = 20.
    # if ("PHOTZP_X" in hdulist[0].header):
    #     zp = hdulist[0].header["PHOTZP_X"]
    # elif ("MAGZERO" in hdulist[0].header):
    #     # most likely a stack
    #     zp = hdulist[0].header['MAGZERO']
    # print "Using zeropoint of %.3f" % (zp)
    # delta_zp = zp - 26.
    # zp_scaling = math.pow(10., 0.4*delta_zp)
    #
    # print "\n\nAdding %s to %s" % (model_fits_fn, infits)
    #
    # print "Convolving with seeing:"
    # if ('SEEING' in hdulist[0].header):
    #     seeing = hdulist[0].header['SEEING']
    # else:
    #     seeing = 1.2
    #
    # gauss_sigma = seeing / 2.3548
    # smoothed_model = scipy.ndimage.filters.gaussian_filter(
    #     input=model, sigma=gauss_sigma,
    #     order=0,
    #     output=None,
    #     mode='reflect')
    #
    # for ext in hdulist:
    #
    #     if (not is_image_extension(ext)):
    #         continue
    #
    #     wcs = astWCS.WCS(ext.header, mode='pyfits')
    #     data = ext.data
    #
    #     # check if cutout is in this OTA
    #     xy = numpy.array(wcs.wcs2pix(ra, dec))
    #     if (xy[0] > -mx and
    #         xy[0] < data.shape[1]+mx  and
    #         xy[1] > -my and
    #         xy[1] < data.shape[0]+my):
    #     # if (xy[0] > -model.shape[1] and
    #     #     xy[0] < data.shape[1]+model.shape[1]  and
    #     #     xy[1] > -model.shape[0] and
    #     #     xy[1] < data.shape[0]+model.shape[0]):
    #
    #         print xy, data.shape, model.shape
    #
    #         # the model covers at least part of this image
    #         _left = int(xy[0] - mx)
    #         _right = _left + model.shape[1]
    #         _bottom = int(xy[1] - my)
    #         _top = _bottom + model.shape[0]
    #         print ext.name, _left, _right, _top, _bottom
    #
    #         #
    #         # handle edge cases
    #         #
    #         _mx1, _my1, _mx2, _my2 = 0,0,model.shape[1], model.shape[0]
    #
    #         if (_left < 0):
    #             _mx1 = -1*_left
    #             _left = 0
    #         if (_right > data.shape[1]):
    #             _mx2 = data.shape[1] - _left
    #             _right = data.shape[1]
    #         if (_bottom  < 0):
    #             _my1 = -1*_bottom
    #             _bottom = 0
    #         if (_top > data.shape[0]):
    #             _my2 = data.shape[0] - _bottom
    #             _top = data.shape[0]
    #
    #         #
    #         # Now add galaxy to the image
    #         #
    #         print _mx1,_mx2,_my1,_my2," ==>> ", _left,_right,_bottom,_top
    #
    #         ext.data[_bottom:_top, _left:_right] += (user_scale*zp_scaling)*smoothed_model[_my1:_my2, _mx1:_mx2]
    #         out_hdulist.append(pyfits.ImageHDU(data=ext.data, header=ext.header))
    #
    # out_hdulist = pyfits.HDUList(out_hdulist)
    # out_fn = (basename+".fits") % (_fn+1) #"model_added.%d.fits" % (_fn+1)
    # clobberfile(out_fn)
    # print "writing %s" % (out_fn)
    # out_hdulist.writeto(out_fn, clobber=True)
    #



if __name__ == "__main__":

    model = create_model()
    print(model.shape)
