#!/usr/bin/env python3


import os
import sys
import numpy
import pyfits
import subprocess
import argparse

import astLib.astWCS as astWCS

import simulate


if __name__ == "__main__":

    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--nsims", dest="n_simulations",
                         default=1, type=int,
                         help="total number of simulated images (real+models) to generate from each input image")
    cmdline.add_argument("--perimage", dest="n_per_image",
                         default=1, type=int,
                         help="number of models to insert into each image")
    cmdline.add_argument("--mag", type=float, default=20.,
                         help="integrated magnitude of model")
    cmdline.add_argument("--axisratio", dest="axis_ratio", type=float, default=1.0,
                         help="axis ratio")
    cmdline.add_argument("--sersic", dest="sersic_index", type=float, default=1.0,
                         help="sersic index")
    cmdline.add_argument("--sma", dest="halflight_semimajor_axis", type=float, default=25.,
                         help="half-light semi-major axis [pixels]")
    cmdline.add_argument("--posangle", dest="position_angle", type=float, default=0.0,
                         help="position angle [degrees]")
    cmdline.add_argument("--outdir", dest="output_directory", type=str, default="sims/",
                         help="output directory to hold simulated images")
    cmdline.add_argument("--autodir", dest="auto_output_directory", default=False,
                         action='store_true',
                         help="auto-generate output directory to hold simulated images")
    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()


    #print(args)

    # make sure the output directory exists
    output_directory = args.output_directory
    if (args.auto_output_directory):
        # we automatically generate the output directory name
        _dir = "mag%05.2f_sersic%04.2f_sma%05.1f_axis%04.2f_pa%03.0f" % (
            args.mag, args.sersic_index, args.halflight_semimajor_axis,
            args.axis_ratio, args.position_angle)
        output_directory = os.path.join(args.output_directory, _dir)
    if (not os.path.isdir(output_directory)):
        print("Creating output directory: %s" % (output_directory))
        os.makedirs(output_directory)

        # write a quick summary file saving run paramters
        # for easier use during completeness checking
        numpy.savetxt(os.path.join(output_directory, ".params"),
                      numpy.array([args.n_per_image,
                                   args.mag,
                                   args.halflight_semimajor_axis,
                                   args.axis_ratio,
                                   args.position_angle,
                                   args.sersic_index,
                                   ]))

    #
    # Now do some magic: create model galaxies and insert into all input images
    # until we have inserted the requested number of models
    #

    for i_fn, input_fn in enumerate(args.input_images):

        print("Inserting models into %s" % (input_fn))
        in_hdulist = pyfits.open(input_fn)
        input_img = in_hdulist[0].data
        print(input_img.shape)

        # start WCS class to convert X/Y to ra/dec
        wcs = astWCS.WCS(in_hdulist[0].header, mode='pyfits')

        _, bn = os.path.split(input_fn)

        mag_zeropoint = 26.
        inst_magnitude = args.mag - mag_zeropoint
        print(inst_magnitude, args.mag)

        model = simulate.create_model(
            sersic_n=args.sersic_index,
            inst_magnitude=inst_magnitude,
            halflight_semimajor_axis=args.halflight_semimajor_axis,
            axis_ratio=args.axis_ratio,
            position_angle=0.,
            background=0.,
        )
        total_flux = numpy.sum(model)
        mag = -2.5*numpy.log10(total_flux)
        print(mag, mag+mag_zeropoint)

        for i_simulation in range(args.n_simulations):
            print(input_fn, i_simulation)

            sim_image = input_img.copy()

            # get random positions
            pos_xy = numpy.random.random((args.n_per_image, 2))
            pos_xy[:, 0] *= input_img.shape[1]
            pos_xy[:, 1] *= input_img.shape[0]
            pos_xy = pos_xy.astype(numpy.int)
            # print("\n\n\n",pos_xy)

            #
            # also convert all x/y to ra/dec for easier plotting
            #
            ra_dec = numpy.array(wcs.pix2wcs(pos_xy[:,0]+1, pos_xy[:,1]+1))

            sim_models = numpy.zeros_like(sim_image)

            for i_model in range(args.n_per_image):

                [cx, cy] = pos_xy[i_model]
                # print("\n",i_simulation, i_model, cx,cy)

                # Now add the model to the input data
                sx1, sx2, sy1, sy2 = 0, model.shape[1], 0, model.shape[0]

                tx1 = cx - model.shape[1]//2
                tx2 = cx + model.shape[1]//2
                ty1 = cy - model.shape[0]//2
                ty2 = cy + model.shape[1]//2
                if (tx1 < 0):
                    sx1 += -1*tx1
                    tx1=0
                if (tx2 > input_img.shape[1]):
                    sx2 -= (tx2 - input_img.shape[1])
                    tx2 = input_img.shape[1]
                if (ty1 < 0):
                    sy1 += -1*ty1
                    ty1=0
                if (ty2 > input_img.shape[0]):
                    sy2 -= (ty2 - input_img.shape[0])
                    ty2 = input_img.shape[0]

                #
                # tx1 = int(numpy.max([0, cx-model.shape[1]//2]))
                # tx2 = int(numpy.min([input_img.shape[1], cx+model.shape[1]//2]))
                # ty1 = int(numpy.max([0, cy-model.shape[0]//2]))
                # ty2 = int(numpy.min([input_img.shape[0], cy+model.shape[0]//2]))
                #
                # sx1 = model.shape[1]//2 + cx - (cx-tx1)
                #
                # # int(model.shape[1]//2-tx2+cx)
                # sx2 = sx1 + (tx2-tx1)
                # sy1 = int(model.shape[0]//2-ty2+cy)
                # sy2 = sy1 + (ty2-ty1)

                # print(tx1, tx2, ty1, ty2, " "*20, sx1, sx2, sy1, sy2)

                # Add the model to the input image
                sim_models[ty1:ty2, tx1:tx2] += model[sy1:sy2, sx1:sx2]


            # out_fn = "%s.%03d.mfits" % (bn[:-5], i_simulation+1)
            # output_img_fn = os.path.join(output_directory, out_fn)
            # print("Writing to %s" % (output_img_fn))
            #
            # out_hdulist = pyfits.PrimaryHDU(data=sim_models.astype(numpy.float32),
            #                                 header=in_hdulist[0].header)
            # out_hdulist.writeto(output_img_fn, clobber=True)


            sim_image += sim_models
            out_fn = "%s.%03d.fits" % (bn[:-5], i_simulation+1)
            output_img_fn = os.path.join(output_directory, out_fn)
            print("Writing to %s" % (output_img_fn))

            out_hdulist = pyfits.PrimaryHDU(data=sim_image.astype(numpy.float32),
                                            header=in_hdulist[0].header)
            out_hdulist.writeto(output_img_fn, clobber=True)

            # also create a log-file reporting all parameters of the inserted
            # model galaxies so we can check for detection efficiencies
            _log_fn = "%s.%03d.log" % (bn[:-5], i_simulation+1) #"sim_%03d.log" % (i_simulation+1)
            log_fn = os.path.join(output_directory, _log_fn)
            combined = numpy.array([
                numpy.arange(args.n_per_image),
                pos_xy[:,0]+1, pos_xy[:,1]+1,
                ra_dec[:,0], ra_dec[:,1],
                numpy.ones((args.n_per_image)) * args.mag,
                numpy.ones((args.n_per_image)) * args.axis_ratio,
                numpy.ones((args.n_per_image)) * args.halflight_semimajor_axis,
                numpy.ones((args.n_per_image)) * args.sersic_index,
                numpy.ones((args.n_per_image)) * args.position_angle,
            ]).T
            numpy.savetxt(log_fn, combined)
