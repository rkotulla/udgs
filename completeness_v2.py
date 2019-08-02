#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import time

import astropy.io.fits as pyfits
import numpy
import distutils.spawn



def range_to_list(arg):

    try:
        val = float(arg)
        return numpy.array([val])
    except ValueError:
        # this means it's not a simple number
        pass

    # check if it in the format a..b:s
    if (arg.find("..")>0 and arg.find(":")):
        # we got this format
        try:
            a = float(arg.split("..")[0])
            b = float(arg.split("..")[1].split(":")[0])
            s = float(arg.split(":")[1])
            # print(a,b,s)
            return numpy.arange(a, b+s, s)
        except ValueError:
            print("Unable to understand format of %s" % (arg))
            return None
    elif (arg.find(",")):
        # this is a list of values
        try:
            l = [float(f) for f in arg.split(",")]
            return numpy.array(l)
        except ValueError:
            print("Illegal format or value in %s" % (arg))

    return None


if __name__ == "__main__":

    # get path of this executable
    fn = os.path.abspath(__file__)
    dirname,_ = os.path.split(fn)
    config_dir = os.path.join(dirname, "config")

    # setup command line parameters
    cmdline = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    cmdline.add_argument("--conf", dest="sex_conf", default=os.path.join(config_dir, "sex4psfex.conf"),
                         help="source extractor config filename")
    cmdline.add_argument("--params", dest="sex_params", default=os.path.join(config_dir, "sex4psfex.param"),
                         help="number of models to insert into each image")
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("--sex", dest="sex_exe", default=distutils.spawn.find_executable("sex"),
                         help="location of SExtractor executable")
    cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
                         help="weight map")
    cmdline.add_argument("--psf", dest='psf_image', type=str, default="_image.fits:_psf.fits",
                         help="psf basename")

    cmdline.add_argument("--refcat", dest='ref_cat', default="_image.fits:_image.vot",
                         help="reference catalog")


    cmdline.add_argument("--mag", dest='mag', default="18..20:0.5",
                         help="magnitude-range, from-to:step")
    cmdline.add_argument("--re", dest='r_eff', default="10..50:5",
                         help="effective radius range, from-to:step")
    cmdline.add_argument("--ar", dest='axisratio', default="0.2..1:0.2",
                         help="axis ratio range, from-to:step")
    cmdline.add_argument("--sersic", dest='sersic', default="0.5..4:0.5",
                         help="sersic index range, from-to:step")
    cmdline.add_argument("--pa", dest='posangle', default="0..180:30",
                         help="position angle range, from-to:step")

    cmdline.add_argument("-n", dest='n_models', default=10, type=int,
                         help="number of models for each configuration")


    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    # cmdline.print_help()
    args = cmdline.parse_args()

    construct_weight_fn = False
    if (args.weight_image is not None):
        construct_weight_fn = (args.weight_image.find(":") > 0)

    # First, lets work out all the parameters to iterate over to get an
    # idea what models we need

    print("mags:       ", range_to_list(args.mag))
    print("eff radius: ", range_to_list(args.r_eff))
    print("b/a:        ", range_to_list(args.axisratio))
    print("sersic:     ", range_to_list(args.sersic))
    print("pos-angle:  ", range_to_list(args.posangle))

    file_queue = multiprocessing.JoinableQueue()

    # # feed with files to sextract
    # for fn in args.input_images:
    #     if (os.path.isfile(fn)):
    #
    #         if (construct_weight_fn):
    #             _parts = args.weight_image.split(":")
    #             search = _parts[0]
    #             replace = _parts[1]
    #             weight_fn = fn.replace(search, replace)
    #         else:
    #             weight_fn = args.weight_image
    #
    #         file_queue.put((fn, weight_fn))
    #
    # # insert termination commands
    # for i in range(args.number_processes):
    #     file_queue.put((None))
    #
    # # start all processes
    # processes = []
    # for i in range(args.number_processes):
    #     p = multiprocessing.Process(
    #         target=run_sex_psfex,
    #         kwargs=dict(
    #             file_queue=file_queue,
    #             sex_exe=args.sex_exe,
    #             sex_conf=args.sex_conf, sex_param=args.sex_params,
    #             psfex_conf=args.psfex_conf, psfex_exe=args.psfex_exe,
    #             supersample=args.supersample,
    #         )
    #     )
    #     p.daemon = True
    #     p.start()
    #     processes.append(p)

    # now wait for all work to be done
    file_queue.join()


