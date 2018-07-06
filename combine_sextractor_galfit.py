#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
import multiprocessing
import pyfits



def read_results(hdr, component, parameter, keyname=None):

    value, uncert, flag = numpy.NaN, numpy.NaN, 99

    if (keyname is None):
        keyname = "%d_%s" % (component, parameter)
    try:
        result = hdr[keyname]
    except:
        return [value, uncert, flag]

    # try and see if we have a +/- result as well

    if (type(result) is str):
        _split = result.split("+/-")
        if (len(_split) == 2):
            # yes, there is a +/-
            # check if there is a * meaning this number should be marked as problematic
            if (_split[0].strip()[0] == "*"):
                flag = 2 # problematic
                value = float(_split[0].strip()[1:-1])
                uncert = float(_split[1].strip()[1:-1])
            else:
                flag = 0
                value = float(_split[0])
                uncert = float(_split[1])
        elif (result.startswith("[")):
            value = float(result.split("[")[1].split("]")[0])
            uncert = numpy.NaN
            flag = 1
        else:
            print("Unable to understand Galfit result: %s" % (result))
    else:
        try:
            value = float(result)
        except:
            value = numpy.NaN

    return [value, uncert, flag]



def parallel_combine(catalog_queue, galfit_directory='galfit'):

    # print("Hello from worker")
    while (True):
        cmd = catalog_queue.get()
        if (cmd is None):
            catalog_queue.task_done()
            # print("Shutting down")
            break

        cat_fn, udg_cat, combined_cat = cmd
        print(udg_cat)

        try:
            catalog = numpy.loadtxt(udg_cat)
        except:
            print("Error opening %s" % (udg_cat))
            catalog_queue.task_done()
            continue

        # print("xxx")
        if (catalog.ndim < 2 or catalog.shape[0] <= 0):
            print("Error with catalog %s" % (udg_cat))
            catalog_queue.task_done()
            continue

        with open(udg_cat, "r") as  cf:
            lines = cf.readlines()
            catalog_header = []
            for l in lines:
                if (l.startswith("#")):
                    catalog_header.append(l.strip())
            print(catalog_header)

        #
        # Now we can start the actual work
        #

        _dir, _fn = os.path.split(udg_cat)
        basename, _ = os.path.splitext(_fn)
        galfit_dir = os.path.join(_dir, galfit_directory)
        # print(galfit_dir)

        galfit_data = [None] * catalog.shape[0]
        for i_src, src in enumerate(catalog):

            src_id = int(src[2])

            galfit_fn = "%s.%05d.galfit.fits" % (basename, src_id)
            galfit_fullfn = os.path.join(galfit_dir, galfit_fn)
            print(galfit_fullfn)

            if (not os.path.isfile(galfit_fullfn)):
                print("File does not exist - failed fit? (%s)" % (galfit_fullfn))
                continue

            # open the galfit result FITS file
            if (True): #try:
                hdulist = pyfits.open(galfit_fullfn)
                hdr = hdulist[2].header
                print("opened file %s" % (galfit_fullfn))


                # now check all components - max # of components is 100
                fit_results = []
                param_labels = []

                for key in ['CHISQ', 'NDOF', 'NFREE', 'NFIX', 'CHI2NU']:
                    fit_results.extend(read_results(hdr, None, None, keyname=key))
                    param_labels.append("galfit_stats - %s" % (key))

                for component in range(1,100):
                    try:
                        comp_name = hdr['COMP_%d' % (component)]
                    except KeyError:
                        break

                    if (comp_name == "sersic"):
                        parameters = ['XC', 'YC', 'MAG', 'RE', 'N', 'AR', 'PA']
                    elif (comp_name == "sky"):
                        parameters = ['XC', 'YC', 'SKY', 'DSDX', 'DSDY']

                    for p in parameters:
                        fit_results.extend(read_results(hdr, component, p))
                        param_labels.append("%s_%d: %s" % (comp_name, component, p))

                print(fit_results)

                galfit_data[i_src] = fit_results

        # Now check all results for the longest value chain
        # ideally all results should have the same length
        n_parameters = [0 if fr is None else len(fr) for fr in galfit_data]
        print(n_parameters)
        print(set(n_parameters))

        # add additional header information about the galfit columns
        galfit_header = ["# %3d %s" % (3*i+1+catalog.shape[1], pl) for i,pl in enumerate(param_labels)]
        print(galfit_header)

        catalog_header.extend(galfit_header)
        print(catalog_header)

        param_count = numpy.array(list(set(n_parameters)))
        print(param_count)
        if (numpy.sum((param_count > 0)) > 1):
            print("Number of parameters returned from GALFIT do not match")

        else:
            # We have valid results - all frames either returned 0 parameters
            # (e.g. in the case the fit failed) or the X parameters
            galfit_results = numpy.empty((catalog.shape[0], numpy.max(param_count)))
            print(galfit_results.shape)
            for i, fr in enumerate(galfit_data):
                galfit_results[i,:] = numpy.NaN if fr is None else galfit_data[i]

            # Now merge the two arrays - starting with the numbers we got from
            # source extractor, followed by the ones from galfit
            combined = numpy.append(catalog, galfit_results, axis=1)

            # write output
            print(catalog_header)
            numpy.savetxt(combined_cat, combined,
                          header="\n".join(catalog_header), comments='')

            print("Combined catalog written to %s" % (combined_cat))

            # except:
            #     continue

        catalog_queue.task_done()



if __name__ == "__main__":


    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    # cmdline.add_argument("--conf", dest="sex_conf", default="sex.conf",
    #                      help="source extractor config filename")
    # cmdline.add_argument("--params", dest="sex_params", default="default.param",
    #                      help="number of models to insert into each image")
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    # cmdline.add_argument("--exe", dest="sex_exe", default="sex",
    #                      help="location of SExtractor executable")
    # cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
    #                      help="weight map")
    cmdline.add_argument("--catext", dest="catalog_extension", type=str, default="udgcat",
                         help="file extension for catalog")
    cmdline.add_argument("--subdir", dest="galfit_directory", type=str, default="galfit/",
                         help="output subdirectory to hold galfit feed-files and output")
    cmdline.add_argument("input_catalogs", nargs="+",
                         help="list of input source catalogs")
    #cmdline.print_help()
    args = cmdline.parse_args()


    # Feed parallel workers
    catalog_queue = multiprocessing.JoinableQueue()
    for cat_fn in args.input_catalogs:
        udg_cat = cat_fn[:-4] + "." + args.catalog_extension #".udgcat"
        final_cat = cat_fn[:-4] + ".galcomb.cat2"
        print("Adding %s --> %s" % (cat_fn, udg_cat))
        catalog_queue.put((cat_fn, udg_cat, final_cat))

    # insert termination commands
    for i in range(args.number_processes):
        catalog_queue.put((None))

    # start all processes
    processes = []
    for i in range(args.number_processes):
        print("Starting worker")
        p = multiprocessing.Process(
            target=parallel_combine,
            kwargs=dict(
                catalog_queue=catalog_queue,
                galfit_directory=args.galfit_directory,
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    catalog_queue.join()
