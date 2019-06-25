#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
import multiprocessing
import pyfits
import logging
import logsetup
logsetup.setup_log()

import astropy.table


def read_results(hdr, component, parameter, keyname=None, x1=0, y1=0):

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

    if (parameter == "XC"):
        value += x1 #hdr['SRC_X1']
    if (parameter == "YC"):
        value += y1

    return [value, uncert, flag]


def prepare_for_galfit(catalog, components):

    cols2add = []
    keylist = []

    # TODO: Add code to handle multiple instances of the same model

    for i_model, model in enumerate(components):

        params = None
        if (model.lower() == "sersic"):
            params = ['XC', 'YC', 'MAG', 'RE', 'N', 'AR', 'PA']
        elif (model.lower() == "sky"):
            params = ["XC", 'YC', 'SKY', 'DSDX', 'DSDY']
        else:
            print("ERROR handlling components __%s__" % (model))
            continue

        modelname = model.upper()
        for p in params:
            val = "%s_%s" % (modelname, p)
            err = "%s_%s_ERR" % (modelname, p)
            cols2add.extend([val,err])

            keylist.append(("%d_%s" % (i_model+1, p), (val,err)))


    # also add some standard columns
    std_cols = ['CHISQ', 'NDOF', 'NFREE', 'NFIX', 'CHI2NU']
    cols2add.extend(std_cols)
    keylist.extend([(x,(x,)) for x in std_cols])

    # now add the columns to the catalog
    # catalog.add_columns( cols2add, )
    n_entries = catalog.as_array().shape[0]
    print(n_entries)
    print("Catalog has %d entries" % (n_entries))
    dummy_col = numpy.full((n_entries), fill_value=numpy.NaN, dtype=numpy.float32)

    # print(dummy_col)
    for column_name in cols2add:
        col = astropy.table.Column(name=column_name, data=dummy_col)
        catalog.add_column(col)

    return catalog, keylist, cols2add


def parallel_combine(catalog_queue, galfit_directory='galfit', components=None):

    # print("Hello from worker")
    logger = logging.getLogger("CombineSexGalfit")

    while (True):
        cmd = catalog_queue.get()
        if (cmd is None):
            catalog_queue.task_done()
            # print("Shutting down")
            break

        img_fn, udg_cat, combined_cat = cmd
        logger.debug("Reading %s" % (udg_cat))

        try:
            # catalog = numpy.loadtxt(udg_cat)
            catalog = astropy.table.Table.read(udg_cat)
        except:
            logger.warning("Error opening %s" % (udg_cat))
            catalog_queue.task_done()
            continue

        # print("xxx")
        # if (catalog.ndim < 2 or catalog.shape[0] <= 0):
        #     logger.warning("Error with catalog %s" % (udg_cat))
        #     catalog_queue.task_done()
        #     continue

        # with open(udg_cat, "r") as  cf:
        #     lines = cf.readlines()
        #     catalog_header = []
        #     for l in lines:
        #         if (l.startswith("#")):
        #             catalog_header.append(l.strip())
        #     # print(catalog_header)

        # Now add the components we need to hold the galfit results
        catalog, keylist, cols_added = prepare_for_galfit(catalog, components)

        #
        # Now we can start the actual work
        #
        _dir, _fn = os.path.split(udg_cat)
        basename, _ = os.path.splitext(_fn)
        if (galfit_directory.find(":") > 0):
            _parts = galfit_directory.split(":")
            _search = _parts[0]
            _replace = _parts[1]
            galfit_dir = img_fn.replace(_search, _replace)
        else:
            galfit_dir = os.path.join(_dir, galfit_directory)

        # print(galfit_dir)

        # galfit_data = [None] * catalog.shape[0]
        for i_src, src in enumerate(catalog):

            src_id = int(src['NUMBER'])

            galfit_fn = "%s.%05d.galfit.fits" % (basename, src_id)
            galfit_fullfn = os.path.join(galfit_dir, galfit_fn)
            print(galfit_fullfn)

            if (not os.path.isfile(galfit_fullfn)):
                logger.error("File does not exist - failed fit? (%s)" % (galfit_fullfn))
                continue

            # open the galfit result FITS file
            if (True): #try:
                try:
                    hdulist = pyfits.open(galfit_fullfn)
                    hdr = hdulist[2].header
                    x1 = hdulist[1].header["SRC_X1"] if "SRC_X1" in hdulist[1].header else 0.0
                    y1 = hdulist[1].header["SRC_Y1"] if "SRC_Y1" in hdulist[1].header else 0.0
                except:
                    print("There was a problem getting SRC_X1")
                    continue

                # print(x1,y1)
                # logger.debug("opened file %s" % (galfit_fullfn))

                # Now read all Galfit results from file header and insert
                # into the catalog
                header = hdulist[2].header
                for (header_key, catalog_key) in keylist:

                    if (len(catalog_key) == 1):
                        # there is only a fixed value
                        # print(catalog_key, header_key)
                        catalog[catalog_key[0]][i_src] = header[header_key]
                    else:
                        # there are two values (value and uncertainty)
                        value_key, error_key = catalog_key
                        error = numpy.NaN

                        fits_value = header[header_key]

                        if (fits_value.startswith("[")):
                            # value is of form 2_XC = '[54.0000]'
                            value = float(fits_value.split("[")[1].split("]")[0])

                        elif (fits_value.find("+/-") > 0):
                            # value looks like
                            # 1_XC    = '55.3318 +/- 0.4551'
                            items = [f.strip() for f in fits_value.split("+/-")]
                            # check if any of the items is problematic
                            problematic = numpy.array([f.startswith("*") for f in items]).any()
                            if (problematic):
                                # print(items)
                                pass
                            else:
                                value = float(items[0])
                                error = float(items[1])

                        else:
                            print("Unable to understand Galfit Result: %s = %s (%s)" % (
                                header_key, fits_value, galfit_fullfn
                            ))
                            continue

                        catalog[value_key][i_src] = value
                        catalog[error_key][i_src] = error
                        # print("Setting %s[%d] --> %f %f" % (header_key, i_src, value, error))

                # # now check all components - max # of components is 100
                # fit_results = []
                # param_labels = []
                #
                # for key in ['CHISQ', 'NDOF', 'NFREE', 'NFIX', 'CHI2NU']:
                #     fit_results.extend(read_results(hdr, None, None, keyname=key))
                #     param_labels.append("galfit_stats - %s" % (key))
                #
                # for component in range(1,100):
                #     try:
                #         comp_name = hdr['COMP_%d' % (component)]
                #     except KeyError:
                #         break
                #
                #     if (comp_name == "sersic"):
                #         parameters = ['XC', 'YC', 'MAG', 'RE', 'N', 'AR', 'PA']
                #     elif (comp_name == "sky"):
                #         parameters = ['XC', 'YC', 'SKY', 'DSDX', 'DSDY']
                #
                #     for p in parameters:
                #         fit_results.extend(read_results(hdr, component, p, x1=x1, y1=y1))
                #         param_labels.append("%s_%d: %s" % (comp_name, component, p))
                #
                # # print(fit_results)
                # if (len(fit_results) != 51):
                #     continue
                #
                # galfit_data[i_src] = fit_results

        # # Now check all results for the longest value chain
        # # ideally all results should have the same length
        # n_parameters = [0 if fr is None else len(fr) for fr in galfit_data]
        # # print(n_parameters)
        # # print(set(n_parameters))
        #
        # # add additional header information about the galfit columns
        # galfit_header = ["# %3d %s" % (3*i+1+catalog.shape[1], pl) for i,pl in enumerate(param_labels)]
        # # print(galfit_header)
        #
        # catalog_header.extend(galfit_header)
        # # print(catalog_header)
        #
        # param_count = numpy.array(list(set(n_parameters)))
        # # print(param_count)
        # if (numpy.sum((param_count > 0)) > 1):
        #     logger.error("Number of parameters returned from GALFIT do not match")
        #     print(param_count)
        #
        # else:
        #     # We have valid results - all frames either returned 0 parameters
        #     # (e.g. in the case the fit failed) or the X parameters
        #     galfit_results = numpy.empty((catalog.shape[0], numpy.max(param_count)))
        #     # print(galfit_results.shape)
        #     for i, fr in enumerate(galfit_data):
        #         galfit_results[i,:] = numpy.NaN if fr is None else galfit_data[i]
        #
        #     # Now merge the two arrays - starting with the numbers we got from
        #     # source extractor, followed by the ones from galfit
        #     combined = numpy.append(catalog, galfit_results, axis=1)
        #
        #     # write output
        #     # print(catalog_header)
        #     numpy.savetxt(combined_cat, combined,
        #                   header="\n".join(catalog_header), comments='')
        #
        #     logger.info("Combined catalog (%d sources) written to %s" % (catalog.shape[0], combined_cat))
        #
        #     # except:
        #     #     continue


        # lastly, write the updated catalog
        catalog.info()
        catalog.write(combined_cat, format='votable', overwrite=True)

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
    cmdline.add_argument("--components", dest="components", type=str, default="sersic,sky",
                         help="list of galfit components")
    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()


    # Feed parallel workers
    _parts = args.catalog_extension.split(":")
    _input_cat = _parts[0]
    _output_cat = _parts[1]
    print(args)

    galfit_components = args.components.split(",")

    # Feed parallel workers
    catalog_queue = multiprocessing.JoinableQueue()
    for img_fn in args.input_images:
        udg_cat = img_fn[:-5] + _input_cat #"." + args.catalog_extension #".udgcat"
        final_cat = img_fn[:-5] + _output_cat #".galcomb.cat2"
        print("Adding %s: %s --> %s" % (img_fn, udg_cat, final_cat))
        catalog_queue.put((img_fn, udg_cat, final_cat))

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
                components=galfit_components
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    catalog_queue.join()

