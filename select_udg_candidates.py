#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
import multiprocessing
import scipy
import scipy.special

import astropy.table

def select(catalog):

    try:
        mag_auto = catalog[:,5]
        a_image = catalog[:,16]
        mu0 = catalog[:,20]
        r50 = catalog[:,21]
    except IndexError:
        return None

    # udg = (mag_auto < 24) & (a_image > 5) & (mu0 > 24) & (r50 > 6)
    udg = (mag_auto < 24) & (mu0 > 24) & (r50 > 6)
    
    # criteria_count = numpy.sum((mag_auto < 24)) + numpy.sum((a_image > 6)) + \
    #                  numpy.sum((mu0 > 24)) + numpy.sum((r50 > 6))
    criteria_count = (mag_auto < 24).astype(numpy.int) + \
                     (a_image > 6).astype(numpy.int) + \
                     (mu0 > 24).astype(numpy.int) + \
                     (r50 > 6).astype(numpy.int)

    # print(criteria_count.shape)
    # udg = (criteria_count >= 3)

    return catalog[udg]


def select_from_galfit(catalog):

    try:
        mag_auto = catalog[:, 5]
        a_image = catalog[:, 16]
        mu0 = catalog[:, 20]

        r50 = catalog[:, 21]
        radius_eff = catalog[:, 52]
        magnitude = catalog[:, 49]
        sersic = catalog[:, 55]
        axis_ratio = catalog[:, 58]

    except IndexError:
        return None

    pixelscale = 0.168

    kappa = 2*sersic - 0.33
    g = scipy.special.gamma(2*sersic)
    m = magnitude
    f_tot = scipy.power(10., -0.4*m)
    sig_e = f_tot/(2*scipy.pi*(radius_eff**2)*((scipy.e)**kappa)*sersic*(kappa**(-2*sersic))*g*axis_ratio)

    sig_0 = sig_e * numpy.exp(kappa)

    #print("sig_e = %s" % sig_e)
    #print("sig_0 = %s" % sig_0)

    surfbrite_0 = -2.5*numpy.log10(sig_0) + 5*numpy.log10(pixelscale)

    #print(surfbrite_0)

    # udg = (mag_auto < 24) & (a_image > 5) & (mu0 > 24) & (r50 > 6)
    udg = (mag_auto < 24) & (mu0 > 24) & (r50 > 6)
    udg = (surfbrite_0 > 23) & (axis_ratio > 0.5) & (radius_eff > 9)

    # criteria_count = numpy.sum((mag_auto < 24)) + numpy.sum((a_image > 6)) + \
    #                  numpy.sum((mu0 > 24)) + numpy.sum((r50 > 6))
    criteria_count = (mag_auto < 24).astype(numpy.int) + \
                     (a_image > 6).astype(numpy.int) + \
                     (mu0 > 24).astype(numpy.int) + \
                     (r50 > 6).astype(numpy.int)

    # print(criteria_count.shape)
    # udg = (criteria_count >= 3)

    # udg = (numpy.isfinite(surfbrite_0)) & (sersic > 0.8)  & (sersic < 6) & (radius_eff > 8)

    return numpy.append(catalog, surfbrite_0.reshape((-1,1)), axis=1)[udg]
    # return catalog[udg]

#This is where I need to input my criteria

def select_maybeUDG(catalog):

    try:
        mean_surface_brightness_flux = (catalog['FLUX_AUTO'] / 2) / catalog['FLUX_RADIUS_50']
        mean_surface_brightness = -2.5*numpy.log10(mean_surface_brightness_flux)
        
        catalog['MEAN_SB'] = mean_surface_brightness
        catalog['MEAN_SB_FLUX'] = mean_surface_brightness_flux

        select = \
            (catalog['CLASS_STAR'] < 0.9) & \
            (catalog['CLASS_STAR'] > 0.00031974) & \
            (catalog['MU_MAX'] < 29.672) & \
            (catalog['MAG_AUTO'] < 35) & \
            (catalog['KRON_RADIUS'] > 0) & \
            (catalog['FWHM_IMAGE'] > 5) & \
            (catalog['FWHM_IMAGE'] < 256) & \
            (catalog['PETRO_RADIUS'] > 0) & \
            (catalog['FLAGS'] < 4) & \
            (catalog['MEAN_SB'] < 0.03) & \
            (catalog['MEAN_SB'] > -7) & \
            (catalog['A_IMAGE'] > 2.94) & \
            (catalog['ELONGATION'] < 2.64) & \
            (catalog['FLUX_AUTO'] > 7) & \
            (catalog['FLUX_RADIUS_50'] > 4.5)
                 
    except:
        pass
        return None

    #print(catalog['NUMBER'].shape, numpy.sum(select))
    final_catalog = catalog[select]
#    if(os.path.isfile(output_fn)):
#        os.remove(output_fn)
#    final_catalog.write(output_fn, format='votable')
    return final_catalog

def parallel_select(catalog_queue, selection=None, redo=False):

    if (selection is None):
        selection = 'sextractor'

    #print(selection)
    if (selection == 'galfit'):
        print("Using Amy's selection function")
        selector_fct = select_from_galfit
    elif (selection == "all"):
        selector_fct = select_all
    elif (selection == "maybeUDGs"):
        print("Using maybeUDGs for selection")
        selector_fct = select_maybeUDG


    # print("Hello from worker")
    while (True):
        cmd = catalog_queue.get()
        if (cmd is None):
            catalog_queue.task_done()
            # print("Shutting down")
            break

        cat_fn, output_fn, output_reg = cmd

        if (os.path.isfile(output_fn) and not redo):
            catalog_queue.task_done()
            continue

        try:
            # catalog = numpy.loadtxt(cat_fn)
            catalog = astropy.table.Table.read(cat_fn)
        except:
            print("Error opening %s" % (cat_fn))
            catalog_queue.task_done()
            continue

        # if (catalog.ndim < 2 or catalog.shape[0] <= 0):
        #     print("Error with catalog %s" % (cat_fn))
        #     catalog_queue.task_done()
        #     continue

        # # Now that we have the data, also read all the header information
        # with open(cat_fn, "r") as  cf:
        #     lines = cf.readlines()
        #     header = []
        #     for l in lines:
        #         if (l.startswith("#")):
        #             header.append(l.strip())
        #     # header = [l if l.startswith("#") for l in lines]
        # # print("header:\n", header)

        udg_candidates = selector_fct(catalog)

        if (udg_candidates is None):
            print("Error with catalog %s" % (cat_fn))
            catalog_queue.task_done()
            continue

        print(cat_fn, catalog.as_array().shape[0], "-->", udg_candidates.as_array().shape[0])

        # numpy.savetxt(output_fn, udg_candidates,
        #               header="\n".join(header), comments='')
        udg_candidates.write(output_fn, format='votable')

        try:
            if (output_reg is not None):
                with open(output_reg, "w") as reg:
                    print(output_reg)
                    hdr = """# Region file format: DS9 version 4.1
                             global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
                             fk5"""
                    reg.write("\n".join([a.strip() for a in hdr.splitlines()]))

                    labels = ["""# text(%f,%+f) font="helvetica 14 normal roman" text={%d}""" % (
                        src['ALPHA_J2000'], src['DELTA_J2000'], src['NUMBER']) for src in udg_candidates]
                    reg.write("\n".join(labels))
        except:
            pass

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

    cmdline.add_argument("--ext", dest="catalog_extension", type=str, default=".cat:.udgcat",
                         help="file extension for catalog")
    cmdline.add_argument("--outreg", dest="region_extension", type=str, default=None,
                         help="file extension for catalog")

    cmdline.add_argument("--select", dest="selection_mode", type=str, default="sextractor",
                         help="selection mode [sextractor/galfit]")
    cmdline.add_argument("--redo", dest="redo", default=False, action='store_true',
                         help="re-run even if output already exists")

    # cmdline.add_argument("--exe", dest="sex_exe", default="sex",
    #                      help="location of SExtractor executable")
    # cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
    #                      help="weight map")
    cmdline.add_argument("input_catalogs", nargs="+",
                         help="list of input source catalogs")
    #cmdline.print_help()
    args = cmdline.parse_args()


    # Feed parallel workers
    _parts = args.catalog_extension.split(":")
    _search = _parts[0]
    _replace = _parts[1]
    print(args)

    catalog_queue = multiprocessing.JoinableQueue()
    for cat_fn in args.input_catalogs:

        output_cat = cat_fn.replace(_search, _replace)
        # output_cat = cat_fn[:-4]+".udgcat"
        if (args.region_extension is not None):
            output_reg = cat_fn[:-4]+args.region_extension
        else:
            output_reg = None

        print("Adding %s --> %s" % (cat_fn, output_cat))
        catalog_queue.put((cat_fn, output_cat, output_reg))

    # insert termination commands
    for i in range(args.number_processes):
        catalog_queue.put((None))

    # start all processes
    processes = []
    for i in range(args.number_processes):
        print("Starting worker")
        p = multiprocessing.Process(
            target=parallel_select,
            kwargs=dict(
                catalog_queue=catalog_queue,
                selection=args.selection_mode,
                redo=args.redo,
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    catalog_queue.join()
            # load catalog
