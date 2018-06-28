#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
import multiprocessing


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



def parallel_select(catalog_queue):

    # print("Hello from worker")
    while (True):
        cmd = catalog_queue.get()
        if (cmd is None):
            catalog_queue.task_done()
            # print("Shutting down")
            break

        cat_fn, output_fn = cmd
        if (os.path.isfile(output_fn)):
            catalog_queue.task_done()
            continue

        try:
            catalog = numpy.loadtxt(cat_fn)
        except:
            print("Error opening %s" % (cat_fn))
            catalog_queue.task_done()
            continue

        if (catalog.ndim < 2 or catalog.shape[0] <= 0):
            print("Error with catalog %s" % (cat_fn))
            catalog_queue.task_done()
            continue

        udg_candidates = select(catalog)
        if (udg_candidates is None):
            print("Error with catalog %s" % (cat_fn))
            catalog_queue.task_done()
            continue

        print(cat_fn, catalog.shape, "-->", udg_candidates.shape)

        numpy.savetxt(output_fn, udg_candidates)
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
    cmdline.add_argument("input_catalogs", nargs="+",
                         help="list of input source catalogs")
    #cmdline.print_help()
    args = cmdline.parse_args()


    # Feed parallel workers
    catalog_queue = multiprocessing.JoinableQueue()
    for cat_fn in args.input_catalogs:
        output_cat = cat_fn[:-4]+".udgcat"
        print("Adding %s --> %s" % (cat_fn, output_cat))
        catalog_queue.put((cat_fn, output_cat))

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
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    catalog_queue.join()
            # load catalog
