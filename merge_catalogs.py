#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing
import threading

import astropy.table




def catalog_reader(file_queue, catalog_queue, thread_id=0):

    while (True):

        fn = file_queue.get()
        if (fn is None):
            file_queue.task_done()
            break

        print("[CPU %2d] Reading catalog from %s ..." % (thread_id, fn))
        cat = astropy.table.Table.read(fn)

        catalog_queue.put(cat)
        file_queue.task_done()


if __name__ == "__main__":


    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("-o", "--output", dest="output", type=str,
                         default="combined_catalog.vot",
                         help="name of output catalog")
    cmdline.add_argument("input", nargs="+",
                         help="list of input images")
    args = cmdline.parse_args()



    combined_catalog = None
    file_queue = multiprocessing.JoinableQueue()
    catalog_queue = multiprocessing.Queue()

    # queue up all files
    for fn in args.input:
        file_queue.put(fn)

    # add termination tags
    for n in range(args.number_processes):
        file_queue.put((None))

    processes = []
    for n in range(args.number_processes):
        p = multiprocessing.Process(
            target=catalog_reader,
            kwargs=dict(file_queue=file_queue,
                        catalog_queue=catalog_queue,
                        thread_id=n+1),
        )

        # This code is more multi-threading instead of multi-processing
        # p = threading.Thread(
        #     target=catalog_reader,
        #     kwargs=dict(file_queue=file_queue,
        #                  catalog_queue=catalog_queue,
        #                 thread_id=n+1),
        # )
        # p.daemon = True

        p.start()

    # weit for processes to finish
    file_queue.join()

    catalogs = []
    for i in range(len(args.input)):
        cat = catalog_queue.get()
        catalogs.append(cat)

    # now merge all catalogs
    print("writing combined catalog")
    combined = astropy.table.vstack(catalogs)
    combined.write(args.output, format='votable', overwrite=True)

    print("all done!")

