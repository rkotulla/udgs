#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import time

import astropy.io.fits as pyfits
import numpy
import ldac2vot

def run_sex(file_queue, sex_exe, sex_conf, sex_param, fix_vot_array=None):

    while (True):

        opts = file_queue.get()
        if (opts is None):
            file_queue.task_done()
            print("Parallel worker received shutdown command")
            break


        img_fn, weight_fn = opts
        ldac_file = img_fn[:-5]+".fitsldac"
        seg_file = img_fn[:-5]+".segments"
        cat_file = img_fn[:-5]+".vot"

        fits_file = img_fn

        if (os.path.isfile(cat_file) and os.path.isfile(ldac_file)):
            file_queue.task_done()
            continue
            
        print("running sex on %s" % (img_fn))
        
        hdu = pyfits.open(img_fn)
        try:
            magzero = 2.5*numpy.log10(hdu[0].header['FLUXMAG0'])
        except:
            magzero = 0

        if (weight_fn is None):
            weight_opts = "-WEIGHT_TYPE NONE"
        else:
            weight_opts = """-WEIGHT_IMAGE "%s" """ % (weight_fn)

        sexcmd = """%s 
        -c %s 
        -PARAMETERS_NAME %s 
        %s 
        -CATALOG_NAME %s 
        -CATALOG_TYPE FITS_LDAC
        -CHECKIMAGE_TYPE SEGMENTATION
        -CHECKIMAGE_NAME %s
        -WEIGHT_THRESH 1e8
        -MAG_ZEROPOINT %.4f 
        %s """ % (
            sex_exe, sex_conf, sex_param,
            weight_opts,
            ldac_file,
            seg_file,
            magzero,
            fits_file)
        # print(" ".join(sexcmd.split()))

        start_time = time.time()
        try:
            # os.system(sexcmd)
            ret = subprocess.Popen(sexcmd.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
            # sextractor_pid = ret.pid
            # print("Started process ID %d" % (sextractor_pid))
            #
            (sex_stdout, sex_stderr) = ret.communicate()
            if (ret.returncode != 0):
                print("return code was not 0")
                print(sex_stdout)
                print(sex_stderr)

            #     break

            # sex = subprocess.run(sexcmd.split(), shell=True, check=True,
            #                      stdout=subprocess.PIPE,
            #                      stderr=subprocess.PIPE,)
        except OSError as e:
            print("Some exception has occured:\n%s" % (str(e)))
        end_time = time.time()
        print("SourceExtractor returned after %.3f seconds" % (end_time - start_time))

        # Now convert the FITS-LDAC catalog to VOTable format
        ldac2vot.fitsldac2vot(ldac_file, vot_fn=cat_file,
                              array_suffix=fix_vot_array)

        file_queue.task_done()



if __name__ == "__main__":

    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--conf", dest="sex_conf", default="sex.conf",
                         help="source extractor config filename")
    cmdline.add_argument("--params", dest="sex_params", default="default.param",
                         help="number of models to insert into each image")
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("--exe", dest="sex_exe", default="sex",
                         help="location of SExtractor executable")
    cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
                         help="weight map")
    cmdline.add_argument("--votfix", dest='fix_vot_arrays', type=str, default=None,
                         help="rename arrays when converting FITS-LDAC to VOTable")
    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()

    construct_weight_fn = False
    if (args.weight_image is not None):
        construct_weight_fn = (args.weight_image.find(":") > 0)

    fix_vot_array = ldac2vot.read_definitions(args.fix_vot_arrays)
    print(fix_vot_array)

    file_queue = multiprocessing.JoinableQueue()

    # feed with files to sextract
    for fn in args.input_images:
        if (os.path.isfile(fn)):

            if (construct_weight_fn):
                _parts = args.weight_image.split(":")
                search = _parts[0]
                replace = _parts[1]
                weight_fn = fn.replace(search, replace)
            else:
                weight_fn = args.weight_image

            file_queue.put((fn, weight_fn))

    # insert termination commands
    for i in range(args.number_processes):
        file_queue.put((None))

    # start all processes
    processes = []
    for i in range(args.number_processes):
        p = multiprocessing.Process(
            target=run_sex,
            kwargs=dict(
                file_queue=file_queue,
                sex_exe=args.sex_exe,
                sex_conf=args.sex_conf, sex_param=args.sex_params,
                fix_vot_array=fix_vot_array
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    # now wait for all work to be done
    file_queue.join()


