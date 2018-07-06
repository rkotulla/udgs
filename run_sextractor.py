#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import time

def run_sex(file_queue, sex_exe, sex_conf, sex_param):

    while (True):

        opts = file_queue.get()
        if (opts is None):
            file_queue.task_done()
            print("Parallel worker received shutdown command")
            break


        img_fn, weight_fn = opts
        cat_file = img_fn[:-5]+".cat"
        seg_file = img_fn[:-5]+".segments"

        fits_file = img_fn

        if (os.path.isfile(cat_file)):
            file_queue.task_done()
            continue
            
        print("running sex on %s" % (img_fn))

        if (weight_fn is None):
            weight_opts = "-WEIGHT_TYPE NONE"
        else:
            weight_opts = """-WEIGHT_IMAGE "%s" """ % (weight_fn)

        sexcmd = """%s 
        -c %s 
        -PARAMETERS_NAME %s 
        %s 
        -CATALOG_NAME %s 
        -CATALOG_TYPE ASCII_HEAD
        -CHECKIMAGE_TYPE SEGMENTATION
        -CHECKIMAGE_NAME %s 
        %s """ % (
            sex_exe, sex_conf, sex_param,
            weight_opts,
            cat_file,
            seg_file,
            fits_file)
        # print(sexcmd)

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
    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()



    file_queue = multiprocessing.JoinableQueue()

    # feed with files to sextract
    for fn in args.input_images:
        if (os.path.isfile(fn)):
            file_queue.put((fn, args.weight_image))

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
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    # now wait for all work to be done
    file_queue.join()


