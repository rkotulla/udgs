#!/usr/bin/env python

import os
import sys
import numpy
import argparse
import scipy.spatial

if __name__ == "__main__":


    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    # cmdline.add_argument("--conf", dest="sex_conf", default="sex.conf",
    #                      help="source extractor config filename")
    # cmdline.add_argument("--params", dest="sex_params", default="default.param",
    #                      help="number of models to insert into each image")
    # cmdline.add_argument("--nprocs", dest="number_processes",
    #                      default=multiprocessing.cpu_count(), type=int,
    #                      help="number of Sextractors to run in parallel")
    # cmdline.add_argument("--exe", dest="sex_exe", default="sex",
    #                      help="location of SExtractor executable")
    # cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
    #                      help="weight map")
    cmdline.add_argument("--rmatch", dest="matching_radius", default=10, type=float,
                         help="matching radius for center positions [pixels]")
    cmdline.add_argument("input_catalogs", nargs="+",
                         help="list of input source catalogs")
    #cmdline.print_help()
    args = cmdline.parse_args()

    total_sim_count = {}
    all_merged = {}

    for cat_fn in args.input_catalogs:

        #
        # Load all catalogs
        #
        udgcat_fn = cat_fn[:-4]+".udgcat"
        log_fn = cat_fn[:-4]+".log"

        if (not os.path.isfile(udgcat_fn)):
            print("Unable to find %s" % (udgcat_fn))
            continue

        print("Checking completeness for %s" % (udgcat_fn))

        udgcat = numpy.loadtxt(udgcat_fn)
        log = numpy.loadtxt(log_fn)

        dir_name, _ = os.path.split(udgcat_fn)
        if (dir_name not in total_sim_count):
            total_sim_count[dir_name] = 0
        total_sim_count[dir_name] += log.shape[0]

        #
        # Now match catalogs and identify the galaxies we have recovered
        #
        # print(udgcat.shape)

        # print("Creating ref KDtree")
        ref_tree = scipy.spatial.cKDTree(log[:, 1:3]) # use x/y here
        src_tree = scipy.spatial.cKDTree(udgcat[:, 3:5])

        # d,i = ref_tree.query(udgcat[:, 3:5],
        #                      k=1, p=2, distance_upper_bound=args.matching_radius)

        d,i = src_tree.query(log[:, 1:3],
                             k=1, p=2, distance_upper_bound=args.matching_radius)

        valid_match = numpy.isfinite(d)
        # print(valid_match)
        # print(i)

        match_log = log[valid_match]
        match_srccat = udgcat[i[valid_match]]

        # print(numpy.sum(valid_match), match_log.shape, match_srccat.shape)

        # merge both catalogs
        merged = numpy.append(match_log, match_srccat, axis=1)
        # print(merged.shape)

        if (dir_name not in all_merged):
            all_merged[dir_name] = merged
        else:
            all_merged[dir_name] = numpy.append(all_merged[dir_name], merged, axis=0)

    # print(all_merged.shape)
    #
    # folder_name, _ = os.path.split(args.input_catalogs[0])
    # print(folder_name)
    # matches_fn = "%s.matches" % (folder_name)
    # numpy.savetxt(matches_fn, all_merged)
    #
    # summary_fn = "%s.summary" % (folder_name)
    # # save quick summary

    print(total_sim_count)

    #
    # Now get some summary information for all simulation runs
    #
    summary = []
    for dir_name in total_sim_count:
        param_file = os.path.join(dir_name, ".params")
        params = numpy.loadtxt(param_file)

        summary.append([params[0], params[1], params[2], params[3],
                        params[4], params[5],
                        total_sim_count[dir_name],
                        all_merged[dir_name].shape[0]])
    summary = numpy.array(summary)
    numpy.savetxt("completeness.summary", summary)