#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
import scipy.spatial


def read_header(fn):

    columns = []
    with open(fn, "r") as f:
        lines = f.readlines()
        columns = [str(l.strip()) for l in lines]

    return columns


def get_headers(columns, n_columns, startat=0):

    hdr = []
    for n in range(n_columns):
        label = columns[n] if  n < len(columns) else "???"
        hdr.append("# column %3d: %s" % (n+1+startat, label))

    return hdr


if __name__ == "__main__":


    # setup command line parameters
    cmdline = argparse.ArgumentParser(
        epilog="""\
This tool is meant to match catalogs. It assumes that the
Ra/Dec coordinates used for matching are the first two 
columns in all catalogs."""
    )
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
    cmdline.add_argument("--rmatch", dest="matching_radius", default=1.0, type=float,
                         help="matching radius for center positions [arcsec]")

    cmdline.add_argument("--ra", dest="ra", default=194.953054, type=float,
                         help="center position of field in RA [deg]")
    cmdline.add_argument("--dec", dest="dec", default=27.980694, type=float,
                         help="center position of field in DEC [deg]")

    cmdline.add_argument("cat1", nargs=1,
                         help="first input source catalogs")
    cmdline.add_argument("cat2", nargs=1,
                         help="second input source catalogs")

    cmdline.add_argument("--header1", dest="header1", default=None, type=str,
                         help="column information for catalog #1")
    cmdline.add_argument("--header2", dest="header2", default=None, type=str,
                         help="column information for catalog #2")

    cmdline.add_argument("--out", dest="output_fn", default=None, type=str,
                         help="output filename")

    #cmdline.print_help()
    args = cmdline.parse_args()
    print(args)


    # load 1st catalog
    print(args.cat1)
    cat1 = numpy.loadtxt(args.cat1[0])
    print("First catalog:  %d entries" % (cat1.shape[0]))
    columns1 = read_header(args.header1)

    # load 2nd catalog
    cat2 = numpy.loadtxt(args.cat2[0])
    print("Second catalog: %d entries" % (cat2.shape[0]))
    columns2 = read_header(args.header2)


    # convert coordinates into projected differences
    cos_dec = numpy.cos(numpy.radians(args.dec))
    radec1 = (cat1[:, 0:2] - [args.ra, args.dec]) * [cos_dec, 1.0]
    radec2 = (cat2[:, 0:2] - [args.ra, args.dec]) * [cos_dec, 1.0]
    matching_radius = args.matching_radius / 3600.
    # dec1 = cat1[:,1] - args.dec
    # dec2 = cat2[:,1] - args.dec
    # ra1 = (cat1[:,0] - args.ra) * cos_dec
    # ra2 = (cat2[:,0] - args.ra) * cos_dec

    print("Creating ref KDtree")
    ref_tree = scipy.spatial.cKDTree(radec1) #cat1[:, 1:3]) # use x/y here
    src_tree = scipy.spatial.cKDTree(radec2) #cat2[:, 1:3])

    # d,i = ref_tree.query(udgcat[:, 3:5],
    #                      k=1, p=2, distance_upper_bound=args.matching_radius)

    d,i = src_tree.query(radec1,
                         k=1, p=2, distance_upper_bound=matching_radius)

    valid_match = numpy.isfinite(d)
    # print(valid_match)
    # print(i)

    match1 = cat1[valid_match]
    match2 = cat2[i[valid_match]]
    print("match dimensions:", match1.shape, match2.shape)
    print("input dimensions:", cat1.shape, cat2.shape)
    print("column data:", len(columns1), len(columns2))

    # print(numpy.sum(valid_match), match_log.shape, match_srccat.shape)

    # merge both catalogs
    merged = numpy.append(match1, match2, axis=1)
    print(merged.shape)
    # print(merged.shape)



    # merge the column headers
    headers = get_headers(columns1, match1.shape[1], startat=0)
    # print(headers)

    headers2 = get_headers(columns2, match2.shape[1], startat=match1.shape[1])
    # print(headers2)

    headers.extend(headers2)

    # headers.extend(get_headers(columns2, match2.shape[1], startat=match1.shape[1]+1))

    if (args.output_fn is None):
        out_f = sys.stdout
    else:
        print("Writing output to %s" % (args.output_fn))
        out_f = open(args.output_fn, "wb")

    # out_f.write(asbytes("\n".join(headers)+"\n"))
    # #print(merged[0])
    # for i in range(merged.shape[0]):
    #     print(i)
    #     numpy.savetxt(out_f, merged[i:i+1])
    numpy.savetxt(out_f, merged, header="\n".join(headers), comments='')


    out_f.close()

