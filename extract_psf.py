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
import shutil

def run_sex_psfex(file_queue, sex_exe, sex_conf, sex_param,
                  psfex_exe, psfex_conf, supersample=2):

    while (True):

        opts = file_queue.get()
        if (opts is None):
            file_queue.task_done()
            print("Parallel worker received shutdown command")
            break


        img_fn, weight_fn = opts
        img_basename = img_fn[:-5]
        ldac_file = img_basename+".ldac4psfex"
        print(ldac_file)
        fits_file = img_fn

        # if (os.path.isfile(ldac_file)):
        #     file_queue.task_done()
        #     continue
            
        print("running sex on %s" % (img_fn))
        
        hdu = pyfits.open(img_fn)
        try:
            magzero = 2.5*numpy.log10(hdu[0].header['FLUXMAG0'])
        except:
            magzero = 0

        if (weight_fn is None):
            weight_opts = "-WEIGHT_TYPE NONE"
        else:
            weight_opts = """-WEIGHT_TYPE MAP_VAR -WEIGHT_IMAGE "%s" """ % (weight_fn)


        if (os.path.isfile(ldac_file)):
            print("Re-using existing LDAC catalog")
        else:
            sexcmd = """%s 
            -c %s 
            -PARAMETERS_NAME %s 
            %s 
            -CATALOG_NAME %s 
            -CATALOG_TYPE FITS_LDAC
            -WEIGHT_THRESH 1e8
            -MAG_ZEROPOINT %.4f 
            %s """ % (
                sex_exe, sex_conf, sex_param,
                weight_opts,
                ldac_file,
                magzero,
                fits_file)
            print(" ".join(sexcmd.split()))

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

        #
        # Now we have the FITS-LDAC catalog with all the little cutouts around stars that
        # we need to extract a valid PSF model using PSFEX
        #

        output_basename = img_basename[:-6]+"_psf" #"xxx_psf"
        check_opts = """
            -CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS
            -CHECKIMAGE_NAME %(bn)s.chi.fits,%(bn)s.proto.fits,%(bn)s.samp.fits,%(bn)s.resi.fits,%(bn)s.snap.fits
            -CHECKPLOT_TYPE FWHM,ELLIPTICITY,COUNTS,COUNT_FRACTION,CHI2,RESIDUALS
            -CHECKPLOT_NAME %(bn)s.fwhm,%(bn)s.ellipticity,%(bn)s.counts,%(bn)s.countfrac,%(bn)s.chi2,%(bn)s.resi
            """ % dict(bn=output_basename)

        psf_output_size = int(64 * supersample)
        if ((psf_output_size % 2) == 0):
            # size is even, we want an odd number
            psf_output_size += 1

        psfex_command = """
        %(psfex_exe)s 
        -c %(psfex_conf)s 
        %(checks)s
        -OUTCAT_TYPE FITS_LDAC
        -OUTCAT_NAME %(bn)s.out.cat 
        -PSF_SUFFIX .psf
        -WRITE_XML Y
        -XML_NAME %(bn)s.psfexlog.xml
        -PSF_SAMPLING %(pixelscale).3f
        -PSF_PIXELSIZE %(pixelscale).3f
        -PSF_SIZE %(psfsize)d,%(psfsize)d 
        %(ldac_catalog)s
        """ % dict(
            psfex_exe=psfex_exe,
            psfex_conf=psfex_conf,
            checks=check_opts,
            ldac_catalog=ldac_file,
            bn=output_basename,
            pixelscale=1./ supersample,
            psfsize=psf_output_size,
        )
        print(psfex_command)
        print(" ".join(psfex_command.split()))

        os.system(" ".join(psfex_command.split()))

        #
        # Now we need to rename a bunch of files to undo the PSFex file naming convention
        #
        for part in ["chi", "proto", "resi", "samp", "snap"]:
            bad_name = "%s.%s_%s.fits" % (output_basename, part, img_basename)
            good_name = "%s.%s.fits" % (output_basename, part)
            print("Copying file %s --> %s" % (bad_name, good_name))
            try:
                shutil.copyfile(bad_name, good_name)
                os.remove(bad_name)
            except:
                pass

        #
        # and finally, extract just the image of the PSF to be used for galfit
        #
        psf_input = "%s.proto.fits" % (output_basename)
        psf_output = "%s.fits" % (output_basename)
        psf_hdu = pyfits.open(psf_input)
        psf_size = psf_hdu[0].header['NAXIS2']
        psf_img = psf_hdu[0].data[:psf_size, :psf_size]
        psf_sum = numpy.sum(psf_img)
        psf_hdu[0].data = psf_img / psf_sum
        psf_hdu[0].header['SUPERSMP'] = supersample
        psf_hdu.writeto(psf_output, overwrite=True)

        # done with this image
        file_queue.task_done()



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
    cmdline.add_argument("--psfexconf", dest="psfex_conf", default=os.path.join(config_dir, "psfex.conf"),
                         help="number of models to insert into each image")
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("--sex", dest="sex_exe", default=distutils.spawn.find_executable("sex"),
                         help="location of SExtractor executable")
    cmdline.add_argument("--psfex", dest="psfex_exe", default=distutils.spawn.find_executable("psfex"),
                         help="location of SExtractor executable")
    cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
                         help="weight map")
    cmdline.add_argument("--psf", dest='psf_image', type=str, default="_image.fits:_psf.fits",
                         help="psf basename")
    cmdline.add_argument("--supersample", dest='supersample', type=float, default=2.,
                         help="supersample PSF model")
    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    # cmdline.print_help()
    args = cmdline.parse_args()

    construct_weight_fn = False
    if (args.weight_image is not None):
        construct_weight_fn = (args.weight_image.find(":") > 0)

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
            target=run_sex_psfex,
            kwargs=dict(
                file_queue=file_queue,
                sex_exe=args.sex_exe,
                sex_conf=args.sex_conf, sex_param=args.sex_params,
                psfex_conf=args.psfex_conf, psfex_exe=args.psfex_exe,
                supersample=args.supersample,
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    # now wait for all work to be done
    file_queue.join()


