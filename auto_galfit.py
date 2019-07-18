#!/usr/bin/env python3

import os
import sys
import numpy
import astropy.io.fits as pyfits

import argparse
import multiprocessing
import time
import subprocess

import logging
import shutil
logging.basicConfig(filename='debug.log',level=logging.DEBUG)

import plot_galfit_results

import astropy.table
import shutil

def parallel_config_writer(queue, galfit_queue,
                           n_galfeeds, n_galfit_queuesize):

    counter = 0
    while (True):

        cmd = queue.get()
        if (cmd is None):
            queue.task_done()
            # galfit_queue.put((None))
            break

        image_fn, feedme_fn, weight_fn, src_info, \
            img_header, galfit_output, galfit_logfile, \
            segmentation_input_fn, galfit_dir, basename, \
            max_size, psf_file, psf_supersample, \
            magzero = cmd

        if (os.path.isfile(feedme_fn)):
            print("Skipping existing feed-file %s" % (feedme_fn))
            queue.task_done()
            if (n_galfit_queuesize is not None):
                with n_galfit_queuesize.get_lock():
                    n_galfit_queuesize.value += 1
            if (n_galfeeds is not None):
                with n_galfeeds.get_lock():
                    n_galfeeds.value += 1
            galfit_queue.put((feedme_fn, galfit_output, galfit_logfile))
            continue

        #
        # Create the galfit feedme file, the cutouts for the
        # image, weight-image, and segmentation mask
        #
        fwhm = src_info['FWHM_IMAGE']
        src_id = int(src_info['NUMBER'])
        size = 3 * fwhm
        if (max_size > 0 and size > max_size): size = max_size

        x, y = src_info['X_IMAGE']-1, src_info['Y_IMAGE']-1
        x1 = int(numpy.max([0, x - size]))
        x2 = int(numpy.min([x + size, img_header['NAXIS1']]))
        y1 = int(numpy.max([0, y - size]))
        y2 = int(numpy.min([y + size, img_header['NAXIS2']]))

        # open the input images and create cutouts
        img_out_fn = "%s/%s.%05d.image.fits" % (galfit_dir, basename, src_id)
        segm_out_fn = "%s/%s.%05d.segm.fits" % (galfit_dir, basename, src_id)
        weight_out_fn = "%s/%s.%05d.sigma.fits" % (galfit_dir, basename, src_id)
        constraints_opt = "%s.%05d.constraints" % (basename, src_id)
        constraints_fn = "%s/%s" % (galfit_dir, constraints_opt)
        psf_out_fn = "%s/%s.%05d.psf.fits" % (galfit_dir, basename, src_id)

        img_hdu = pyfits.open(image_fn)
        img = img_hdu[0].data[y1:y2, x1:x2]
        phdu = pyfits.PrimaryHDU(data=img)
        phdu.header['SRC_X1'] = x1
        phdu.header['SRC_Y1'] = y1
        phdu.writeto(img_out_fn, clobber=True)
        img_hdu.close()

        _, _img = os.path.split(img_out_fn)
        _weight, _bpm = 'none', 'none'
        _, _out = os.path.split(galfit_output)

        if (weight_fn is not None):
            wht_hdu = pyfits.open(weight_fn)
            wht = wht_hdu[0].data[y1:y2, x1:x2]
            pyfits.PrimaryHDU(data=wht, header=phdu.header).writeto(weight_out_fn, clobber=True)
            wht_hdu.close()
            _, _weight = os.path.split(weight_out_fn)

        if (segmentation_fn is not None):
            try:
                segm_hdu = pyfits.open(segmentation_fn)
                segm = segm_hdu[0].data[y1:y2, x1:x2].astype(numpy.int)
                segm[segm == src_id] = 0
                pyfits.PrimaryHDU(data=segm, header=phdu.header).writeto(segm_out_fn, clobber=True)
                segm_hdu.close()
                _, _bpm = os.path.split(segm_out_fn)
            except IOError:
                pass
        else:
            print("Unable to generate source mask from segmentation file (%s)" % (segmentation_fn))

        if (psf_file is not None and os.path.isfile(psf_file)):
            # copy the PSF model
            shutil.copyfile(psf_file, psf_out_fn)
            _, galfit_psf_option = os.path.split(psf_out_fn)

            if (psf_supersample <= 0):
                # this is meant to be auto-scaled
                try:
                    psf_hdu = open(psf_out_fn)
                    psf_supersample = psf_hdu[0].header['SUPERSMP']
                    psf_hdu.close()
                except:
                    psf_supersample = 1.
        else:
            galfit_psf_option = 'none'


        #
        # Generate the constraints file
        #
        dx = numpy.min([2., 3 * numpy.sqrt(src_info['ERRX2WIN_IMAGE']) + 1.])
        dy = numpy.min([2., 3 * numpy.sqrt(src_info['ERRY2WIN_IMAGE']) + 1.])
        with open(constraints_fn, "w") as cf:
            constraints = """
                1   x   %(dx).2f %(dx).2f
                1   y   %(dy).2f %(dy).2f    
                        
            """ % {
                'dx': dx,
                'dy': dy,
            }
            cf.write("\n".join([c.strip() for c in constraints.splitlines(keepends=False)]))


        galfit_info = {
            'imgfile': _img, #img_out_fn, #image_fn,
            'srcid': src_id,
            'x1': 0, #x1,
            'x2': x2-x1, #x2
            'y1': 0, #y1,
            'y2': y2-y1, #y2,
            'pixelscale': 0.18,
            'weight_image': _weight, #weight_out_fn if weight_fn is not None else 'none',
            'galfit_output': _out, #galfit_output,
            'bpm': _bpm, #segm_out_fn if segmentation_fn is not None else 'none',
            'psf': galfit_psf_option,
            'psf_supersample': int(psf_supersample),
            'magzero': magzero,
            'constraints': constraints_opt,
        }

        head_block = """
            A) %(imgfile)s         # Input data image (FITS file)
            B) %(galfit_output)s   # Output data image block
            C) %(weight_image)s                # Sigma image name (made from data if blank or "none") 
            D) %(psf)s   #        # Input PSF image and (optional) diffusion kernel
            E) %(psf_supersample)d                   # PSF fine sampling factor relative to data 
            F) %(bpm)s                # Bad pixel mask (FITS image or ASCII coord list)
            G) %(constraints)s                # File with parameter constraints (ASCII file) 
            H) %(x1)d %(x2)d %(y1)d %(y2)d   # Image region to fit (xmin xmax ymin ymax)
            I) 100    100          # Size of the convolution box (x y)
            J) %(magzero).3f              # Magnitude photometric zeropoint 
            K) %(pixelscale).3f %(pixelscale).3f            # Plate scale (dx dy)    [arcsec per pixel]
            O) regular             # Display type (regular, curses, both)
            P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        """ % (galfit_info)
            # print(head_block)

        posangle = 90 - src_info['THETA_IMAGE']
        src_info = {
            'x': x-x1,
            'y': y-y1,
            'magnitude': src_info['MAG_AUTO'],  # +magzero,
            'halflight_radius': src_info['FLUX_RADIUS_50'],
            'sersic_n': 1.5, #src[7],
            'axis_ratio': 1./src_info['ELONGATION'],  #sextractur uses a/b, galfit needs b/a
            'position_angle': posangle,

        }
        object_block = """
            # Object number: 1
             0) sersic                 #  object type
             1) %(x)d  %(y)d  1 1  #  position x, y
             3) %(magnitude).3f     1          #  Integrated magnitude	
             4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
             5) %(sersic_n).3f      1          #  Sersic index n (de Vaucouleurs n=4) 
             6) 0.0000      0          #     ----- 
             7) 0.0000      0          #     ----- 
             8) 0.0000      0          #     ----- 
             9) %(axis_ratio).3f      1          #  axis ratio (b/a)  
            10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
             Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
            
            # # Object number: 1
            #  0) devauc                 #  object type
            #  1) %(x)d  %(y)d  1 1  #  position x, y
            #  3) %(magnitude).3f     1          #  Integrated magnitude	
            #  4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
            #  9) %(axis_ratio).3f      1          #  axis ratio (b/a)  
            # 10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
            #  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
            # 
            # # Object number: 1
            #  0) expdisk                 #  object type
            #  1) %(x)d  %(y)d  1 1  #  position x, y
            #  3) %(magnitude).3f     1          #  Integrated magnitude	
            #  4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
            #  9) %(axis_ratio).3f      1          #  axis ratio (b/a)  
            # 10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
            #  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
            
            # Object number: 2
             0) sky                    #  object type
             1) 0.0000      1          #  sky background at center of fitting region [ADUs]
             2) 0.0000      0          #  dsky/dx (sky gradient in x)
             3) 0.0000      0          #  dsky/dy (sky gradient in y)
             Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
                
        """ % src_info
        # print(object_block)

        # feedme_fn = "feedme.%d" % (int(src[4]))
        # feedme_fn = "%s.src%05d.galfeed" % (config_basename, src_id)
        with open(feedme_fn, "w") as feedfile:
            feedfile.write("\n".join([l.strip() for l in head_block.splitlines()]))
            feedfile.write("\n".join([l.strip() for l in object_block.splitlines()]))


        # Now also prepare the segmentation mask, if available

        # queue up a new execution of galfit
        if (n_galfit_queuesize is not None):
            with n_galfit_queuesize.get_lock():
                n_galfit_queuesize.value += 1
        if (n_galfeeds is not None):
            with n_galfeeds.get_lock():
                n_galfeeds.value += 1
        galfit_queue.put((feedme_fn, galfit_output, galfit_logfile))
        counter += 1

        queue.task_done()
        continue

    # print("Prepared %d galfeeds, queue-size = %d %d" % (
    #     n_galfeeds.value, n_galfit_queuesize.value, counter))


dryrun = False


def parallel_run_galfit(galfit_queue, galfit_exe='galfit',
                        make_plots=True, redo=False,
                        n_galfit_complete=None, n_total_galfit_time=None,
                        n_galfit_queuesize=None, n_galfeeds=None,
                        galfit_timeout=60):

    logger = logging.getLogger("GalfitWorker")

    logger.debug("Galfit worker started")
    logger.debug("%s %s %s %s" % (str(n_galfit_queuesize), str(n_galfit_complete), str(n_total_galfit_time), str(n_galfeeds)))
    # return

    counter = 0
    while (True):

        cmd = galfit_queue.get()
        if (cmd is None):
            print("Received shutdown command")
            galfit_queue.task_done()
            break

        if (dryrun):
            print(cmd)
            galfit_queue.task_done()
            continue

        feedme_fn, galfit_output_fn, logfile = cmd
        _cwd, _feedfile = os.path.split(feedme_fn)

        # if (n_galfit_queuesize is not None):
        #     with n_galfit_queuesize.get_lock():
        #         n_galfit_queuesize.value -= 1

        counter += 1
        # print("get message %d - %d" % (counter, n_galfit_queuesize.value))
        # galfit_queue.task_done()
        # continue

        if (os.path.isfile(galfit_output_fn) and not redo):
            print("Skipping galfit run for completed file (%s)" % (galfit_output_fn))
            if (n_galfit_queuesize is not None):
                with n_galfit_queuesize.get_lock():
                    n_galfit_queuesize.value -= 1
            if (n_galfit_complete is not None):
                with n_galfit_complete.get_lock():
                    n_galfit_complete.value += 1
            galfit_queue.task_done()
            continue

        galfit_cmd = "%s %s" % (galfit_exe, _feedfile) #feedme_fn)
        # galfit_queue.task_done()
        # continue

        start_time = time.time()
        returncode = -99999999
        try:
            # os.system(sexcmd)
#             time.sleep(.2)
            with subprocess.Popen(galfit_cmd.split(),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  cwd=_cwd) as galfit_process:
                try:
                    _stdout, _stderr = galfit_process.communicate(input=None, timeout=galfit_timeout)


                    returncode = galfit_process.returncode
                    if (galfit_process.returncode != 0):
                        print("return code was not 0 (%d)" % (galfit_process.returncode))
                        print(str(_stdout))
                        print(str(_stderr))
                #print(_stdout)

                    with open(logfile, "wb") as log:
                        log.write(_stdout)
                        #log.write("\n*10STDERR\n========\n")
                        log.write(_stderr)

                except (TimeoutError, subprocess.TimeoutExpired) as e: #TimeoutExpired
                    galfit_process.kill()
                    print("Terminating galfit after timeout")
                    returncode = -9999999

            # ret = subprocess.Popen(galfit_cmd.split(),
            #                        stdout=subprocess.PIPE,
            #                        stderr=subprocess.PIPE)
            # (_stdout, _stderr) = ret.communicate()
            #     break

            # sex = subprocess.run(sexcmd.split(), shell=True, check=True,
            #                      stdout=subprocess.PIPE,
            #                      stderr=subprocess.PIPE,)
        except OSError as e:
            print("Some exception has occured:\n%s" % (str(e)))
        end_time = time.time()
        galfit_time = end_time - start_time
        # print("Galfit returned after %.3f seconds" % (end_time - start_time))
        # print(n_galfit_queuesize, n_galfit_complete, n_total_galfit_time)

        logger.debug("%s ==> %d" % (galfit_cmd, returncode))


        if (n_galfit_queuesize is not None):
            with n_galfit_queuesize.get_lock():
                n_galfit_queuesize.value -= 1
        if (n_galfit_complete is not None):
            with n_galfit_complete.get_lock():
                n_galfit_complete.value += 1
        if (n_total_galfit_time is not None):
            with n_total_galfit_time.get_lock():
                n_total_galfit_time.value += galfit_time
        # if (n_galfit_queuesize is not None and
        #         n_galfit_complete is not None and
        #         n_total_galfit_time is not None and
        #         n_galfeeds is not None):
        #     avg_galfit_time = n_total_galfit_time.value / n_galfit_complete.value
        #     print("Finished %d (of %d) galfit runs, %d (est. %.1f seconds) left" % (
        #         n_galfit_complete.value,
        #         n_galfeeds.value,
        #         n_galfit_queuesize.value,
        #         n_galfit_queuesize.value*avg_galfit_time,
        #     ))

        # if ((make_plots or True) and os.path.isfile(galfit_output_fn)):
        #     try:
        #         plot_galfit_results.plot_galfit_result(
        #             fits_fn=galfit_output_fn,
        #             plot_fn=galfit_output_fn[:-5]+".png",
        #         )
        #     except:
        #         print("Error while making plot")

        galfit_queue.task_done()
        continue

    print("Shutting down galfit parallel worker-process")



# def make_galfit_configs(image_fn, catalog_fn, config_basename, weight_image=None):
#
#
#     print (os.path.isfile(cat_fn))
#     catalog = numpy.loadtxt(cat_fn)
#     print(catalog.shape)
#
#     img_hdu = pyfits.open(image_fn)
#     magzero = 2.5*numpy.log10(img_hdu[0].header['FLUXMAG0'])
#
#     #
#     # Create the galfit feedme file
#     #
#     for src in catalog:
#
#         fwhm = src[9]
#         src_id = int(src[2])
#         size = 3 * fwhm
#         x, y = src[3], src[4]
#         x1 = numpy.max([0, x-size])
#         x2 = numpy.min([x+size])
#         y1 = y - size
#         y2 = y + size
#
#         galfit_info = {
#             'imgfile': img_fn,
#             'srcid': src_id,
#             'x1': x1,
#             'x2': x2,
#             'y1': y1,
#             'y2': y2,
#             'pixelscale': 0.18,
#             'weight_image': weight_image if weight_image is not None else 'none',
#         }
#
#         head_block = """
# A) %(imgfile)s         # Input data image (FITS file)
# B) output_%(srcid)d.fits    # Output data image block
# C) %(weight_image)s                # Sigma image name (made from data if blank or "none")
# D) psf.fits   #        # Input PSF image and (optional) diffusion kernel
# E) 1                   # PSF fine sampling factor relative to data
# F) none                # Bad pixel mask (FITS image or ASCII coord list)
# G) none                # File with parameter constraints (ASCII file)
# H) %(x1)d %(x2)d %(y1)d %(y2)d   # Image region to fit (xmin xmax ymin ymax)
# I) 100    100          # Size of the convolution box (x y)
# J) 26.000              # Magnitude photometric zeropoint
# K) %(pixelscale).3f %(pixelscale).3f            # Plate scale (dx dy)    [arcsec per pixel]
# O) regular             # Display type (regular, curses, both)
# P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
#
# """ % (galfit_info)
#         # print(head_block)
#
#         src_info = {
#             'x': x,
#             'y': y,
#             'magnitude': src[5], #+magzero,
#             'halflight_radius': src[22],
#             'sersic_n': 1.5, #src[7],
#             'axis_ratio': src[17]/src[16],
#             'position_angle': src[18],
#
#         }
#         object_block = """
# # Object number: 1
#  0) sersic                 #  object type
#  1) %(x)d  %(y)d  0 0  #  position x, y
#  3) %(magnitude).3f     1          #  Integrated magnitude
#  4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
#  5) %(sersic_n).3f      1          #  Sersic index n (de Vaucouleurs n=4)
#  6) 0.0000      0          #     -----
#  7) 0.0000      0          #     -----
#  8) 0.0000      0          #     -----
#  9) %(axis_ratio).3f      1          #  axis ratio (b/a)
# 10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
#  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
#
# # # Object number: 1
# #  0) devauc                 #  object type
# #  1) %(x)d  %(y)d  1 1  #  position x, y
# #  3) %(magnitude).3f     1          #  Integrated magnitude
# #  4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
# #  9) %(axis_ratio).3f      1          #  axis ratio (b/a)
# # 10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
# #  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
# #
# # # Object number: 1
# #  0) expdisk                 #  object type
# #  1) %(x)d  %(y)d  1 1  #  position x, y
# #  3) %(magnitude).3f     1          #  Integrated magnitude
# #  4) %(halflight_radius).3f      1          #  R_e (half-light radius)   [pix]
# #  9) %(axis_ratio).3f      1          #  axis ratio (b/a)
# # 10) %(position_angle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
# #  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
#
# # Object number: 2
#  0) sky                    #  object type
#  1) 0.0000      1          #  sky background at center of fitting region [ADUs]
#  2) 0.0000      0          #  dsky/dx (sky gradient in x)
#  3) 0.0000      0          #  dsky/dy (sky gradient in y)
#  Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
#
#
# """ % src_info
#         # print(object_block)
#
#         # feedme_fn = "feedme.%d" % (int(src[4]))
#         feedme_fn = "%s.src%05d.galfeed" % (config_basename, src_id)
#         with open(feedme_fn, "w") as feedfile:
#             feedfile.write(head_block)
#             feedfile.write(object_block)
#


if __name__ == "__main__":

    # setup command line parameters
    cmdline = argparse.ArgumentParser()
    # cmdline.add_argument("--nsims", dest="n_simulations",
    #                      default=1, type=int,
    #                      help="total number of simulated images (real+models) to generate from each input image")
    # cmdline.add_argument("--perimage", dest="n_per_image",
    #                      default=1, type=int,
    #                      help="number of models to insert into each image")
    # cmdline.add_argument("--mag", type=float, default=20.,
    #                      help="integrated magnitude of model")
    # cmdline.add_argument("--axisratio", dest="axis_ratio", type=float, default=1.0,
    #                      help="axis ratio")
    # cmdline.add_argument("--sersic", dest="sersic_index", type=float, default=1.0,
    #                      help="sersic index")
    # cmdline.add_argument("--sma", dest="halflight_semimajor_axis", type=float, default=25.,
    #                      help="half-light semi-major axis [pixels]")
    # cmdline.add_argument("--posangle", dest="position_angle", type=float, default=0.0,
    #                      help="position angle [degrees]")
    # cmdline.add_argument("--outdir", dest="output_directory", type=str, default="sims/",
    #                      help="output directory to hold simulated images")
    # cmdline.add_argument("--autodir", dest="auto_output_directory", default=False,
    #                      action='store_true',
    #                      help="auto-generate output directory to hold simulated images")

    cmdline.add_argument("--catext", dest="catalog_extension", type=str, default="udgcat",
                         help="file extension for catalog")
    cmdline.add_argument("--subdir", dest="galfit_directory", type=str, default="galfit/",
                         help="output subdirectory to hold galfit feed-files and output")
    cmdline.add_argument("--weight", dest="weight_file", type=str, default=None,
                         help="weight file")
    cmdline.add_argument("--nprocs", dest="number_processes",
                         default=multiprocessing.cpu_count(), type=int,
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("--galfit", dest="galfit_exe", type=str, default="galfit",
                         help="location of galfit executable")

    cmdline.add_argument("--plot", dest="plot_results", default=False,
                         action='store_true',
                         help="create plots from GALFIT results")

    cmdline.add_argument("--maxsize", dest="max_size", default=-1, type=int,
                         help="maximum cutout size for fitting")

    cmdline.add_argument("--psf", dest="psf", default=None, type=str,
                         help="filename of PSF model")
    cmdline.add_argument("--psfres", dest="psf_supersample", default=1, type=float,
                         help="super-sample factor of PSF model")

    cmdline.add_argument("--timeout", dest="galfit_timeout", default=60, type=float,
                         help="maximum tme allowed for a galfit run")

    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()

    print(args)

    # initialize work queues
    src_queue = multiprocessing.JoinableQueue()
    galfit_queue = multiprocessing.JoinableQueue()

    # start some counters for progress info and book-keeping
    n_galfeeds = multiprocessing.Value('i', 0, lock=True)#multiprocessing.Lock())
    n_galfit_queuesize = multiprocessing.Value('i', 0, lock=True)#lock=multiprocessing.Lock())
    n_galfit_complete = multiprocessing.Value('i', 0, lock=True)#lock=multiprocessing.Lock())
    n_total_galfit_time = multiprocessing.Value('f', 0.0, lock=True) #lock=multiprocessing.Lock())
    n_galfeeds.value = 0
    n_galfit_queuesize.value = 0
    n_galfit_complete.value = 0
    n_total_galfit_time.value = 0.


    total_feed_count = 0
    for fn in args.input_images:

        # open input image
        hdulist = pyfits.open(fn)

        try:
            magzero = 2.5 * numpy.log10(hdulist[0].header['FLUXMAG0'])
        except:
            magzero = 0

        # construct all filenames we need
        bn,_ = os.path.splitext(fn)
        #basename = os.path.split(bn)
        basedir,basename = os.path.split(bn)
        catalog_fn = "%s.%s" % (bn, args.catalog_extension)

        segmentation_fn = "%s.segments" % (bn)
        if (not os.path.isfile(segmentation_fn)):
            segmentation_fn = None

        # load catalog
        if (not os.path.isfile(catalog_fn)):
            print("Unable to open catalog %s" % (catalog_fn))
            continue

        # catalog = numpy.loadtxt(catalog_fn)
        try:
            catalog = astropy.table.Table.read(catalog_fn)
            print("done loading catalog")
        except:
            print("Unable to open catalog %s" % (catalog_fn))
            continue



        for src in catalog:
            src_id = int(src['NUMBER'])
            feedme_fn = "%s.%05d.galfeed" % (basename, src_id)

            if (args.galfit_directory.find(":") >= 0):
                _parts = args.galfit_directory.split(":")
                _search = _parts[0]
                _replace = _parts[1]
                galfit_dir = fn.replace(_search, _replace)
            else:
                galfit_dir = os.path.join(basedir, args.galfit_directory)

            if (not os.path.isdir(galfit_dir)):
                print("Creating directory: %s" % (galfit_dir))
                os.makedirs(galfit_dir)


            # also work out the appropriate weight filename
            if (args.weight_file.find(":") >= 0):
                _parts = args.weight_file.split(":")
                _search = _parts[0]
                _replace = _parts[1]
                weight_file = fn.replace(_search, _replace)
            else:
                weight_file = args.weight_file

            #
            # also construct the appropriate filename for the PSF
            #
            if (args.psf is not None):
                if (args.psf.find(":") >= 0):
                    # we need to do some replacing here
                    _parts = args.psf.split(":")
                    _search = _parts[0]
                    _replace = _parts[1]
                    psf_file = fn.replace(_search, _replace)
                else:
                    psf_file = args.psf
            else:
                psf_file = None

            #
            # Read PSF supersampling from PSF-file if available or default to 1
            #
            psf_supersample = 1.
            if (args.psf_supersample == 0 and psf_file is not None):
                try:
                    psf_hdu = pyfits.open(psf_file)
                    psf_supersample = 1./psf_hdu[0].header['PSF_SAMP']
                    psf_hdu.close()
                except Exception as e:
                    print(e)
                    pass

            feedme_fullfn = os.path.join(galfit_dir, feedme_fn)
            # print(feedme_fullfn)

            # _segmentation_cutout_fn = "%s.%05d.segment" % (basename, src_id)
            # segmentation_cutout_fn = os.path.join(galfit_dir, _segmentation_cutout_fn)
            #
            # _segmenation_input_fn =

            galfit_fn = "%s.%05d.galfit.fits" % (basename, src_id)
            galfit_fullfn = os.path.join(galfit_dir, galfit_fn)


            galfit_logfn = "%s.%05d.galfit.log" % (basename, src_id)
            galfit_fulllogfn = os.path.join(galfit_dir, galfit_logfn)

            src_queue.put((
                fn,                 # image_fn,
                feedme_fullfn,      # feedme_fn,
                weight_file,        # weight_fn,
                src,                # src_info,
                hdulist[0].header,  # img_header
                galfit_fullfn,      # galfit output filename
                galfit_fulllogfn,
                segmentation_fn,
                galfit_dir,
                basename,
                args.max_size,
                psf_file, psf_supersample,
                magzero,
            ))
            # print(src)

            total_feed_count += 1

    print("Configuring %d galfit runs using %d CPU-cores" % (total_feed_count, args.number_processes))
    # sys.exit(0)

    #
    # Now that we have filled up the work-queue, get to work
    #
    feedme_workers = []
    for i in range(args.number_processes):
        p = multiprocessing.Process(
            target=parallel_config_writer,
            kwargs=dict(queue=src_queue,
                        galfit_queue=galfit_queue,
                        n_galfeeds=n_galfeeds,
                        n_galfit_queuesize=n_galfit_queuesize,
                        ),
        )
        p.daemon = True
        p.start()
        feedme_workers.append(p)
        src_queue.put((None))
    print("Feed-me creation workers started")
    src_queue.join()
    print("Done creating all GALFIT config files, joining the fun!")
    # print(n_galfeeds, n_galfit_queuesize)
    # time.sleep(1)

    # counter = 0
    # while(True):
    #     cmd = galfit_queue.get()
    #     galfit_queue.task_done()
    #     if (cmd is not None):
    #         counter += 1
    #     print(counter, n_galfeeds.value)



    galfit_workers = []
    for i in range(args.number_processes):
        p = multiprocessing.Process(
            target=parallel_run_galfit,
            kwargs=dict(galfit_queue=galfit_queue,
                        galfit_exe=args.galfit_exe,
                        make_plots=args.plot_results,
                        redo=False,
                        n_galfit_complete=n_galfit_complete,
                        n_total_galfit_time=n_total_galfit_time,
                        n_galfit_queuesize=n_galfit_queuesize,
                        n_galfeeds=n_galfeeds,
                        galfit_timeout=args.galfit_timeout,
                        )
        )
        p.daemon = True
        p.start()
        galfit_workers.append(p)
        # galfit_queue.put((None))

    print("Galfit workers started")

    # wait for all config files to be written
    # for i in range(30):
    #     print n_galfit_queuesize.lock
    #     time.sleep(10)

    while (n_galfit_queuesize.value > 0):
        try:
            avg_galfit_time = n_total_galfit_time.value / n_galfit_complete.value
        except ZeroDivisionError:
            avg_galfit_time = 10.
        sys.stdout.write("\rFinished %d (of %d) galfit runs, %d (est. %.1f seconds) left" % (
            n_galfit_complete.value,
            n_galfeeds.value,
            n_galfit_queuesize.value,
            n_galfit_queuesize.value * avg_galfit_time,
        ))
        sys.stdout.flush()
        time.sleep(1)

    # galfit_queue.join()
    print("done with all work!")

    # img_fn = sys.argv[1]
    # cat_fn = sys.argv[2]
    #
    # weight_image = None
    # try:
    #     weight_image = sys.argv[3]
    # except IndexError:
    #     pass
    #
    # bn,_ = os.path.splitext(cat_fn)
    # make_galfit_configs(image_fn=img_fn, catalog_fn=cat_fn, config_basename=bn,
    #                     weight_image=weight_image)
