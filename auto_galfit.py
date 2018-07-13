#!/usr/bin/env python3

import os
import sys
import numpy
import pyfits

import argparse
import multiprocessing
import time
import subprocess


import plot_galfit_results


def parallel_config_writer(queue, galfit_queue):

    while (True):

        cmd = queue.get()
        if (cmd is None):
            queue.task_done()
            galfit_queue.put((None))
            break

        image_fn, feedme_fn, weight_fn, src_info, img_header, galfit_output, galfit_logfile, segmentation_input_fn, galfit_dir, basename = cmd
        if (os.path.isfile(feedme_fn)):
            queue.task_done()
            galfit_queue.put((feedme_fn, galfit_output, galfit_logfile))
            continue

        magzero = 2.5 * numpy.log10(img_header['FLUXMAG0'])

        #
        # Create the galfit feedme file, the cutouts for the
        # image, weight-image, and segmentation mask
        #
        fwhm = src_info[9]
        src_id = int(src_info[2])
        size = 3 * fwhm
        x, y = src_info[3]-1, src_info[4]-1
        x1 = int(numpy.max([0, x - size]))
        x2 = int(numpy.min([x + size, img_header['NAXIS1']]))
        y1 = int(numpy.max([0, y - size]))
        y2 = int(numpy.min([y + size, img_header['NAXIS2']]))

        # open the input images and create cutouts
        img_out_fn = "%s/%s.%05d.image.fits" % (galfit_dir, basename, src_id)
        segm_out_fn = "%s/%s.%05d.segm.fits" % (galfit_dir, basename, src_id)
        weight_out_fn = "%s/%s.%05d.sigma.fits" % (galfit_dir, basename, src_id)

        img_hdu = pyfits.open(image_fn)
        img = img_hdu[0].data[y1:y2, x1:x2]
        phdu = pyfits.PrimaryHDU(data=img)
        phdu.header['SRC_X1'] = x1
        phdu.header['SRC_Y1'] = y1
        phdu.writeto(img_out_fn)
        img_hdu.close()

        if (weight_fn is not None):
            wht_hdu = pyfits.open(weight_fn)
            wht = wht_hdu[0].data[y1:y2, x1:x2]
            pyfits.PrimaryHDU(data=wht, header=phdu.header).writeto(weight_out_fn)
            wht_hdu.close()

        if (segmentation_fn is not None):
            segm_hdu = pyfits.open(segmentation_fn)
            segm = segm_hdu[0].data[y1:y2, x1:x2].astype(numpy.int)
            segm[segm == src_id] = 0
            pyfits.PrimaryHDU(data=segm, header=phdu.header).writeto(segm_out_fn)
            segm_hdu.close()


        galfit_info = {
            'imgfile': img_out_fn, #image_fn,
            'srcid': src_id,
            'x1': 0, #x1,
            'x2': x2-x1, #x2
            'y1': 0, #y1,
            'y2': y2-y1, #y2,
            'pixelscale': 0.18,
            'weight_image': weight_out_fn if weight_fn is not None else 'none',
            'galfit_output': galfit_output,
            'bpm': segm_out_fn if segmentation_fn is not None else 'none',
        }

        head_block = """
            A) %(imgfile)s         # Input data image (FITS file)
            B) %(galfit_output)s   # Output data image block
            C) %(weight_image)s                # Sigma image name (made from data if blank or "none") 
            D) psf.fits   #        # Input PSF image and (optional) diffusion kernel
            E) 1                   # PSF fine sampling factor relative to data 
            F) %(bpm)s                # Bad pixel mask (FITS image or ASCII coord list)
            G) none                # File with parameter constraints (ASCII file) 
            H) %(x1)d %(x2)d %(y1)d %(y2)d   # Image region to fit (xmin xmax ymin ymax)
            I) 100    100          # Size of the convolution box (x y)
            J) 26.000              # Magnitude photometric zeropoint 
            K) %(pixelscale).3f %(pixelscale).3f            # Plate scale (dx dy)    [arcsec per pixel]
            O) regular             # Display type (regular, curses, both)
            P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        """ % (galfit_info)
            # print(head_block)

        src_info = {
            'x': x-x1,
            'y': y-y1,
            'magnitude': src_info[5],  # +magzero,
            'halflight_radius': src_info[22],
            'sersic_n': 1.5,  # src[7],
            'axis_ratio': src_info[17] / src_info[16],
            'position_angle': src_info[18],

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


        galfit_queue.put((feedme_fn, galfit_output, galfit_logfile))

        queue.task_done()
        continue


dryrun = False


def parallel_run_galfit(galfit_queue, galfit_exe='galfit', make_plots=True, redo=False):

    print("Galfit worker started")

    while (True):

        cmd = galfit_queue.get()
        if (cmd is None):
            galfit_queue.task_done()
            break

        if (dryrun):
            print(cmd)
            galfit_queue.task_done()
            continue

        feedme_fn, galfit_output_fn, logfile = cmd

        if (os.path.isfile(galfit_output_fn) and not redo):
            galfit_queue.task_done()
            continue

        galfit_cmd = "%s %s" % (galfit_exe, feedme_fn)
        print(galfit_cmd)
        # galfit_queue.task_done()
        # continue

        start_time = time.time()
        try:
            # os.system(sexcmd)
            ret = subprocess.Popen(galfit_cmd.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
            (_stdout, _stderr) = ret.communicate()
            if (ret.returncode != 0):
                print("return code was not 0")
                print(_stdout)
                print(_stderr)
            #print(_stdout)

            with open(logfile, "wb") as log:
                log.write(_stdout)
                #log.write("\n*10STDERR\n========\n")
                log.write(_stderr)

            #     break

            # sex = subprocess.run(sexcmd.split(), shell=True, check=True,
            #                      stdout=subprocess.PIPE,
            #                      stderr=subprocess.PIPE,)
        except OSError as e:
            print("Some exception has occured:\n%s" % (str(e)))
        end_time = time.time()
        print("Galfit returned after %.3f seconds" % (end_time - start_time))


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

    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    #cmdline.print_help()
    args = cmdline.parse_args()


    src_queue = multiprocessing.JoinableQueue()
    galfit_queue = multiprocessing.JoinableQueue()

    total_feed_count = 0
    for fn in args.input_images:

        # open input image
        hdulist = pyfits.open(fn)

        # construct all filenames we need
        bn,_ = os.path.splitext(fn)
        #basename = os.path.split(bn)
        basedir,basename = os.path.split(bn)
        catalog_fn = "%s.%s" % (bn, args.catalog_extension)

        segmentation_fn = "%s.segments" % (bn)
        if (not os.path.isfile(segmentation_fn)):
            segmentation_fn = None

        # load catalog
        catalog = numpy.loadtxt(catalog_fn)

        for src in catalog:
            src_id = int(src[2])
            feedme_fn = "%s.%05d.galfeed" % (basename, src_id)
            galfit_dir = os.path.join(basedir, args.galfit_directory)
            if (not os.path.isdir(galfit_dir)):
                print("Creating directory: %s" % (galfit_dir))
                os.makedirs(galfit_dir)
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
                args.weight_file,   # weight_fn,
                src,                # src_info,
                hdulist[0].header,  # img_header
                galfit_fullfn,      # galfit output filename
                galfit_fulllogfn,
                segmentation_fn,
                galfit_dir, basename,
            ))
            # print(src)

            total_feed_count += 1

    print("Configuring %d galfit runs using %d CPU-cores" % (total_feed_count, args.number_processes))

    #
    # Now that we have filled up the work-queue, get to work
    #
    feedme_workers = []
    for i in range(8): #args.number_processes):
        p = multiprocessing.Process(
            target=parallel_config_writer,
            kwargs=dict(queue=src_queue,
                        galfit_queue=galfit_queue,
                        ),
        )
        p.daemon = True
        p.start()
        feedme_workers.append(p)
        src_queue.put((None))
    print("Feed-me creation workers started")
    src_queue.join()
    print("Done creating all GALFIT config files, joining the fun!")

    galfit_workers = []
    for i in range(args.number_processes):
        p = multiprocessing.Process(
            target=parallel_run_galfit,
            kwargs=dict(galfit_queue=galfit_queue,
                        galfit_exe=args.galfit_exe,
                        make_plots=args.plot_results,
                        redo=False,
                        )
        )
        p.daemon = True
        p.start()
        galfit_workers.append(p)
    print("Galfit workers started")

    # wait for all config files to be written

    galfit_queue.join()
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