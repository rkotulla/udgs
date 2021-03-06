#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import time
import itertools
import pandas

import astropy
import astropy.table
import astropy.io.fits as pyfits
import numpy
import distutils.spawn
import shutil

import run_sextractor

import ldac2vot

fix_vot_array_data = ldac2vot.read_definitions("FLUX_RADIUS:50,80")


def range_to_list(arg):

    try:
        val = float(arg)
        return numpy.array([val]),0
    except ValueError:
        # this means it's not a simple number
        pass

    # check if it in the format a..b:s
    scatter = 0
    if (arg.find("..")>0 and arg.find(":")):
        # we got this format
        try:
            a = float(arg.split("..")[0])
            b = float(arg.split("..")[1].split(":")[0])
            s = float(arg.split(":")[1])
            try:
                scatter = float(arg.split(":")[2])
            except:
                scatter = 0
            # print(a,b,s)
            return numpy.arange(a, b+s, s), scatter
        except ValueError:
            print("Unable to understand format of %s" % (arg))
            return None
    elif (arg.find(",")):
        # this is a list of values
        try:
            l = [float(f) for f in arg.split(",")]
            return numpy.array(l), scatter
        except ValueError:
            print("Illegal format or value in %s" % (arg))

    return None


def set_or_replace(input_fn, param):

    if (param.find(":") > 0):
        # This is a search & replace
        find = param.split(":")[0]
        replace = param.split(":")[1]
        out = input_fn.replace(find, replace)
    else:
        out = param
    return out


def completeness_worker(file_queue, src_img, img_size, psf_file,
                        output_dir, singles_dir, singles_size,
                        galfit_exe,
                        sex_exe, sex_conf, sex_param,
                        weight_fn,
                        ):

    print("Worker started")
    # print(src_img.info())

    #
    # Find out some basics about the image and the PSF model
    #
    psf_hdu = pyfits.open(psf_file)
    psf_sampling = psf_hdu[0].header['SUPERSMP']

    img_hdr = src_img[0].header
    img_data = src_img[0].data
    pixelscale = img_hdr['CD2_2'] * 3600.
    magzero = 2.5*numpy.log10(img_hdr['FLUXMAG0']) if 'FLUXMAG0' in img_hdr else 27.0
    print("Pixelscale:", pixelscale, ", PSF-sampling:", psf_sampling)

    (img_x, img_y) = img_size

    minisize = 2*(singles_size//2)+1
    halfsize = (minisize-2)//2
    while (True):

        job = file_queue.get()
        if (job is None):
            file_queue.task_done()
            break

        # Generate one full-frame image to receive all individual model images
        model_assembly = numpy.zeros_like(src_img[0].data)
        print(model_assembly.shape)

        chunk_id = job['chunk_id']

        ###########
        #
        # Generate all individual models
        #
        ###########
        total_galfit_time = 0.
        for i, src in job['sources'].iterrows():
            # print(src['id'])


            #
            # Work out how we can fit each individual model image into the
            # full frame, making sure to truncate frames along the edges
            #
            x, y = src['cx']-1, src['cy']-1
            x1 = int(numpy.max([0, x - halfsize]))
            x2 = int(numpy.min([x + halfsize, img_hdr['NAXIS1']]))
            y1 = int(numpy.max([0, y - halfsize]))
            y2 = int(numpy.min([y + halfsize, img_hdr['NAXIS2']]))

            #
            # Generate all filenames we will need
            #
            feedme_fn = os.path.join(singles_dir, "model_%06d.feedme" % (src['id']))
            model_only_fn = os.path.join(singles_dir, "model_%06d.raw.fits" % (src['id']))
            galfit_logfile = "model_%06d.galfit.log" % (src['id'])

            _, feedme_bn = os.path.split(feedme_fn)
            _,_model_only_fn = os.path.split(model_only_fn)
            _, psf_bn = os.path.split(psf_file)

            # print("Generating feed-me file: %s" % (feedme_fn))
            # job['xxx'].info()


            ff = open(feedme_fn, "w")
            ff_header = """
                A) 
                B) %(model_only_fn)s   # Output data image block
                C)                 # Sigma image name (made from data if blank or "none")
                D) %(psf_file)s   #        # Input PSF image and (optional) diffusion kernel
                E) %(psf_sampling)d                   # PSF fine sampling factor relative to data
                F)                 # Bad pixel mask (FITS image or ASCII coord list)
                G)                 # File with parameter constraints (ASCII file)
                H) 0 %(img_x)d 0 %(img_y)d   # Image region to fit (xmin xmax ymin ymax)
                I) 200    200          # Size of the convolution box (x y)
                J) %(magzero).3f              # Magnitude photometric zeropoint
                K) %(dx).3f %(dy).3f            # Plate scale (dx dy)    [arcsec per pixel]
                O) regular             # Display type (regular, curses, both)
                P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
            """ % dict(
                model_only_fn=_model_only_fn,
                psf_file=psf_bn,
                psf_sampling=psf_sampling,
                dx=pixelscale, dy=pixelscale,
                img_x=(x2-x1), #minisize,
                img_y=(y2-y1), #minisize,
                magzero=magzero,
            )
            ff.write("\n".join([l.strip() for l in ff_header.splitlines(keepends=False)]))

            #
            # Now generate all the model definitions
            #
            # params = job['params']
            # positions = job['positions']
            # n_galaxies = params.shape[0]
            # for gal in range(n_galaxies):
            src_def = """
            
                # Object number: %(i)d
                0) sersic                 #  object type
                1) %(x).3f %(y).3f  1 1  #  position x, y
                3) %(mag).3f     1          #  Integrated magnitude
                4) %(reff).3f      1          #  R_e (half-light radius)   [pix]
                5) %(sersic).3f      1          #  Sersic index n (de Vaucouleurs n=4)
                6) 0.0000      0          #     -----
                7) 0.0000      0          #     -----
                8) 0.0000      0          #     -----
                9) %(axisratio).3f      1          #  axis ratio (b/a)
                10) %(posangle).3f    1          #  position angle (PA) [deg: Up=0, Left=90]
                Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
                
                ###-X1: %(x1)d
                ###-Y1: %(y1)d
            """ % dict(
                i=src['id'],
                x=src['cx']-x1, #'(minisize-1)//2, #positions[gal,0],
                y=src['cy']-y1, #'(minisize-1)//2, #positions[gal,1],
                mag=src['final_mag'],
                reff=src['final_r_eff'],
                sersic=src['final_sersic'],
                axisratio=src['final_axisratio'],
                posangle=src['final_posangle'],
                x1=x1, y1=y1,
            )
            ff.write("\n".join(l.strip() for l in src_def.splitlines(keepends=False)))

            ff.close()

            ################################
            #
            # Now we have everything in place we need to run galfit and
            # create the simulated model image
            #
            ################################

            galfit_cmd = "%s %s" % (galfit_exe, feedme_bn) #feedme_fn)
            galfit_timeout = 300

            start_time = time.time()
            returncode = -99999999
            try:
                with subprocess.Popen(galfit_cmd.split(),
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      cwd=singles_dir) as galfit_process:
                    try:
                        _stdout, _stderr = galfit_process.communicate(
                            input=None, timeout=galfit_timeout)

                        returncode = galfit_process.returncode
                        if (galfit_process.returncode != 0):
                            print("return code was not 0 (%d)" % (
                                galfit_process.returncode))
                            print(str(_stdout))
                            print(str(_stderr))
                        # print(_stdout)

                        with open(galfit_logfile, "wb") as log:
                            log.write(_stdout)
                            # log.write("\n*10STDERR\n========\n")
                            log.write(_stderr)

                    except (TimeoutError,
                            subprocess.TimeoutExpired) as e:  # TimeoutExpired
                        galfit_process.kill()
                        print("Terminating galfit after timeout")
                        returncode = -9999999

                        # problem_str = "%s ::: %s\n" % (
                        # feedme_fn, " ".join(galfit_cmd.split()))
                        # problems_queue.put(problem_str)

            except OSError as e:
                print("Some exception has occured:\n%s" % (str(e)))
            end_time = time.time()
            galfit_time = end_time - start_time
            total_galfit_time += galfit_time
            # print("Galfit returned after %.3f seconds" % (galfit_time))
            # print(n_galfit_queuesize, n_galfit_complete, n_total_galfit_time)

            try:
                model_hdu = pyfits.open(model_only_fn)
            except (FileNotFoundError, OSError) as e:
                print("Bad model file:", e)
            model_img = model_hdu[0].data
            model_assembly[y1:y2, x1:x2] += model_img
            model_hdu.close()

            # logger.debug("%s ==> %d" % (galfit_cmd, returncode))

        print("Generated %d artificial galaxies in %.2f seconds" % (
            len(job['sources'].index), total_galfit_time))
        ###############################
        #
        # Now we have all models in one batch completed
        #
        ###############################
        chunk_models_fn = os.path.join(output_dir, "modelchunk_%05d.raw.fits" % (chunk_id))
        pyfits.PrimaryHDU(data=model_assembly).writeto(chunk_models_fn, overwrite=True)

        #
        # Also add the model images to the actual input images to generate
        # the full simulated data used for completeness analysis
        #
        comp_image_fn = os.path.join(output_dir, "modelchunk_%05d.fits" % (chunk_id))
        completeness_input = img_data + model_assembly
        completeness_hdu = pyfits.PrimaryHDU(
            data=completeness_input,
            header=img_hdr,
        )
        completeness_hdu.writeto(comp_image_fn, overwrite=True)


        ###############################
        #
        # TODO: implement...
        # Next up: Run SourceExtractor on the newly generated simulated image,
        # load the resulting catalog and find for sources that are not part of
        # the reference catalog generated from only the input image
        #
        ###############################
        raw_catalog = run_sextractor.run_sex(
            file_queue=None,
            sex_exe=sex_exe,
            sex_conf=sex_conf,
            sex_param=sex_param,
            fix_vot_array=fix_vot_array_data,
            single_frame=(comp_image_fn, weight_fn)
        )
        # raw_catalog.info()


        # TODO: write script to convert FITS catalogs to some format
        #  ds9 can handle to make testing etc easier, but that doesn't take
        #  forever to read & write

        # print(job['params'])
        file_queue.task_done()

    print("Worker shutting down")



if __name__ == "__main__":

    # get path of this executable
    fn = os.path.abspath(__file__)
    dirname,_ = os.path.split(fn)
    config_dir = os.path.join(dirname, "config")

    # setup command line parameters
    cmdline = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    cmdline.add_argument("--conf", dest="sex_conf", default=os.path.join(config_dir, "sex.conf"),
                         help="source extractor config filename")
    cmdline.add_argument("--params", dest="sex_params", default=os.path.join(config_dir, "sex.param"),
                         help="number of models to insert into each image")
    cmdline.add_argument("--nprocs", dest="number_processes", type=int,
                         default=4, #multiprocessing.cpu_count(),
                         help="number of Sextractors to run in parallel")
    cmdline.add_argument("--sex", dest="sex_exe", default=distutils.spawn.find_executable("sex"),
                         help="location of SExtractor executable")
    cmdline.add_argument("--galfit", dest="galfit_exe", default=distutils.spawn.find_executable("galfit"),
                         help="location ofGalfit executable")
    cmdline.add_argument("--weight", dest='weight_image', type=str, default=None,
                         help="weight map")
    cmdline.add_argument("--psf", dest='psf_image', type=str, default="_image.fits:_psf.fits",
                         help="psf basename")

    cmdline.add_argument("--refcat", dest='ref_cat', default="_image.fits:_image.vot",
                         help="reference catalog")


    cmdline.add_argument("--mag", dest='mag', default="21..24:0.5",
                         help="magnitude-range, from-to:step")
    cmdline.add_argument("--re", dest='r_eff', default="5..50:5",
                         help="effective radius range, from-to:step")
    cmdline.add_argument("--ar", dest='axisratio', default="0.2..1:0.2",
                         help="axis ratio range, from-to:step")
    cmdline.add_argument("--sersic", dest='sersic', default="0.5..4:0.5",
                         help="sersic index range, from-to:step")
    cmdline.add_argument("--pa", dest='posangle', default="0..180:30:30",
                         help="position angle range, from-to:step")
    cmdline.add_argument("--imgmargin", dest='img_margin', default=100, type=float,
                         help="position angle range, from-to:step")

    cmdline.add_argument("-n", "--nmodels", dest='n_models', default=10, type=int,
                         help="number of models for each configuration")
    cmdline.add_argument("-p", "--perframe", dest='models_per_frame', default=25, type=int,
                         help="number of models, across configurations, to add simultaneously")

    cmdline.add_argument("--dir", dest='dirname', default="_image.fits:_completeness/", type=str,
                         help="filename of completeness log")
    cmdline.add_argument("--logfile", dest='logfile', default="_image.fits:_completeness/log.fits", type=str,
                         help="filename of completeness log")

    cmdline.add_argument("input_images", nargs="+",
                         help="list of input images")
    # cmdline.print_help()
    args = cmdline.parse_args()

    construct_weight_fn = False
    if (args.weight_image is not None):
        construct_weight_fn = (args.weight_image.find(":") > 0)

    # First, lets work out all the parameters to iterate over to get an
    # idea what models we need
    model_mags, scatter_mags = range_to_list(args.mag)
    model_radius, scatter_radius = range_to_list(args.r_eff)
    model_axisratio, scatter_axisratio = range_to_list(args.axisratio)
    model_sersic, scatter_sersic = range_to_list(args.sersic)
    model_posangle, scatter_posangle = range_to_list(args.posangle)

    print("mags:       ", model_mags)
    print("eff radius: ", model_radius)
    print("b/a:        ", model_axisratio)
    print("sersic:     ", model_sersic)
    print("pos-angle:  ", model_posangle)

    total_number_of_models = \
        model_mags.shape[0] \
        * model_radius.shape[0] \
        * model_axisratio.shape[0] \
        * model_sersic.shape[0] \
        * model_posangle.shape[0] \
        * args.n_models

    print("total # of models: ", total_number_of_models)
    print("# of frames:", total_number_of_models//args.models_per_frame)

    # now generate the full list of model options
    params = itertools.product(model_mags, model_radius, model_axisratio, model_sersic, model_posangle)
    # print(list(params))
    full_params = numpy.array(list(params) * args.n_models)
    # print(full_params.shape)

    pandas_columns = ['mags', 'r_eff', 'axisratio', 'sersic', 'posangle',
                      'cx', 'cy',
                      'd_mag', 'd_r_reff', 'd_axisratio', 'd_sersic', 'd_posangle',
                      'final_mag', 'final_r_eff', 'final_axisratio', 'final_sersic', 'final_posangle',
    ]

    pandas_params = pandas.DataFrame(full_params, columns=['mags', 'r_eff', 'axisratio', 'sersic', 'posangle'])
    # pandas_params.info()

    for filename in args.input_images:

        model_parameters = full_params.copy()

        #
        # Open the input file and get some basic information
        #
        hdulist = pyfits.open(filename)
        hdr = hdulist[0].header

        img_x = hdr['NAXIS1']
        img_y = hdr['NAXIS2']
        magzero = -2.5*numpy.log10(hdr['FLUXMAG0']) if 'FLUXMAG0' in hdr else 27.0

        positions = (numpy.random.random((full_params.shape[0], 2)) * \
                     [img_x-2*args.img_margin, img_y-2*args.img_margin]) + [args.img_margin, args.img_margin]


        # get output directory name
        output_dirname = set_or_replace(filename, args.dirname)
        singles_dirname = os.path.join(output_dirname, "singles")
        if (not os.path.isdir(output_dirname)):
            os.makedirs(output_dirname)
        if (not os.path.isdir(singles_dirname)):
            os.makedirs(singles_dirname)

        #
        # Figure out the other files needed for this one
        #
        psf_file = set_or_replace(filename, args.psf_image)
        weight_fn = set_or_replace(filename, args.weight_image)

        #
        # also copy the PSF image to the singles directory
        _,bn = os.path.split(psf_file)
        shutil.copy(psf_file, os.path.join(singles_dirname, bn))

        #
        # Now add the scatter to each of the datapoints
        #
        param_scatter = (numpy.random.random(size=full_params.shape) - 0.5) * \
                     [scatter_mags, scatter_radius, scatter_axisratio, scatter_sersic, scatter_posangle]
        # print(param_scatter[:20])

        final_params = full_params + param_scatter

        # pandas_scatter = pandas.DataFrame(param_scatter, columns=['d_mag', 'd_r_reff', 'd_axisratio', 'd_sersic', 'd_posangle'])
        # model_parameters.append(pandas_scatter)

        #
        # Prepare some shuffeling to get a more even sampling of model galaxies
        # in each of the fake frames
        #
        permutater = numpy.random.permutation(full_params.shape[0])
        # print(permutater[:25])
        # combine all parameters into a single array
        param_data = numpy.hstack((full_params, positions, param_scatter, final_params))
        # print(param_data.shape)
        #
        # apply shuffeling
        param_data_shuffled = param_data[permutater]

        #
        # Now convert to pandas array to make handling columns easier
        #
        pandas_params = pandas.DataFrame(param_data_shuffled, columns=pandas_columns)
        pandas_params['id'] = numpy.arange(param_data_shuffled.shape[0], dtype=numpy.int)
        # pandas_params.info()

        #
        # add unique ID to each source

        #
        # Before we do anything else, we need to record what we inserted and where
        #
        model_log_filename = set_or_replace(filename, args.logfile)
        print("Writing completeness log to %s" % (model_log_filename))
        # pandas_params.to_csv(model_log_filename)
        ap_params = astropy.table.Table.from_pandas(pandas_params)
        # ap_params.info()
        ap_params.write(model_log_filename, format='fits', overwrite=True)

        model_params = (full_params + param_scatter)[permutater]

        #
        # Start the completeness workers before handing out work
        # that way they can get started right away
        #
        file_queue = multiprocessing.JoinableQueue()
        processes = []
        for i in range(args.number_processes):
            p = multiprocessing.Process(
                target=completeness_worker,
                kwargs=dict(
                    file_queue=file_queue,
                    src_img=hdulist,
                    img_size=(img_x,img_y),
                    psf_file=psf_file,
                    output_dir=output_dirname,
                    singles_dir=singles_dirname,
                    singles_size=300,
                    galfit_exe=args.galfit_exe,
                    sex_exe=args.sex_exe,
                    sex_conf=args.sex_conf,
                    sex_param=args.sex_params,
                    weight_fn=weight_fn,
                )
            )
            p.daemon = True
            p.start()
            processes.append(p)

        #
        # Split the sample into chunks based on the number of galaxies to be
        # inserted into each original input frame
        #
        chunksize = args.models_per_frame
        n_chunks = int(numpy.ceil(total_number_of_models / chunksize))
        n_chunks=4
        for chunk in range(n_chunks):

            file_queue.put(dict(
                positions=positions[chunk*chunksize:(chunk+1)*chunksize],
                params=model_params[chunk*chunksize:(chunk+1)*chunksize],
                chunk_id=chunk,
                sources=pandas_params[chunk*chunksize:(chunk+1)*chunksize],
            ))




        #
        # Now wait for all work to be completed before going on to the
        # next input frame
        #

        # insert termination commands
        for i in range(args.number_processes):
            file_queue.put((None))
        # now wait for all work to be done
        file_queue.join()


