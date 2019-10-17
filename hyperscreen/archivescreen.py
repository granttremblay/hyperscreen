#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function

from hyperscreen import hypercore
import os
import sys
import time
import glob
import argparse


import multiprocessing
import pickle
from functools import partial

from astropy.io import fits

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")


def reportCard(evt1_object, savepath, show=True, save=True, rasterized=True, dpi=150, verbose=False):  # pragma: no cover

    obs = evt1_object

    if verbose is True:
        print("Doing {}, {}".format(obs.obsid, obs.detector))

    reportCard_savepath = os.path.join(
        savepath, '{}_{}_{}_hyperReport.pdf'.format(obs.obsid, obs.target.replace(' ', '_'), obs.detector))

    with PdfPages(reportCard_savepath) as pdf:

        # MAKE PAGE 1

        fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey='row')

        obs.boomerang(mask=obs.data['Hyperbola test passed'], ax=axes[0, 0], create_subplot=True,
                      show=False, title='Legacy Hyperbola Test', cmap='magma', rasterized=rasterized)
        obs.boomerang(ax=axes[0, 1], create_subplot=True,
                      show=False, title='Test2', cmap='inferno', rasterized=rasterized)

        obs.image(ax=axes[1, 0], detcoords=True, show=False,
                  create_subplot=True, title="Test1", rasterized=rasterized)
        obs.image(ax=axes[1, 1], detcoords=True, show=False,
                  create_subplot=True, title="Test2", rasterized=rasterized)

        fig.suptitle('ObsID {} | {} | {}'.format(
            obs.obsid, obs.target, obs.detector))

        if save is True:
            pdf.savefig(fig)

        # MAKE PAGE 2

        fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey='row')

        obs.boomerang(mask=obs.data['Hyperbola test passed'], ax=axes[0, 0], create_subplot=True,
                      show=False, title='Legacy Hyperbola Test', cmap='magma', rasterized=rasterized)
        obs.boomerang(ax=axes[0, 1], create_subplot=True,
                      show=False, title='Test2', cmap='inferno', rasterized=rasterized)

        obs.image(masked_x=obs.data['detx'][obs.data['Hyperbola test failed']], masked_y=obs.data['dety'][obs.data['Hyperbola test failed']], ax=axes[1, 0], detcoords=True, show=False,
                  create_subplot=True, title="Test1", rasterized=rasterized)
        obs.image(ax=axes[1, 1], detcoords=True, show=False,
                  create_subplot=True, title="Test2", rasterized=rasterized)

        fig.suptitle('ObsID {} | {} | {}'.format(
            obs.obsid, obs.target, obs.detector))

        if save is True:
            pdf.savefig(fig)

        if verbose is True:
            print("Created {}".format(reportCard_savepath))

        if show is True:
            plt.show()

    plt.close()


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('-a', '--archivepath', help='Absolute PATH to Archive of EVT1 Files',
                        default='/Users/grant/Science/HRC_Database/EVT1_Files/')

    parser.add_argument('-c', '--cluster', action='store_true',
                        help='Point to the HRC Database stored on the Smithsonian Hydra Cluster')

    parser.add_argument('-p', '--picklename', help='Name of the Pickle you would like to create',
                        default='hyperscreen_master_pickle.pkl')

    parser.add_argument('-r', '--reportcard', help='Make Report Card .pdf files while screening Archive? Defaults to True.',
                        default=True)

    parser.add_argument('-s', '--savepath', help='Absolute PATH to location in which to save all outputs from this script, including .pdf files of plots. If not specified, this location will default to your Desktop.',
                        default=os.path.join(os.environ['HOME'], 'Desktop/HyperScreen_Results/'))

    # parser.add_argument('--showplots', help='Absolute PATH to location in which to save all outputs from this script, including .pdf files of plots. If not specified, this location will default to your Desktop.',
    #                     default=os.path.join(os.environ['HOME'], 'Desktop/HyperScreen_Results/'))

    parser.add_argument('-t', '--testdata', action='store_true',
                        help='Use the supplied test data as the input archive path')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Make HyperScreen chatty on stdout.')

    parser.add_argument('-w', '--windowstest', action='store_true',
                        help='Point to my Windows database')

    return parser.parse_args(argv)


def setPaths(args, verbose=False):  # pragma: no cover
    '''
    Set paths for this hyperscreen test
    '''

    savepath = args.savepath
    if not os.path.exists(savepath):
        print("Creating Directory {} in which HyperScreen results will be saved.".format(
            savepath))
        os.makedirs(savepath)
    if verbose is True:
        print("Plots will be saved in {}".format(savepath))

    if args.testdata is True:
        archivepath = os.path.abspath(os.path.dirname(
            os.path.abspath(__file__)) + '/../tests/data/')
    elif args.windowstest is True:
        archivepath = '/mnt/c/Users/grant/HRCOps/Datalake/'
    elif args.cluster is True:
        archivepath = '/pool/sao/gtremblay/HRC/'
    else:
        archivepath = args.archivepath

    if not os.path.isdir(archivepath):
        sys.exit('Supplied archive Path ({}) not found'.format(archivepath))

    return savepath, archivepath


def inventoryArchive(archivepath, limit=None, verbose=False, sort=False):
    ''' Parse the archive'''

    # Check to make sure the HRC database path is right
    if (sys.version_info > (3, 0)):
        # Python 3 code
        evt1_files = glob.glob(archivepath + '**/*evt1*', recursive=True)
    else:
        # Python 2 code
        # Python <3.5 glob can't walk directories recursively
        import fnmatch
        evt1_files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(
            archivepath) for f in fnmatch.filter(files, '*evt1*')]

    if len(evt1_files) == 0:
        sys.exit(
            'ERROR: No EVT1 files round in supplied archive path ({})'.format(archivepath))

    if limit is None:
        master_list = evt1_files
    else:
        if not isinstance(limit, int):
            sys.exit('ERROR: limit passed to inventoryArchive must be an integer.')
        if verbose is True:
            print('Limiting archive crawl to {} observations'.format(limit))
        master_list = evt1_files[:limit]

    if sort is True:
        # Let's split out by detector to keep things clean.
        hrcI_files = []
        hrcS_files = []

        for evt1_file in master_list:
            # Grab header for the second extension:
            hdr = fits.getheader(evt1_file, 1)
            if hdr['DETNAM'] == 'HRC-I':
                hrcI_files.append(evt1_file)
            elif hdr['DETNAM'] == 'HRC-S':
                hrcS_files.append(evt1_file)

            if verbose is True:
                print("Sorting {} | ObsID {}, {}, {} ksec".format(
                    evt1_file.split('/')[-1], hdr['OBS_ID'], hdr['DETNAM'], round(hdr['EXPOSURE'] / 1000, 2)))

        return evt1_files, hrcI_files, hrcS_files
    elif sort is False:
        print('Archive parsed. {} EVT1 files found.'.format(len(master_list)))
        return master_list


def poolScreen(evt1file, verbose=False, savepath=None, make_reportCard=True, show=False):  # pragma: no cover

    obs = hypercore.HRCevt1(evt1file)

    if verbose is True:
        print("Gathering HyperScreen performance statistics for {} | {}, {} ksec".format(
            obs.obsid, obs.detector, round(obs.exptime/1000.,)))

    if make_reportCard is True:
        reportCard(obs, show=show, savepath=savepath)
        if verbose is True:
            print("Report Card generated for {} | {}, {} ksec".format(
                obs.obsid, obs.detector, round(obs.exptime/1000.,)))

    try:
        results_dict = obs.hyperscreen()
        if make_reportCard is True:
            reportCard(obs, show=False, savepath=savepath)
        return results_dict
    except:
        print("ERROR on {} ({} | {} ksec | {:,} events | {:,} good time events), pressing on".format(
            obs.obsid, obs.detector, round(obs.exptime/1000, 2), obs.numevents, obs.goodtimeevents))


def screenArchive(evt1_file_list, savepath=None, verbose=False, make_reportCard=True, create_pickle=False, picklename=None, show=False):  # pragma: no cover

    p = multiprocessing.Pool()

    # This is how you pass a keyword argument to a pool.Map
    kwargs = {'verbose': verbose,  # be chatty
              # save the products, like report cards and hyperscreen results list-'o-dicts
              'savepath': savepath,
              'make_reportCard': make_reportCard,  # make report cards?
              'show': show}  # show these? *** DEFINITELY a bad idea if you're screening more than 10 evt1 files! ***

    # Passing kwargs to poolScreen requires wrapping with partial()
    hyperscreen_dicts = p.map(partial(poolScreen, **kwargs), evt1_file_list)

    # Now close the Pool
    p.close()
    p.join()

    pickle_set = create_pickle is True and picklename is not None
    pickle_unspecified = create_pickle is True and picklename is None

    if pickle_set:
        with open(savepath + picklename, 'wb') as gherkin:
            pickle.dump(hyperscreen_dicts, gherkin)
            if verbose:
                print("Pickled List of HyperScreen Result Dictonaries. Saved to {}".format(
                    savepath+picklename))
    elif pickle_unspecified:
        raise Exception(
            'create_pickle is True but picklename is None (i.e. unspecified). Please give a pickle name!')

    return hyperscreen_dicts


def main():  # pragma: no cover
    """Console script for hyperscreen."""

    args = getArgs()

    savepath, archivepath = setPaths(args)
    evt1_files = inventoryArchive(
        archivepath, limit=None, verbose=args.verbose, sort=False)

    hyperscreen_dicts = screenArchive(
        evt1_files, savepath=savepath, verbose=args.verbose, make_reportCard=args.reportcard, create_pickle=False, picklename=args.picklename)

    # improvement=[]
    # exptime=[]

    # hyperscreen.styleplots()

    # for observation in evt1_files:
    #     obs=hyperscreen.HRCevt1(observation)
    #     # obs.boomerang(mask=obs.data['Hyperbola test passed'])
    #     try:
    #         results=obs.hyperscreen()
    #         improvement.append(results['Percent improvement'])
    #         exptime.append(obs.exptime / 1000)
    #     except:
    #         print("ERROR on {} ({}, {} ksec), pressing on".format(
    #             obs.obsid, obs.detector, round(obs.exptime/1000), 2))
    #         continue

    # fig, ax=plt.subplots()
    # ax.plot(exptime, improvement, linewidth=0, marker='.')
    # ax.set_xlabel("Exposure Time (ksec)")
    # ax.set_ylabel("Percent Improvement over Legacy Test")

    # plt.show()

    # reportCard(obs, savepath=savepath, verbose=verbose, rasterized=True)

    # hyperscreen.styleplots()

    # for evt1_file in evt1_files:
    #     reportCard(evt1_file, savepath=savepath)

    # print(hyperscreen_dicts)
    # for evt1_file in evt1_files:
    #     obs = hyperscreen.HRCevt1(evt1_file)
    #     tapscreen_results_dict = obs.hyperscreen()
    #     # obs.image(show=False, detcoords=True,
    #     #           savepath="/Users/grant/Desktop/image_test/{}.pdf".format(obs.obsid), create_subplot=False)
if __name__ == "__main__":

    start_time = time.time()
    main()  # pragma: no cover
    runtime = round((time.time() - start_time) / 60, 3)
    sys.exit(print("Finished in {} minutes".format(runtime)))
