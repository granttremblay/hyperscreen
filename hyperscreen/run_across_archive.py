#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function

import os
import sys
import time
import glob
import argparse


import multiprocessing
import pickle
from functools import partial


import hyperscreen as hyperscreen


from astropy.io import fits

import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


def reportCard(evt1_object, savepath, show=True, save=True, rasterized=True, dpi=150, verbose=False):

    obs = evt1_object

    if verbose is True:
        print("Doing {}, {}".format(obs.obsid, obs.detector))

    reportCard_savepath = os.path.join(
        savepath, '{}_{}_{}_hyperReport.pdf'.format(obs.obsid, obs.target, obs.detector))
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
        fig.savefig(reportCard_savepath, rasterized=rasterized, dpi=dpi)

    if verbose is True:
        print("Created {}".format(reportCard_savepath))

    if show is True:
        plt.show()

    # for obs in evt1_files[4]:
    #     results = clean(obs)
    #     print(results)


def multiprocess_clean(evt1_file):
    try:
        obs = hyperscreen.HRCevt1(evt1_file)
        hyperscreen_dict = obs.hyperscreen()
        print("Finished ObsID {}, {}, {}, {} ksec".format(
            obs.obsid, obs.target, obs.detector, obs.exptime))
        return hyperscreen_dict
    except:
        print("PROBLEM with {}, skipping for now.".format(
            evt1_file.split('/')[-1]))
        return None


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('-a', '--archivepath', help='Absolute PATH to Archive of EVT1 Files',
                        default='/Users/grant/Science/HRC_Database/EVT1_Files/')

    parser.add_argument('-h', '--hydra', action='store_true',
                        help='Point to the HRC Database stored on the Smithsonian Hydra Cluster')   

    parser.add_argument('-p', '--picklename', help='Name of the Pickle you would like to create',
                        default='hyperscreen_master_pickle.pkl')

    parser.add_argument('-s', '--savepath', help='Absolute PATH to location in which to save all outputs from this script, including .pdf files of plots. If not specified, this location will default to your Desktop.',
                        default=os.path.join(os.environ['HOME'], 'Desktop/HyperScreen_Results/'))

    parser.add_argument('-t', '--testdata', action='store_true',
                        help='Use the supplied test data as the input archive path')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Make HyperScreen chatty on stdout.')

    parser.add_argument('-w', '--windowstest', action='store_true',
                        help='Point to my Windows database')

    return parser.parse_args(argv)


def setPaths(args, verbose=False):
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
        archivepath = '../tests/data/'
    elif args.windowstest is True:
        archivepath = '/mnt/c/Users/grant/HRCOps/Datalake/'
    elif args.hydra is True:
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


def poolClean(evt1file, verbose=True):

    obs = hyperscreen.HRCevt1(evt1file)

    if verbose is True:
        print("Gathering HyperScreen performance statistics for {} | {}, {} ksec".format(
            obs.obsid, obs.detector, round(obs.exptime/1000.,)))

    try:
        results_dict = obs.hyperscreen()
        return results_dict
    except:
        print("ERROR on {} ({}, {} ksec), pressing on".format(
            obs.obsid, obs.detector, round(obs.exptime/1000), 2))


def screen_and_pickle_archive(evt1_file_list, savepath, picklename):

    p = multiprocessing.Pool()
    hyperscreen_dicts = p.map(poolClean, evt1_file_list)
    p.close()
    p.join()

    with open(savepath + picklename, 'wb') as gherkin:
        pickle.dump(hyperscreen_dicts, gherkin)
        print("Pickled List of HyperScreen Result Dictonaries. Saved to {}".format(
            savepath+picklename))


def main():
    """Console script for hyperscreen."""

    args = getArgs()

    verbose = args.verbose

    savepath, archivepath = setPaths(args)
    evt1_files = inventoryArchive(
        archivepath, limit=None, verbose=verbose, sort=False)

    screen_and_pickle_archive(
        evt1_files, savepath=savepath, picklename=args.picklename)

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
    sys.exit(main())
