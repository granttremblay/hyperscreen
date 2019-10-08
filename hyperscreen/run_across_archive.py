#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function

import multiprocessing
import hyperscreen as hyperscreen
import os
import sys
import time
import glob
import argparse

import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


def reportCard(hrcEVT1_fits_file, savepath):

    obs = hyperscreen.HRCevt1(hrcEVT1_fits_file)

    print("Doing {}, {}".format(obs.obsid, obs.detector))

    reportCard_savepath = os.path.join(
        savepath, '{}_{}_{}_hyperReport.pdf'.format(obs.obsid, obs.target, obs.detector))
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    obs.boomerang(ax=axes[0, 0], create_subplot=True,
                  show=False, title='Test1', cmap='magma')
    obs.boomerang(ax=axes[0, 1], create_subplot=True,
                  show=False, title='Test2', cmap='inferno')

    obs.image(ax=axes[1, 0], detcoords=True, show=False,
              create_subplot=True, title="Test1")
    obs.image(ax=axes[1, 1], detcoords=True, show=False,
              create_subplot=True, title="Test2")

    fig.savefig(reportCard_savepath)

    test = obs.hyperscreen()
    print("Created {}".format(reportCard_savepath))

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


def main():
    """Console script for hyperscreen."""
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--archivepath', help='Absolute PATH to Archive of EVT1 Files',
                        default='/Users/grant/Science/HRC_Database/EVT1_Files/')

    parser.add_argument('-s', '--savepath', help='Absolute PATH to location in which to save all outputs from this script, including .pdf files of plots. If not specified, this location will default to your Desktop.',
                        default=os.path.join(os.environ['HOME'], 'Desktop/HyperScreen_Results/'))

    parser.add_argument('-t', '--testdata', action='store_true',
                        help='Use the supplied test data as the input archive path')

    parser.add_argument('-w', '--windowstest', action='store_true',
                        help='Point to my Windows database')

    args = parser.parse_args()

    savepath = args.savepath
    if not os.path.exists(savepath):
        print("Creating Directory {} in which HyperScreen results will be saved.".format(
            savepath))
        os.makedirs(savepath)

    print("Savepath is {}".format(savepath))

    if args.testdata is True:
        archive_path = '../tests/data/'
    elif args.windowstest is True:
        archive_path = '/mnt/c/Users/grant/HRCOps/Datalake/'
    else:
        archive_path = args.archivepath

    if not os.path.isdir(archive_path):
        sys.exit('Supplied archive Path ({}) not found'.format(archive_path))

    # Check to make sure the HRC database path is right
    if (sys.version_info > (3, 0)):
        # Python 3 code
        evt1_files = glob.glob(archive_path + '**/*evt1*', recursive=True)
    else:
        # Python 2 code
        # Python <3.5 glob can't walk directories recursively
        import fnmatch
        evt1_files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(
            archive_path) for f in fnmatch.filter(files, '*evt1*')]

    if len(evt1_files) == 0:
        sys.exit(
            'No EVT1 files round in supplied archive path ({})'.format(archive_path))

    # hyperscreen.styleplots()

    # for evt1_file in evt1_files:
    #     reportCard(evt1_file, savepath=savepath)

    p = multiprocessing.Pool()
    hyperscreen_dicts = p.map(
        multiprocess_clean, evt1_files[:100])
    p.close()
    p.join()

    print(hyperscreen_dicts)

    # for evt1_file in evt1_files:
    #     obs = hyperscreen.HRCevt1(evt1_file)
    #     tapscreen_results_dict = obs.hyperscreen()
    #     # obs.image(show=False, detcoords=True,
    #     #           savepath="/Users/grant/Desktop/image_test/{}.pdf".format(obs.obsid), create_subplot=False)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
