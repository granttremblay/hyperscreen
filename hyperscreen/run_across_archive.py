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


# def clean(evt1_file):
#     obs = hyperscreen.HRCevt1(evt1_file)
#     print("Doing {}, {}, {} events".format(
#         obs.obsid, obs.detector, obs.numevents))
#     tapscreen_results_dict = obs.tapscreen()
#     hyperscreen.image(x[survival_mask], y[survival_mask], title='{} | {} | {}'.format(
#         obs.obsid, obs.target, obs.numevents), show=False, savepath="/Users/grant/Desktop/hyperplots/{}.pdf".format(obs.obsid))
#     #print(tapscreen_results_dict["Percent improvement"])
#     return tapscreen_results_dict


def main():
    """Console script for hyperscreen."""
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--archivepath', help='Absolute PATH to Archive of EVT1 Files',
                        default='/Users/grant/Science/HRC_Database/EVT1_Files/')

    parser.add_argument('-s', '--savepath', help='Absolute PATH to location in which to save all outputs from this script, including .pdf files of plots. If not specified, this location will default to your Desktop.',
                        default=os.path.join(os.environ['HOME'], 'Desktop/'))

    parser.add_argument('-t', '--testdata', action='store_true',
                        help='Use the supplied test data as the input archive path')

    parser.add_argument('-w', '--windowstest', action='store_true',
                        help='Point to my Windows database')

    args = parser.parse_args()

    savepath = args.savepath

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

    # p = multiprocessing.Pool()
    # p.map(clean, evt1_files[:4])
    # p.close()
    # p.join()

    # for evt1_file in evt1_files:
    #     obs = hyperscreen.HRCevt1(evt1_file)
    #     tapscreen_results_dict = obs.hyperscreen()
    #     # obs.image(show=False, detcoords=True,
    #     #           savepath="/Users/grant/Desktop/image_test/{}.pdf".format(obs.obsid), create_subplot=False)

    obs = hyperscreen.HRCevt1(evt1_files[1])

    hyperscreen.styleplots()

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
              create_subplot=True, title="Test2", cmap='magma')

    plt.show()
    fig.savefig(reportCard_savepath)
    print("Created {}".format(reportCard_savepath))

    # for obs in evt1_files[4]:
    #     results = clean(obs)
    #     print(results)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
