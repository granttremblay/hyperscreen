
"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function


from hyperscreen import evtscreen
from hyperscreen import hypercore
import gc
import os
import sys
import time
import glob
import argparse


import multiprocessing

import json

import pickle
from functools import partial

from astropy.io import fits

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('-r', '--results_dir', help='Absolute PATH to directory of archivescreen results.',
                        default='/Users/grant/Desktop/hyperscreen_results/')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Make HyperScreen chatty on stdout.')

    return parser.parse_args(argv)


def inventoryJSONs(results_dir, verbose=False):
        # Check to make sure the HRC database path is right
    if (sys.version_info > (3, 0)):
        # Python 3 code
        json_files = glob.glob(results_dir + '*.json', recursive=True)
    else:
        # Python 2 code
        # Python <3.5 glob can't walk directories recursively
        import fnmatch
        json_files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(
            results_dir) for f in fnmatch.filter(files, '*.json')]

    if len(json_files) == 0:
        sys.exit(
            'ERROR: No JSON HyperScreen Result files round in supplied archive path ({})'.format(results_dir))

    if verbose is True:
        print('{} JSON HyperScreen Result Files found'.format(len(json_files)))
    return json_files


def parseJSONs(json_files, verbose=False):

    exptimes = []
    improvements = []
    legacy_percent = []
    hyperscreen_percent = []

    for json_file in json_files:
        if verbose is True:
            print("Parsing {}".format(json_file.split('/')[-1]))
        with open(json_file) as json_data:
            data = json.load(json_data)
            exptimes.append(data['Exposure Time'])
            improvements.append(data['Percent improvement'])
            legacy_percent.append(data['Percent rejected by Hyperbola'])
            hyperscreen_percent.append(data['Percent rejected by Tapscreen'])

    trends_dict = {'Exposure Times': exptimes,
                   'Improvement': improvements,
                   'Legacy Hyperbola Test Percent': legacy_percent,
                   'Hyperscreen Percent': hypercore}

    return trends_dict


def main():  # pragma: no cover
    """[summary]
    """

    hypercore.styleplots()

    args = getArgs()

    results_dir = args.results_dir

    if os.path.exists(results_dir):
        pass
    else:
        raise Exception(
            'Supplied results directory ({}) does not exist.'.format(results_dir))

    json_files = inventoryJSONs(results_dir, verbose=args.verbose)
    trends_dict = parseJSONs(json_files, verbose=args.verbose)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(trends_dict['Exposure Times'],
            trends_dict['Improvement'], linewidth=0, marker='.')

    plt.show()


if __name__ == "__main__":

    start_time = time.time()
    main()  # pragma: no cover
    runtime = round((time.time() - start_time) / 60, 3)
    sys.exit(print("Finished in {} minutes".format(runtime)))
