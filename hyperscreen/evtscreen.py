#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function

from hyperscreen import hypercore
import os
import sys
from shutil import copyfile
import time
import glob
import argparse

from astropy.io import fits


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('input_fits_file', nargs='?')

    parser.add_argument('-c', '--comparison_products', action='store_true',
                        help='Make additional HyperScreen result images (rejected events & difference map)')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Make HyperScreen chatty on stdout.')

    # parser.add_argument('-b', '--backup_dir', help='Absolute PATH to backup of EVT1 Files',
    #                     default=None)

    return parser.parse_args(argv)


def screenHRCevt1(input_fits_file, hyperscreen_results_dict=None, comparison_products=True, savepath=None, verbose=True, backup=True):

    obsid = fits.getheader(input_fits_file)['OBS_ID']

    # Get the root string to use as our naming convention
    file_name = input_fits_file.split('/')[-1]  # Split off the path
    # Split off the .fits (or .fits.gz)
    file_root = file_name.split('.fits')[0].split('_')[0]
    file_path = os.path.realpath(os.path.dirname(input_fits_file))
    print(file_path)
    if savepath is None:
        savepath = file_path
    backup_dir = os.path.join(
        savepath, '{}_hyperscreen_report'.format(obsid))

    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
        if verbose is True:
            print("Made HyperScreen Products Directory {}".format(backup_dir))

    # Check if the passed input file is a string ending with .fits or .fits.gz
    if isinstance(input_fits_file, str):
        filetype_match = input_fits_file.split(
            '.')[-1] == 'fits' or input_fits_file.split('.')[-1] == 'gz'

        if filetype_match is False:
            raise Exception(
                'ERROR: Input given ({}) is not recognized as a .fits[.gz] or HRCevt1 object.'.format(input_fits_file))

    if hyperscreen_results_dict is None:
        # Then you need to make it!
        if verbose is True:
            print("Applying HyperScreen algorithm to {}".format(file_name))
        obs = hypercore.HRCevt1(input_fits_file)
        hyperscreen_results = obs.hyperscreen()
    else:
        hyperscreen_results = hyperscreen_results_dict

    # if backup is True:
    #     if not os.path.exists(backup_dir):
    #         os.makedirs(backup_dir)
    #         if verbose is True:
    #             print("Made backup directory {}".format(backup_dir))

    #     if verbose is True:
    #         print("Backing up (copying) input fits file {} to {}".format(
    #             input_fits_file.split('/')[-1], backup_dir))
    #     copyfile(input_fits_file, os.path.join(
    #         file_path, '{}_original_event_list.fits'.format(obsid)))

    # hyperscreen_fits_file = file_path + '/hyperscreen_' + file_name
    hyperscreen_fits_file = os.path.join(file_path, 'hyperscreen_' + file_name)

    rejected_events_file = os.path.join(
        backup_dir, '{}_hyperscreen_rejected_events.fits'.format(obsid))
    # difference_map_file = file_path + '/hyperscreen_DIFFERENCE_MAP_' + file_name

    survival_mask = hyperscreen_results['All Survivals (boolean mask)']
    failure_mask = hyperscreen_results['All Failures (boolean mask)']

    with fits.open(input_fits_file) as hdul:
        original_data = hdul[1].data
        hdul[1].data = hdul[1].data[survival_mask]
        if verbose is True:
            print("Masking data with HyperScreen Results")
        hdul.writeto(hyperscreen_fits_file, overwrite=True)

        if verbose is True:
            print(
                "Wrote new HyperScreen-filtered evt1 file to {}".format(hyperscreen_fits_file))

        if comparison_products is True:
            hdul[1].data = original_data[failure_mask]
            hdul.writeto(rejected_events_file, overwrite=True)
            if verbose is True:
                print("Wrote Rejected Events Map {}".format(rejected_events_file))

    original_evt1_file_path = os.path.join(
        file_path, '{}_original_event_list.fits'.format(obsid))
    os.rename(input_fits_file, original_evt1_file_path)
    print("Backed up original evt1 file to {}".format(original_evt1_file_path))


def main():

    args = getArgs()

    screenHRCevt1(args.input_fits_file, verbose=True,
                  comparison_products=args.comparison_products)


if __name__ == "__main__":

    start_time = time.time()
    main()  # pragma: no cover
    runtime = round((time.time() - start_time) / 60, 3)
    sys.exit(print("Finished in {} minutes".format(runtime)))
