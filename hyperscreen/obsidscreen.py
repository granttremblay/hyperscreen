#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""

from __future__ import division
from __future__ import print_function


import os
import sys
from shutil import copyfile
import time
import glob
import argparse

from astropy.io import fits

from hyperscreen import hypercore
from hyperscreen import evtscreen


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('obsid_directory', nargs='?', type=os.path.abspath)

    return parser.parse_args(argv)


def main():

    args = getArgs()

    evt1_file_path = glob.glob(
        args.obsid_directory + '/secondary/*evt1*', recursive=True)

    if len(evt1_file_path) == 0:
        raise Exception(
            "ERROR - can't find an EVT1 file in the /secondary/ directory - make sure you pass the top-level ObsID directory (e.g. ~/Destkop/1501/)")

    if len(evt1_file_path) > 1:
        raise Exception("ERROR: more than one evt1 file found!")

    evt1_file = evt1_file_path[0]
    savepath = os.path.dirname(evt1_file)

    hdr = fits.getheader(evt1_file, 1)
    print(" -------    Welcome to HyperScreen  --------")
    print("ObsID {} | {} | {} | {} ksec".format(
        hdr['OBS_ID'], hdr['DETNAM'], hdr['OBJECT'], round(hdr['EXPOSURE']/1000, 2)))

    hdr = fits.getheader(evt1_file_path[0])
    evtscreen.screenHRCevt1(evt1_file_path[0], verbose=True,
                            comparison_products=True, savepath=savepath)


if __name__ == "__main__":

    start_time = time.time()
    main()  # pragma: no cover
    runtime = round((time.time() - start_time) / 60, 3)
    sys.exit(print("Finished in {} minutes".format(runtime)))
