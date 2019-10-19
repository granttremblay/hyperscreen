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

from astropy.io import fits


def getArgs(argv=None):
    parser = argparse.ArgumentParser(
        description='Test HyperScreen across an archive of HRC Observations')

    parser.add_argument('inputFile', nargs='?', type=argparse.FileType('rt'))

    return parser.parse_args(argv)


def screenHRCevt1(args):

    
    inputFile = args.inputFile

    print(inputFile)
    obs = hypercore.HRCevt1(inputFile)

    hyperscreen_results_dict = obs.hyperscreen()
    mask = hyperscreen_results['All Survivals (boolean mask)']

    with fits.open(inputFile) as hdul:
        hdu = hdul[1]
        hdul.info()
        # hdu.data = hdu.data[mask]
        # hdu.writeto('test_clean.fits')

def main():

    args = getArgs()

    screenHRCevt1(args)


if __name__ == "__main__":

    start_time = time.time()
    main()  # pragma: no cover
    runtime = round((time.time() - start_time) / 60, 3)
    sys.exit(print("Finished in {} minutes".format(runtime)))
