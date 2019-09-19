#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""
import os
import sys
import time
import glob
import argparse

# from astropy.io import fits
# from astropy.table import Table

# import pandas as pd
# import numpy as np
# np.seterr(divide='ignore')

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# # from matplotlib.backends.backend_pdf import PdfPages

import warnings
warnings.filterwarnings("ignore",category =RuntimeWarning)

import hyperscreen as hc
import multiprocessing


def clean(evt1_file):
    obs = hc.HRCevt1(evt1_file, as_dataframe=True)
    print("Doing {}, {}, {} events".format(obs.obsid, obs.detector, obs.numevents))  
    tapscreen_results_dict = obs.tapscreen()
    #hc.image(x[survival_mask],y[survival_mask], title='{} | {} | {}'.format(obs.obsid, obs.target, obs.numevents), show=False, savepath="/Users/grant/Desktop/hyperplots/{}.pdf".format(obs.obsid))
    #print(tapscreen_results_dict["Percent improvement"])


def main():
    """Console script for hyperscreen."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', help='Test Argument')
    args = parser.parse_args()

    print(args.test)

    archive_path='/Users/grant/Science/HRC_Database/EVT1_Files/'
    if not os.path.isdir(archive_path):
        sys.exit("Supplied archive Path ({}) not found".format(archive_path))

    # Check to make sure the HRC database path is right
    evt1_files = glob.glob(archive_path + '**/*evt1*', recursive=True)
    if len(evt1_files) == 0:
        sys.exit("No EVT1 files round in supplied archive path ({})".format(archive_path))

    # p = multiprocessing.Pool()
    # p.map(hyperclean, evt1_files[:4])
    # p.close()
    # p.join()

    for obs in evt1_files[:1]:
        clean(obs)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
