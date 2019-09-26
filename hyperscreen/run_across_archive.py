#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Console script for hyperscreen."""
import multiprocessing
import hyperscreen
import os
import sys
import time
import glob
import argparse


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

    args = parser.parse_args()

    archive_path = args.archivepath
    if not os.path.isdir(archive_path):
        sys.exit("Supplied archive Path ({}) not found".format(archive_path))

    # Check to make sure the HRC database path is right
    evt1_files = glob.glob(archive_path + '**/*evt1*', recursive=True)
    if len(evt1_files) == 0:
        sys.exit(
            "No EVT1 files round in supplied archive path ({})".format(archive_path))

    # p = multiprocessing.Pool()
    # p.map(clean, evt1_files[:4])
    # p.close()
    # p.join()

    for evt1_file in evt1_files[:1]:
        obs = hyperscreen.HRCevt1(evt1_file)
        tapscreen_results_dict = obs.hyperscreen()
        obs.boomerang()
    # for obs in evt1_files[4]:
    #     results = clean(obs)
    #     print(results)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
