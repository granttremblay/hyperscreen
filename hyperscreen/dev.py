#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import glob
import argparse


import hyperscreen

archive_path = '/Users/grant/Science/HRC_Database/EVT1_Files/'
evt1_files = glob.glob(archive_path + '**/*evt1*', recursive=True)

for observation in evt1_files:
    evt1 = hyperscreen.HRCevt1(observation)
    instrument = evt1.detector
    size = evt1.numevents
    obsid = evt1.obsid
    print('{} is {} with {} events'.format(obsid, instrument, size))
