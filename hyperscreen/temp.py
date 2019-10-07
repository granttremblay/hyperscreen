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

