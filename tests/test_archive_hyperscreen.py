#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
pytest unit tests for hyperscreen and related packages. 
"""

from __future__ import division
from __future__ import print_function

import sys
import os
from contextlib import contextmanager

import astropy
from astropy.io import fits

import pandas as pd

import matplotlib.pyplot as plt

import pytest

from hyperscreen import archive_hyperscreen
    
'''
This test module uses pytest Fixtures defined in conftest.py
'''

# I STILL GET AN 'HYPERSCREEN' HAS NO ATTRIBUTE 'HRCEVT1' WHEN CALLING THIS WITH PYTEST. 
# I THINK IT HAS SOMETHING TO DO WITH MY __INIT__.PY IN THE HYPERSCREEN MODULE. 
# I STILL DON'T REALLY UNDERSTAND IMPORT STRUCTURE

# def test_main():
#     results_dict = archive_hyperscreen.screenArchive([os.path.dirname(os.path.abspath(__file__))+'/data/hrcI_evt1_testfile.fits.gz'], create_pickle=False)
#     assert len(results_dict) > 1
    
