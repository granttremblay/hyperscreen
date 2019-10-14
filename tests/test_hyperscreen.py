#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PyTest unit tests for hyperscreen and related packages. 
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

from hyperscreen import hyperscreen

# myPath = os.path.dirname(os.path.abspath(__file__))
# sys.path.insert(0, myPath + '/../')

# @pytest.fixture(scope="module")
# def fits_object():
#     print("Reading raw HRC-I EVT1 file with astropy.io.fits")
#     rawfile = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcI_evt1_testfile.fits.gz')
#     fits_object = fits.open(rawfile)
#     yield fits_object
#     fits_object.close()
#     print("Closing the raw fits file")
    

@pytest.fixture(scope="module")
def hrcI_evt1():
    print("Loading test HRC-I EVT1 File as a pytest Fixture")
    hrcI_file = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcI_evt1_testfile.fits.gz')
    hrcI_evt1 = hyperscreen.HRCevt1(hrcI_file)
    return hrcI_evt1

@pytest.fixture(scope="module")
def hrcS_evt1():
    print("Loading test HRC-S EVT1 File as a pytest Fixture")
    hrcS_file = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcS_evt1_testfile.fits.gz')
    hrcS_evt1 = hyperscreen.HRCevt1(hrcS_file)
    return hrcS_evt1

@contextmanager
def assert_plot_figures_added():
    num_figures_before = plt.gcf().number
    yield
    num_figures_after = plt.gcf().number
    assert num_figures_before < num_figures_after


# def test_hrcI_evt1(hrcI_evt1):
#     assert len(hrcI_evt1.data['fp_u']) == len(hrcI_evt1.data['time'])

# def test_hrcS_evt1(hrcS_evt1):
#     assert len(hrcS_evt1.data['fp_u']) == len(hrcS_evt1.data['time'])

# def test_gti_masks(hrcI_evt1, hrcS_evt1):
#     match_I = hrcI_evt1.badtimeevents == hrcI_evt1.numevents - hrcI_evt1.goodtimeevents
#     match_S = hrcS_evt1.badtimeevents == hrcS_evt1.numevents - hrcS_evt1.goodtimeevents
#     assert match_I
#     assert match_S

def test_HRCevt1(hrcI_evt1, hrcS_evt1):
    assert len(hrcI_evt1.data['fp_u']) == len(hrcI_evt1.data['time'])
    assert len(hrcS_evt1.data['fp_u']) == len(hrcS_evt1.data['time'])

    match_I = hrcI_evt1.badtimeevents == hrcI_evt1.numevents - hrcI_evt1.goodtimeevents
    match_S = hrcS_evt1.badtimeevents == hrcS_evt1.numevents - hrcS_evt1.goodtimeevents
    assert match_I
    assert match_S

def test_astropy_return():
    hrcI_file = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcI_evt1_testfile.fits.gz')
    hrcI_evt1_df = hyperscreen.HRCevt1(hrcI_file)
    hrcI_evt1_table = hyperscreen.HRCevt1(hrcI_file, as_astropy_table=True)
    assert isinstance(hrcI_evt1_df.data, pd.DataFrame)
    assert type(hrcI_evt1_table.data) is astropy.table.table.Table



def test_hyperscreen(hrcI_evt1, hrcS_evt1):
    hyperscreen_hrcI = hrcI_evt1.hyperscreen()
    print(hyperscreen_hrcI)
    assert isinstance(hyperscreen_hrcI, dict)
    # This will give a warning. You can have pytest ignore with pytest --disable-warnings

def test_boomerang(hrcI_evt1, hrcS_evt1):
    with assert_plot_figures_added():
        hrcI_evt1.boomerang(mask=hrcI_evt1.data['Hyperbola test failed'], show=False)
        hrcS_evt1.boomerang(show=False)



    
