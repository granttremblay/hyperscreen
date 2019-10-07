#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
py.test unit tests for hyperscreen and related packages. 
Currently, these are pretty dumb. 
"""

import sys
import os

import pytest

from hyperscreen import hyperscreen


@pytest.fixture()
def hrcI_evt1():
    hrcI_file = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcI_evt1_testfile.fits.gz')
    hrcI_evt1 = hyperscreen.HRCevt1(hrcI_file)
    return hrcI_evt1

@pytest.fixture()
def hrcS_evt1():
    hrcS_file = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/data/hrcS_evt1_testfile.fits.gz')
    hrcS_evt1 = hyperscreen.HRCevt1(hrcS_file)
    return hrcS_evt1

def test_hrcI_evt1(hrcI_evt1):
    assert len(hrcI_evt1.data['fp_u']) == len(hrcI_evt1.data['time'])

def test_hrcS_evt1(hrcS_evt1):
    assert len(hrcS_evt1.data['fp_u']) == len(hrcS_evt1.data['time'])

def test_gti_masks(hrcI_evt1, hrcS_evt1):
    match_I = hrcI_evt1.badtimeevents == hrcI_evt1.numevents - hrcI_evt1.goodtimeevents
    match_S = hrcS_evt1.badtimeevents == hrcS_evt1.numevents - hrcS_evt1.goodtimeevents
    assert match_I
    assert match_S

# def test_legacy_hyperbola_status_bits(hrcI_evt1, hrcS_evt1)