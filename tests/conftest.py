#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Configure PyTest fixtures for the test modules
"""

import os
import pytest

from hyperscreen import hyperscreen


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