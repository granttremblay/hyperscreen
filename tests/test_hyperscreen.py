#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `hyperscreen` package."""



import sys
import os
# Set the path explicitly #
sys.path.insert(0, os.path.abspath(__file__+"/../.."))
from hyperscreen import hyperscreen

import pytest


@pytest.fixture()
def loadHRCI():
    hrcI = hyperscreen.HRCevt1('data/hrcI_evt1_testfile.fits.gz')
    return hrcI

def test_structure(loadHRCI):
    assert hrcI.numevents > 0

