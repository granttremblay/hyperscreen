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

from hyperscreen import archivescreen

'''
This test module uses pytest Fixtures defined in conftest.py
'''

# I STILL GET AN 'HYPERSCREEN' HAS NO ATTRIBUTE 'HRCEVT1' WHEN CALLING THIS WITH PYTEST.
# I THINK IT HAS SOMETHING TO DO WITH MY __INIT__.PY IN THE HYPERSCREEN MODULE.
# I STILL DON'T REALLY UNDERSTAND IMPORT STRUCTURE


def test_inventoryArchive():
    archivepath = os.path.dirname(os.path.abspath(__file__))+'/data/'
    master_list = archivescreen.inventoryArchive(archivepath, sort=True, verbose=True)
    print(master_list)
    assert len(master_list) > 1


def test_getArgs():
    parser = archivescreen.getArgs(['-v', '-c', '--savepath=/hello/', '--archivepath=/hi/there/'])
    assert parser.verbose is True
    assert parser.cluster is True
    assert parser.savepath == '/hello/'
    assert parser.archivepath == '/hi/there/'

def test_reportCard(hrcI_evt1):
    from .test_hypercore import assert_plot_figures_added
    with assert_plot_figures_added():
        archivescreen.reportCard(hrcI_evt1, savepath='blah', show=False, save=False)

