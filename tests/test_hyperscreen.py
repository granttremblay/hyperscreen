#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `hyperscreen` package."""

import pytest
import os
from hyperscreen import hyperscreen


HRC_I_TESTDATA = os.path.join(os.path.dirname(
    __file__), 'data/hrcI_evt1_testfile.fits.gz')
HRC_S_TESTDATA = os.path.join(os.path.dirname(
    __file__), 'data/hrcS_evt1_testfile.fits.gz')

print(HRC_I_TESTDATA)


# class MyTest(unittest.TestCase)

#    def setUp(self):
#         self.testfile = open(TESTDATA_FILENAME)
#         self.testdata = self.testfile.read()

#     def tearDown(self):
#         self.testfile.close()

#     def test_something(self):
#         ....
