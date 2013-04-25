#!/usr/bin/env python
# testing ProBin script
from nose.tools import assert_almost_equal, assert_equal
import probin
import os
import sys
import fileinput

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestProBin(object):
    def setUp(self):
        pass
    def tearDown(self):
        pass
    def test_main(self):
        assert_equal(True,True)
    def test_none(self):
        assert_equal(None,None)
