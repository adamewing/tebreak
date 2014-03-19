#!/usr/bin/env python

import rcseq
import unittest

class TestRCseq(unittest.TestCase):
    def test_rc(self):
        ''' reverse complement function '''
        self.assertEqual(rcseq.rc('ATCGatcg'), 'cgatCGAT')
