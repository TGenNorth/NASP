# coding: utf-8

__author__ = 'jtravis'

import unittest
from nasp.nasp_objects import GenomeCollection

class GenomeCollectionTestCase(unittest.TestCase):
    def setUp(self):
        self.genome = GenomeCollection()