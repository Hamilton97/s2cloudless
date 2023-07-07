import unittest
import ee

ee.Initialize()
from s2cloudless import (S2CloudlessBuilder, 
                         S2CloudCollection,
                         S2CloudProbabilityCollection,
                         S2SRCollection)


class TestS2Cloudless(unittest.TestCase):
    def setUp(self) -> None:
        self.aoi = None
        self.date_range = None
        return super().setUp()
    
    def test_S2CloudlessBuilder(self):
        builder = S2CloudlessBuilder()
        self.assertIsInstance(builder, S2CloudlessBuilder)
    
    def test_S2CloudProbabilityCollection(self):
        col = S2CloudProbabilityCollection()
        self.assertIsInstance(col, S2CloudProbabilityCollection)
    
    def test_S2SRCollection(self):
        col = S2SRCollection()
        self.assertIsInstance(col, S2SRCollection)
    
    def test_S2CloudCollection(self):
        cld_prb = S2CloudProbabilityCollection()
        sr_sr = S2SRCollection()
        col = S2CloudCollection(cld_prb, sr_sr)
        self.assertIsInstance(col, S2CloudCollection)