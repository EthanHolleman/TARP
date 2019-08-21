import unittest
from Glycine_Remap.pairing import linker
from Glycine_Remap.pairing import get_path


class Test_Pairing(unittest.TestCase):
    def test_linker(self):
        GOLD_FILE = 'Gmr30.fna'
        ANS = tuple(['consensus_Gmr30SOLO.fna', 'consensus_Gmr30INTACT.fna'])
        result = linker(GOLD_FILE)
        self.assertIs(ANS, result)

    def test_get_path(self):
        GOLD_DIR = 'test/ing/this/dir'
        ANS = 'test/ing/this'
        result = get_path(GOLD_DIR)
        self.assertIs(ANS, result)
