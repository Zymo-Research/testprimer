from __future__ import absolute_import

import os.path
import types
import StringIO
import unittest

from Bio import SeqIO
from Bio.Seq import Seq

from testprimer.pcr import Template, Primer, PrimerPool, PCRArray, simple_match


class TestTemplate(unittest.TestCase):
    
    def test_attr(self):
        dataset = """\
>GAXI01000525.151.1950 Eukaryota;Opisthokonta;Holozoa;Eumetazoa;Bilateria
.....AACGGUUU---UUUAA-GGGUGAAA.......
"""
        handle = StringIO.StringIO(dataset)
        self.seqrecord = SeqIO.parse(handle, 'fasta').next()
        template = Template(self.seqrecord)

        self.assertTrue(template.seq,
                        ".....AACGGUUU---UUUAA-GGGUGAAA.......")
        self.assertTrue(template.description,
                        ">GAXI01000525.151.1950 Eukaryota;Opisthokonta;Holozoa;Eumetazoa;Bilateria ")
        self.assertTrue(template.id,
                        ">GAXI01000525.151.1950")
        self.assertTrue(template.id,
                        "Eukaryota;Opisthokonta;Holozoa;Eumetazoa;Bilateria")


class TestPrimer(unittest.TestCase):

    def test_class(self):
        self.assertTrue(issubclass(Primer, Seq))


class TestPrimerPool(unittest.TestCase):

    def test_attr(self):
        pp = PrimerPool("515f", 13, 14, [Primer("AA"), Primer("GG")])
        self.assertEqual(pp.name, "515f")
        self.assertEqual(pp.start, 13)
        self.assertEqual(pp.end, 14) 
        self.assertEqual(pp.primers, [Primer("AA"), Primer("GG")])

    def test_input_types(self):
        with self.assertRaises(TypeError):
            pp = PrimerPool("515f", 13, 14, Primer("AA"))
        with self.assertRaises(TypeError):
            pp = PrimerPool("515f", 13, 14, ["AA", "GG"])    
        with self.assertRaises(TypeError):
            pp = PrimerPool("515f", "13", 14,
                            [Primer("AA"), Primer("GG")])
        with self.assertRaises(TypeError):
            pp = PrimerPool("515f", 13, "14",
                            [Primer("AA"), Primer("GG")])


class TestPCRArray(unittest.TestCase):
    
    def setUp(self):
        self.fw_primer_pool_path = os.path.join(
            os.path.dirname(__file__),
            'data',
            'ForwardPrimerPool.txt'
        )
        self.rv_primer_pool_path = os.path.join(
            os.path.dirname(__file__),
            'data',
            'ReversePrimerPool.txt'
        )
        self.template_path = os.path.join(
            os.path.dirname(__file__),
            'data',
            'SILVA_test.fasta'
        )
        self.pcrarray = PCRArray(
            self.template_path,
            self.fw_primer_pool_path, 
            self.rv_primer_pool_path
        )

    def test_parse_primer_pool(self):
        pp = PCRArray.parse_primer_pool(self.fw_primer_pool_path)
        self.assertTrue(isinstance(pp, PrimerPool))
        self.assertEqual(pp.name, '515f')
        self.assertEqual(pp.start, 11895)
        self.assertEqual(pp.end, 13861)
        self.assertEqual(pp.primers, [Primer('GTGCCAGCAGTCGCGGTAA'),
                                      Primer('GTGCCAGCAGGAGCGGTAA')])

    def test_fw_primer_pool(self):
        self.assertTrue(
            isinstance(self.pcrarray.fw_primer_pool, PrimerPool)
        )
        self.assertEqual(len(self.pcrarray.fw_primer_pool.primers), 2)

    def test_rv_primer_pool(self):
        self.assertTrue(
            isinstance(self.pcrarray.rv_primer_pool, PrimerPool)
        )
        self.assertEqual(len(self.pcrarray.rv_primer_pool.primers), 3)

    def test_iter(self):
        self.assertTrue(
            isinstance(self.pcrarray.iter(), types.GeneratorType)
        )

    def test_to_df(self):
        result = self.pcrarray.to_df()
        self.assertEqual(result["fw_match"].tolist(),
                         [False, False, False, False])
        self.assertEqual(result["rv_match"].tolist(),
                         [False, False, True, True])
        self.assertEqual(result["is_amplified"].tolist(),
                         [False, False, False, False])
        self.assertEqual(result["id"].tolist(), ["GAXI01000525.151.1950", 
                                                 "GAXI01000526.151.1950",
                                                 "GAXI01005455.1.1233",
                                                 "GAXI01006199.29.1525"])
        self.assertEqual(
            result["taxonomy"].tolist(),
            ["Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Arthropoda;Hexapoda;Ellipura;Collembola;Tetrodontophora bielanensis (giant springtail)",
             "Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Arthropoda;Hexapoda;Ellipura;Collembola;Tetrodontophora bielanensis (giant springtail)",
             "Bacteria;Bacteroidetes;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Chryseobacterium;Tetrodontophora bielanensis (giant springtail)",
             "Bacteria;Chlamydiae;Chlamydiae;Chlamydiales;Simkaniaceae;Candidatus Rhabdochlamydia;Tetrodontophora bielanensis (giant springtail)"]
        )

    def test_to_sql(self):
        pass


class TestMatch(unittest.TestCase):

    def test_simple_match(self):
        self.assertTrue(simple_match('AA', 'AA'))
        self.assertTrue(simple_match('aa', 'aa'))
        self.assertTrue(simple_match('Aa', 'AA'))
        self.assertTrue(simple_match('aa', 'AA'))


if __name__ == '__main__':
    unittest.main()
