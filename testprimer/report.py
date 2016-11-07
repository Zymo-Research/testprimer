from __future__ import division
import os.path
import sqlite3
from collections import defaultdict, Counter

import pandas as pd


class Analysis:
    """Analysis.

    Parameters
    ----------
    sql_path : str
        SQL path containing in silico PCR result. Table must include columns
        'id', 'fw_match', 'rv_match', 'is_amplified' and 'taxonomy'.
    analyzer : Analyzer
        Analyzer must have methods 'run', 'filter' and 'output'.
    out_dir : str
        Directory the report(s) output to.
        
    Attributes
    ----------
    df : pandas.DataFrame
        Table 'testprimer' from in silico PCR SQL.
    """

    def __init__(self, sql_path, analyzer, out_dir):
        if not os.path.exists(sql_path):
            raise
        self._sql_path = sql_path
        self._analyzer = analyzer
        self.out_dir = out_dir

        with sqlite3.connect(self._sql_path) as conn:
            self.df = pd.read_sql("SELECT * FROM testprimer;", conn)

    def execute(self):
        result = self._analyzer.run(self.df)
        filtered = self._analyzer.filter(result)
        return self._analyzer.output(filtered, self.out_dir)


class TaxaCoverage:

    def run(self, df):
        data = defaultdict(Counter)
        for index, r in df.iterrows():
            is_amplified = bool(r['is_amplified'])
            taxa = r['taxonomy'].split(';')
            for i in range(len(taxa)):
                taxon = ';'.join(taxa[:i+1])
                data[taxon].update([is_amplified])
        
        coverage = pd.DataFrame(data).T.fillna(0).reset_index() \
                     .rename(columns={'index':'taxonomy', True:'match', False:'mismatch'})
        coverage['coverage'] = coverage['match'] / (coverage['match'] + coverage['mismatch'])
        return coverage

    def filter(self, coverage):
        """Result filter
        
        Only display domain level coverage and phylum level coverage.
        """
        domain = coverage[coverage['taxonomy'].str.count(';')==0]
        phylum = coverage[coverage['taxonomy'].str.count(';')==1]
        return [domain, phylum]

    def output(self, filtered, out_dir):
        writer = pd.ExcelWriter(os.path.join(out_dir, 'coverage.xlsx'))
        domain.to_excel(writer, 'domain', index=False)
        phylum.to_excel(writer, 'phylum', index=False)
        writer.save()
