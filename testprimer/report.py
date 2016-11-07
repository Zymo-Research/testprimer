from __future__ import division
import os.path
import sqlite3
from collections import defaultdict, Counter

import pandas as pd


class Analysis:

    def __init__(self, sql_path, analyzer):
        if not os.path.exists(sql_path):
            raise
        self._sql_path = sql_path
        self._analyzer = analyzer

        with sqlite3.connect(self._sql_path) as conn:
            self.df = pd.read_sql("SELECT * FROM testprimer;", conn)

    def execute(self):
        result = self._analyzer.run(self.df)
        return self._analyzer.filter(result)


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
        
        Only display phylum level coverage report.
        """
        return coverage[coverage['taxonomy'].str.count(';')==1]
