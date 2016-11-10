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
        Path to SQL file containing in silico PCR result. Table must include
        columns 'id', 'fw_match', 'rv_match', 'is_amplified' and 'taxonomy'.
    analyzer : Analyzer
        Analyzer must have methods 'run', 'filter' and 'output'.
    out_dir : str
        Directory the report(s) output to.
        
    Attributes
    ----------
    df : pandas.DataFrame
        Table 'testprimer' from SQL file.
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

    OUTPUT_FILENAME = 'coverage.xlsx'

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
        """Result selector.
        
        Only display domain level, phylum level and human disease related
        pathogen coverage.
        """
        domain = coverage[coverage['taxonomy'].str.count(';')==0]
        phylum = coverage[coverage['taxonomy'].str.count(';')==1]

        # human disease related pathogens
        genus = coverage[(coverage['taxonomy'].str.startswith('Bacteria')) & (coverage['taxonomy'].str.count(';')==5)]
        with open('pathogens.txt', 'r') as handle:
            pathogenlist = map(lambda x: x.strip(), handle.readlines())

        data = defaultdict(list)
        for candidate in pathogenlist:
            row = genus[genus['taxonomy'].str.endswith(candidate)]
            data['pathogen'].append(candidate)
            
            if row.shape[0] != 0:    
                data['taxonomy'].append(row.iloc[0]['taxonomy'])
                data['mismatch'].append(row.iloc[0]['mismatch'])
                data['match'].append(row.iloc[0]['match'])
                data['coverage'].append(row.iloc[0]['coverage'])
            else:
                data['taxonomy'].append(None)
                data['mismatch'].append(None)
                data['match'].append(None)
                data['coverage'].append(None)
        pathogen = pd.DataFrame(data, columns=['pathogen', 'taxonomy', 'mismatch', 'match', 'coverage'])

        return [domain, phylum, pathogen]

    def output(self, filtered, out_dir):
        writer = pd.ExcelWriter(os.path.join(out_dir, OUTPUT_FILENAME))
        domain, phylum, pathogen = filtered
        domain.to_excel(writer, 'domain', index=False)
        phylum.to_excel(writer, 'phylum', index=False)
        pathogen.to_excel(writer, 'pathogen', index=False)
        writer.save()
