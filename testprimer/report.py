from __future__ import division
import os.path
import sqlite3

import pandas as pd


class Analysis:

    def __init__(self, sqlpath, get_stats):
        self._sqlpath = sqlpath
        self._get_stats = get_stats
        if not os.path.exists(self._sqlpath):
            raise

    def df(self):
        with sqlite3.connect(self._sqlpath) as conn:
            df = pd.read_sql("SELECT * FROM testprimer;", conn)
        return df

    def process(self):
        pass
        

def taxa_coverage(df):
    from collections import defaultdict, Counter

    d = defaultdict(Counter)
    for index, r in df.iterrows():
        is_amplified = bool(r['is_amplified'])
        taxa = r['taxonomy'].split(';')
        for i in range(len(taxa)):
            taxon = ';'.join(taxa[:i+1])
            d[taxon].update([is_amplified])
    return d


# class SQLParser:

    # def __init__(self, sql_path):
        # self.sql_path = sql_path

    # def parse(self):
        # with sqlite3.connect(self.sql_path) as conn:
            # df = pd.read_sql("SELECT * FROM testprimer;", conn)
        # return df


# class Analysis:

    # def __init__(self, df):
        # self._df = df

    # @property
    # def coverage(self):
        # d = defaultdict(Counter)
        # for i, r in df.iterrows():
            # is_amplified = bool(r['is_amplified'])
            # ranks = r['taxonomy'].split(';')
            # for i in range(len(ranks)):
                # rank = ';'.join(ranks[:i+1])
                # d[rank].update([is_amplified])
        # return d
