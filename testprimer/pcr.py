import os.path
import sqlite3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd


class Template:
    """PCR template.

    Parameters
    ----------
    seqrecord : Bio.SeqRecord
        Mostly used with element in Bio.SeqIO.parse iterator.

    Attributes
    ----------
    seq : Bio.Seq
        Unambiguous DNA/RNA sequence. Gap supported.
    description : str
        FASTA header line.
    id : str
        Unique sequence identifier.
    taxonomy : str
        Semicolon-separated, hierarchical taxonomic classification.
    """

    def __init__(self, seqrecord):
        if not isinstance(seqrecord, SeqRecord):
            raise TypeError
        self._seqrecord = seqrecord
        self.seq = self._seqrecord.seq
        self.description = self._seqrecord.description
        self.id = self._seqrecord.id
        self.taxonomy = self.description.split(' ', 1)[1]#.split(';')


class Primer(Seq):
    """PCR primer.
    
    A copy of Bio.Seq class.
    """
    pass
    

class PrimerPool:
    """PCR primer pool.

    A collection of DNA primers for one strand. Must share the same
    starting and ending coordinates.
    
    Parameters
    ----------
    name : str
        Primer pool name.
    start : int
        1-based genomic coordinate indicating where the primer starts.
    end : int
        1-based genomic coordinate indicating where the primer ends.
    primers : list of Primer or list of Seq, or None, default None
        Primers for one strand. Specified upon instance creation or by
        instance method later.
    """

    def __init__(self, name, start, end, primers=None):
        self.name = name
        self.start = start
        self.end = end
        if not primers:
            self.primers = primers
        else:
            if not isinstance(primers, list):
                raise TypeError
            else:
                self.primers = primers

    def parse(self, seqs_str):
        """Function to add/update primers.

        Add primers to primer pool by parsing the given string in which
        primers are separated by line break. If primers already exist,
        they will be replaced.

        Parameters
        ----------
        seqs_str : str
            Primers formatted in one string separated by '\n'.   

        Examples
        --------
        GTGCCAGCAGTCGCGGTAA
        GTGCCAGCAGGAGCGGTAA
        GTGCCACCAGCCGCGGTAA
        GTGCCAGAAGTCTCGGTAA
        GTGCCAGAAGCGTCGGTAA
        GTGCCAGAAGCCTCGGTAA
        """

        self.primers = []
        for seq in seqs_str.strip().split('\n'):
            self.primers.append(Primer(seq.strip()))


class PCRMatch:

    def __init__(self, id, taxonomy, fw_match, rv_match, is_amplified):
        self.id = id
        self.taxonomy = taxonomy
        self.fw_match = fw_match
        self.rv_match = rv_match
        self.is_amplified = is_amplified


class PCR:
    """in silico PCR.
    
    Simulate PCR on a template sequence with forward primer pool and reverse
    primer pool following certain rules defining primer/template match.

    Parameters
    ----------
    template : Template
        PCR template. Support gapped/ungapped DNA/RNA sequence.
    fw_primer_pool : PrimerPool
        Forward primer pool containing one or multiple primers.
    rv_primer_pool : PrimerPool
        Reverse primer pool containing one or multiple primers.
    find_match : function
        Function to determine if primer could bind to template. Only
        return boolean value.
    """

    def __init__(self, template, fw_primer_pool, rv_primer_pool, find_match):
        if not isinstance(template, Template):
            raise TypeError
        if not isinstance(fw_primer_pool, PrimerPool):
            raise TypeError
        if not isinstance(rv_primer_pool, PrimerPool):
            raise TypeError
        self.template = template
        self.fw_primer_pool = fw_primer_pool
        self.rv_primer_pool = rv_primer_pool
        self._find_match = find_match

    def cut_forward_primer_from_template(self):
        start = self.fw_primer_pool.start - 1
        end = self.fw_primer_pool.end
        return self.template.seq[start:end]

    def cut_reverse_primer_from_template(self):
        start = self.rv_primer_pool.start - 1
        end = self.rv_primer_pool.end
        return self.template.seq[start:end]

    def run(self):
        """Run in silico PCR.

        """

        GAP = '-'  # gap placeholder
        UNKNOWN = '.'  # leading and trailing placeholder 

        # cut primer from template and back transcribe to DNA
        fw_tmplt = self.cut_forward_primer_from_template().back_transcribe()
        rv_tmplt = self.cut_reverse_primer_from_template().back_transcribe()

        # ungap
        fw_tmplt = fw_tmplt.ungap(UNKNOWN).ungap(GAP)
        rv_tmplt = rv_tmplt.ungap(UNKNOWN).ungap(GAP)

        # test match
        for fw in self.fw_primer_pool.primers:
            if self._find_match(str(fw), str(fw_tmplt)):
                fw_match = True
                break
            fw_match = False

        for rv in self.rv_primer_pool.primers:
            if self._find_match(str(rv.reverse_complement()), str(rv_tmplt)):
                rv_match = True
                break
            rv_match = False

        is_amplified = fw_match and rv_match
        return PCRMatch(self.template.id, self.template.taxonomy, 
                        fw_match, rv_match, is_amplified)


class PCRArray:

    def __init__(self, fasta_path, fw_path, rv_path):
        self._fasta_path = fasta_path
        self._fw_path = fw_path
        self._rv_path = rv_path

    def fw_primer_pool(self):
        '''Assuming primer pool in this format:

        #515f 11895-13861
        GTGCCAGCAGTCGCGGTAA
        GTGCCAGCAGGAGCGGTAA
        GTGCCACCAGCCGCGGTAA
        GTGCCAGAAGTCTCGGTAA
        GTGCCAGAAGCGTCGGTAA
        GTGCCAGAAGCCTCGGTAA
        '''
        with open(self._fw_path, 'r') as handle:
            lines = [line.strip() for line in handle if line.strip()]
            header = lines[0].lstrip('#').strip()
            primers = map(Seq, lines[1:])
            name, position = header.split()
            start, end  = map(int, position.split('-'))
        return PrimerPool(name, start, end, primers)

    def rv_primer_pool(self):
        with open(self._rv_path, 'r') as handle:
            lines = [line.strip() for line in handle if line.strip()]
            header = lines[0].lstrip('#').strip()
            primers = map(Seq, lines[1:])
            name, position = header.split()
            start, end  = map(int, position.split('-'))
        return PrimerPool(name, start, end, primers)

    def iter(self):
        for seqrecord in SeqIO.parse(self._fasta_path, 'fasta'):
            template = Template(seqrecord)
            fw_primer_pool = self.fw_primer_pool()
            rv_primer_pool = self.rv_primer_pool()
            pcr = PCR(template, fw_primer_pool, rv_primer_pool, simple_match)
            yield pcr.run()

    def to_df(self):
        data = [pcrmatch.__dict__ for pcrmatch in self.iter()]
        return pd.DataFrame(data)

    def to_sql(self, filename, out_dir):
        df = self.to_df()
        with sqlite3.connect(os.path.join(out_dir, filename)) as conn:
            df.to_sql('testprimer', conn, if_exists='replace', index=False)
        return


def simple_match(seq1, seq2):
    return seq1.upper() == seq2.upper() 


def pcr(fasta_path, fw_path, rv_path, filename, out_dir):
    '''Main module entrance'''

    pcrarray = PCRArray(fasta_path, fw_path, rv_path)
    pcrarray.to_sql(filename, out_dir)
    return
