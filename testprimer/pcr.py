from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd


class Template:

    def __init__(self, seqrecord):
        if not isinstance(seqrecord, SeqRecord):
            raise TypeError
        self._seqrecord = seqrecord
        self.seq = self._seqrecord.seq
        self.description = self._seqrecord.description
        self.id = self._seqrecord.id
        self.taxonomy = self.description.split(' ', 1)[1]#.split(';')


class Primer(Seq):
    '''A copy of Seq class'''
    pass
    

class PrimerPool:

    def __init__(self, name, start, end, primers=None):
        '''
        start, end is 1-based
        '''
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
        self.primers = []
        for seq in seqs_str.strip().split('\n'):
            self.primers.append(Primer(seq.strip()))


# class AmplificationError(Exception):
    # pass


class PCRMatch:

    def __init__(self, id, taxonomy, fw_match, rv_match, is_amplified):
        self.id = id
        self.taxonomy = taxonomy
        self.fw_match = fw_match
        self.rv_match = rv_match
        self.is_amplified = is_amplified


class PCR:

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

    # def to_sql(self, out_dir='/mnt'):
        # df = self.to_df()


def simple_match(seq1, seq2):
    return seq1.upper() == seq2.upper() 


if __name__ == '__main__':
    pass
