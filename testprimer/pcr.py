from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Template:

    def __init__(self, seqrecord):
        if not isinstance(seqrecord, SeqRecord):
            raise TypeError
        self._seqrecord = seqrecord
        self.seq = self._seqrecord.seq
        self.description = self._seqrecord.description
        self.id = self._seqrecord.id
        self.taxonomy = self.description.split(' ', 1)[1].split(';')


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

    def __init__(self, template, fw_match, rv_match, is_amplified):
        self.template = template
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
        return PCRMatch(self.template, fw_match, rv_match, is_amplified)


def simple_match(seq1, seq2):
    return seq1.upper() == seq2.upper() 


if __name__ == '__main__':
    pass
