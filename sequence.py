import numpy as np
import xml.etree.ElementTree as et
from Bio import SeqIO, SeqUtils, Data
from pathlib import Path

base_pairs = ["-", "N", "W", "A", "T", "C", "G"]
null_val = 9


def construct_numerization_dict():
    """
        mapping from base pairs to numbers
    """
    numerization_dict = {}
    numerization_array = np.arange(len(base_pairs))
    for i, bp in enumerate(base_pairs):
        numerization_dict[bp] = numerization_array[i]
    return numerization_dict


def construct_binarization_dict():
    """
        mapping from base pairs to binary words
    """
    binarization_dict = {}
    binarization_array = np.vstack(
        [np.zeros(len(base_pairs) - 1),
         np.eye(len(base_pairs) - 1)])
    for i, bp in enumerate(base_pairs):
        binarization_dict[bp] = binarization_array[i, :]
    return binarization_dict


NUMERIZATION_DICT = construct_numerization_dict()
BINARIZATION_DICT = construct_binarization_dict()


class SequenceCollection():
    """
        SequenceCollection is a list of
        Sequences objects from a file
        TODO
        - add combining collections
        - add multiple files
        - check for multiple reference sequences
    """
    def __init__(self, filepath):
        self.reference = None
        self.collection = self.load_from_fasta()
        self.bin_stack = self.stack_sequences("bin")
        self.num_stack = self.stack_sequences("num")
        self.seq_stack = self.stack_sequences("seq")
        self.mutation_stack = self.construct_mutation_stack(
            self.num_stack, self.reference)

    def load_from_fasta(self):
        """
        a dict {str: id : SeqRecord: record}
            SeqRecord.seq
            SeqRecord.id
        """
        collection = []
        with open(filepath) as handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        # collect sequences
        for id in seq_dict.keys():
            sequence = Sequence(id, seq_dict[id])
            if "CAP" in id:
                # sequence is reference AAV
                self.reference = sequence
            else:
                collection.append(sequence)
        return collection

    def stack_sequences(self, type):
        stack = []
        for sequence in self.collection:
            if type == "bin":
                stack.append(sequence.bin)
            elif type == "num":
                stack.append(sequence.num)
            elif type == "seq":
                stack.append(sequence.seq)
        return np.stack(stack)

    @classmethod
    def construct_mutation_stack(cls, stack, reference):
        # subtraction yields 0s where sequence truly overlaps
        diff = stack - reference.num
        mutations_only = np.ones(diff.shape) * null_val
        mutation_locs = np.argwhere(diff > 0)
        for i, j in mutation_locs:
            mutations_only[i, j] = stack[i, j]
        return mutations_only


class Sequence():
    """
    inputs:
        id: string 
        record: SeqRecord
     
    """
    def __init__(self, id, record):
        self.id = id
        self.record = record  # SeqRecord object
        self.seq = record.seq  # Seq object
        self.bin = self.binarize(self.seq)
        self.num = self.numerize(self.seq)
        try:
            self.amino_acids = self.seq.translate()
        except Data.CodonTable.TranslationError:
            self.amino_acids = None

    @classmethod
    def binarize(cls, sequence):
        """
            add binary representation of Seq
            using binarization_dict e.g.
            A: 1000
            T: 0100
            C: 0010
            G: 0001
        """
        binary_basepairs = []
        for bp in sequence:
            binary_basepairs.append(BINARIZATION_DICT[bp])
        return np.hstack(binary_basepairs)

    @classmethod
    def numerize(cls, sequence):
        """
            add numerical representation of Seq 
            using numerization dict e.g.
            A: 1
            T: 2
            C: 3
            G: 4
        """
        numerical_basepairs = []
        for bp in sequence:
            numerical_basepairs.append(NUMERIZATION_DICT[bp])
        return np.array(numerical_basepairs)


# not used for now
def load_from_xml(filepath):
    output = {}
    tree = et.parse(filepath)
    root = tree.getroot()
    for sequence_element in root:
        for name_element in sequence_element.findall("INSDSeq_locus"):
            sequence = sequence_element.find("INSDSeq_sequence").text
            output[name_element.text] = sequence
    return output


# folder = '/Users/spencerw/Dropbox (UCL)/Murray Lab/Geneious/Alba data example/'
folder = "data"
fasta_filename = 'prelim_capsids.fasta'
filepath = Path(folder) / Path(fasta_filename)
sc = SequenceCollection(filepath)

print(f"Found {len(sc.collection)} records in file")
i = 31
print(f"Sequence {i} contains {len(sc.collection[i].seq)} base pairs")
print(f"Sequence {i} contains {len(sc.collection[i].amino_acids)} amino acids")

print(f"unique amino acids: {set(sc.collection[i].seq.translate())}")
print(f"unique amino acids: {set(sc.collection[i].amino_acids)}")

print(f"Reference contains {len(sc.reference.seq)} base pairs")
print(f"Shape of binary stack {sc.bin_stack.shape}")
print(f"Shape of numerary stack {sc.num_stack.shape}")
print(f"Shape of raw stack {sc.seq_stack.shape}")
print(f"Shape of mutation stack {sc.mutation_stack.shape}")
