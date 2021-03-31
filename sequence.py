import numpy as np
import xml.etree.ElementTree as et
from Bio import SeqIO
from pathlib import Path

base_pairs = ["-", "N", "W", "A", "T", "C", "G"]
null_val = 9


class Sequence(SeqIO):
    def __init__(self):
        super().__init__()
        print(dir(self))


def load_from_xml(filepath):
    output = {}
    tree = et.parse(filepath)
    root = tree.getroot()
    for sequence_element in root:
        for name_element in sequence_element.findall("INSDSeq_locus"):
            sequence = sequence_element.find("INSDSeq_sequence").text
            output[name_element.text] = sequence
    return output


def binarize_sequence(sequence):
    binary_basepairs = []
    for bp in sequence:
        binary_basepairs.append(binarization_dict[bp])
    return np.hstack(binary_basepairs)


def numerize_sequence(sequence):
    numerical_basepairs = []
    for bp in sequence:
        numerical_basepairs.append(numerization_dict[bp])
    return np.array(numerical_basepairs)


def load_from_fasta(filepath):
    """
        returns a dict {id : SeqRecord}
        SeqRecord.seq
        SeqRecord.id
        ...
    """
    with open(filepath) as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def construct_numerization_dict():
    numerization_dict = {}
    numerization_array = np.arange(len(base_pairs))
    for i, bp in enumerate(base_pairs):
        numerization_dict[bp] = numerization_array[i]
    return numerization_dict


def construct_binarization_dict():
    binarization_dict = {}
    binarization_array = np.vstack(
        [np.zeros(len(base_pairs) - 1),
         np.eye(len(base_pairs) - 1)])
    for i, bp in enumerate(base_pairs):
        binarization_dict[bp] = binarization_array[i, :]
    return binarization_dict


def construct_sequence_dict(seq_dict):
    sequences = {}
    reference_sequences = {}
    for id in seq_dict.keys():
        record = seq_dict[id]
        if "CAP" in id:
            print(id)
            reference_sequences[id] = {}
            reference_sequences[id]["seq"] = record.seq
            reference_sequences[id]["bin"] = binarize_sequence(record.seq)
            reference_sequences[id]["num"] = numerize_sequence(record.seq)
        else:
            sequences[id] = {}
            sequences[id]["seq"] = record.seq
            sequences[id]["bin"] = binarize_sequence(record.seq)
            sequences[id]["num"] = numerize_sequence(record.seq)
    return sequences, reference_sequences


def stack_sequences(sequences, type="num"):
    stack = []
    for id in sequences.keys():
        stack.append(sequences[id][type])
    return np.stack(stack)


binarization_dict = construct_binarization_dict()
numerization_dict = construct_numerization_dict()

folder = '/Users/spencerw/Dropbox (UCL)/Murray Lab/Geneious/Alba data example/'
# filename = 'nucleotide_alignment.xml'
# filepath = folder + filename

fasta_filename = 'nucleotide_alignment.fasta'
filepath = Path(folder + fasta_filename)
seq_dict = load_from_fasta(filepath)
print("Found %i records in file" % len(seq_dict.keys()))

sequences, reference_sequences = construct_sequence_dict(seq_dict)
print(f"found {len(reference_sequences.keys())} reference sequences")
print(f"found {len(sequences.keys())} experimental sequences")

# converted sequences
binary_stack = stack_sequences(sequences, type="bin")
numerical_stack = stack_sequences(sequences, type="num")
raw_stack = stack_sequences(sequences, type="seq")

# visualizing mutations only
viz_type = "num"
reference = list(reference_sequences.values())[0][viz_type]
# there are no zeros in the reference, only 3,4,5,6 == A,T,C,G
stack = stack_sequences(sequences, type=viz_type)
# subtraction yields 0s where sequence truly overlaps
diff = stack - reference
mutations_only = np.ones(diff.shape) * null_val
mutation_locs = np.argwhere(diff > 0)
for i, j in mutation_locs:
    mutations_only[i, j] = stack[i, j]

S = Sequence()