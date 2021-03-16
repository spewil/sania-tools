import numpy as np
import xml.etree.ElementTree as et
from Bio import SeqIO

folder = '/Users/spencerw/Dropbox (UCL)/Murray Lab/Geneious/Alba data example/'
filename = 'nucleotide_alignment.xml'
filepath = folder + filename


def load_from_xml(filepath):
    output = {}
    tree = et.parse(filepath)
    root = tree.getroot()
    for sequence_element in root:
        for name_element in sequence_element.findall("INSDSeq_locus"):
            sequence = sequence_element.find("INSDSeq_sequence").text
            output[name_element.text] = sequence
    return output


binarization_array = np.eye(4)
base_pairs = ["A", "T", "C", "G"]
binarization_dict = {}
for i, bp in enumerate(base_pairs):
    binarization_dict[bp] = binarization_array[:, i]

print(binarization_dict)


def binarize_sequence(sequence):
    binary_basepairs = []
    for bp in sequence:
        binary_basepairs.append(binarization_dict[bp])
    return np.hstack(binary_basepairs)


d = load_from_xml(filepath)
num_seqs = len(d.values())
print(d.keys())
# print(binarize_sequence(s).shape)
# binarized_sequences = np.empty()

fasta_filename = 'nucleotide_alignment.fasta'
filepath = folder + filename
seq = SeqIO.parse(filepath, "fasta")
records = list(seq)
print(seq)
print("Found %i records" % len(records))
