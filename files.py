import numpy as np
import xml.etree.ElementTree as et
from Bio import SeqIO
from pathlib import Path

base_pairs = ["-", "N", "W", "A", "T", "C", "G"]


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


binarization_dict = construct_binarization_dict()
numerization_dict = construct_numerization_dict()

folder = '/Users/spencerw/Dropbox (UCL)/Murray Lab/Geneious/Alba data example/'
# filename = 'nucleotide_alignment.xml'
# filepath = folder + filename

fasta_filename = 'nucleotide_alignment.fasta'
filepath = Path(folder + fasta_filename)
seq_dict = load_from_fasta(filepath)
print("Found %i records in file" % len(seq_dict.keys()))


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


sequences, reference_sequences = construct_sequence_dict(seq_dict)
print(f"found {len(reference_sequences.keys())} reference sequences")
print(f"found {len(sequences.keys())} experimental sequences")


def stack_sequences(sequences, type="num"):
    stack = []
    for id in sequences.keys():
        stack.append(sequences[id][type])
    return np.stack(stack)


stack = stack_sequences(sequences)
stack_with_reference = np.vstack(
    [list(reference_sequences.values())[0]["num"], stack])

import matplotlib
import matplotlib.pyplot as plt
from utils import plot
matplotlib.use("TkAgg")

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
img = ax.imshow(stack_with_reference, aspect='auto', cmap='tab10')
ax.set_yticks(
    [i - 0.5 for i in list(range(stack_with_reference.shape[0] + 1))])
ax.set_yticklabels(["AAV6"] + list(range(1, stack_with_reference.shape[0])) +
                   [""])
ax.yaxis.grid(True, which='major', color="k")
ax.set_ylabel("sequence")
ax.set_xlabel("nucleotide")
cbar = fig.colorbar(img, ticks=range(10))
labels = base_pairs + [" ", " ", " "]
cbar.ax.set_yticklabels(labels)

stack = stack_sequences(sequences, type="bin")
stack_with_reference = np.vstack(
    [list(reference_sequences.values())[0]["bin"], stack])

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
img = ax.imshow(stack_with_reference, aspect='auto', cmap='binary')
ax.set_yticks(
    [i - 0.5 for i in list(range(stack_with_reference.shape[0] + 1))])
ax.set_yticklabels(["AAV6"] + list(range(1, stack_with_reference.shape[0])) +
                   [""])
ax.yaxis.grid(True, which='major', color="r")
ax.set_ylabel("sequence")
ax.set_xlabel("nucleotide")
plt.show()