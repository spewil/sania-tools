from collections import Counter
from sequence import SequenceCollection, Sequence, sc
import numpy as np

# for all sequences, how many sequences are identical?
# get unique sequences, count them, plot counts and sequences

# top plot -- bar plot of most frequent sequences


def list_to_word(l):
    w = ""
    for li in l:
        w += li
    return w


# convert sequences to words?
words = [s.seq._data for s in sc.collection]


def count_dict(sequences):
    d = Counter(sequences)
    return sorted(d.items(), key=lambda item: item[1])


sorted_words = count_dict(words)[::-1][:10]
unique_sequences = []
counts = []
for word, count in sorted_words:
    unique_sequences.append((Sequence.numerize(word)))
    counts.append(count)
unique_stack = np.stack(unique_sequences)
unique_mutations = SequenceCollection.construct_mutation_stack(
    unique_stack, sc.reference)

for i in range(5):
    print(np.nonzero(unique_mutations[i] - 9)[0])
    print([
        unique_mutations[i, k] for k in np.nonzero(unique_mutations[i] - 9)[0]
    ])
