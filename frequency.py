from collections import Counter
from sequence import sequences

# for all sequences, how many sequences are identical?
# get unique sequences, count them, plot counts and sequences


def list_to_word(l):
    w = ""
    for li in l:
        w += li
    return w


# convert sequences to words?
words = []
for seq in sequences.keys():
    s = sequences[seq]["seq"]
    words.append(s._data)

counts = Counter(words)

print(counts.values())
print(sum(counts.values()))
