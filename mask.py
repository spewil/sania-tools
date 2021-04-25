import numpy as np

ref = np.ones(shape=(1, 3))
ref[0, 0] = 2
print(ref)
seq = np.ones(shape=(3, 3))
ref = np.vstack([ref for _ in range(seq.shape[0])])

print(ref)
print(seq)

seq[np.where(seq == ref)] = np.nan

print(seq)