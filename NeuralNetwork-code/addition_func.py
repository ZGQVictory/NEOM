import numpy as np
import os
from arg_parse import get_args

args = get_args()


def normalize(freq_mat, length):
    # freq_mat: (aa_type, length)
    for i in range(length):
        loc_sum = np.sum(freq_mat[:, i])
        for j in range(20):
            freq_mat[j, i] = freq_mat[j, i] / loc_sum
    return freq_mat


def get_freq(direction, filename, length):
    aa = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    num_pep = 0
    pos = np.zeros((20, length))  ##line --> 20aa, column --> position
    dictionary = dict(zip(aa, np.linspace(0, 19, 20, dtype=int)))
    with open(os.path.join(direction, filename), "r") as f:
        norm = True  # justify if the sequence is normal (length is 9, all aa in the dictionary)
        for line in f:
            line = line.replace("\n", "")
            line = line.replace(" ", "")
            if len(line) == length and line.isupper() == True:
                for i in range(length):
                    if line[i] in dictionary.keys():
                        norm = True
                        continue
                    else:
                        norm = False
                        break
                if norm == True:
                    num_pep += 1
                    for i in range(length):
                        nrow = dictionary[line[i]]
                        pos[nrow][i] += 1
                else:
                    continue
    return normalize(pos, length), num_pep
