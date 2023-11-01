import numpy as np
from itertools import product
combinations = product("ACGT", repeat=3)
combinations_list = ["".join(comb) for comb in combinations]
stop_codons = ['TAA', 'TAG', 'TGA']
sense_codons = [comb for comb in combinations_list if comb not in stop_codons]

import numpy as np

def compare_strings(strings):
    n = len(strings)
    matrix = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            str1 = np.array(list(strings[i]))
            str2 = np.array(list(strings[j]))
            diff_count = np.sum(str1 != str2)            
            if diff_count == 1:
                matrix[i, j] = 1
                matrix[j, i] = 1
    diag_indices = np.diag_indices(n)
    row_sums = np.sum(matrix, axis=1)
    matrix[diag_indices] = -row_sums           
    return matrix

result = compare_strings(sense_codons)
np.set_printoptions(threshold=np.inf)































from itertools import product
combinations = product("ACGT", repeat=3)
combinations_list = ["".join(comb) for comb in combinations]
stop_codons = ['TAA', 'TAG', 'TGA']
combinations_list = [comb for comb in combinations_list if comb not in stop_codons]
combinations_list.remove('ATG')
combinations_list.insert(0, 'ATG')

import random
combinations_list = [combinations_list]
for _ in range(10):
    new_codon = 'TAA'
    while new_codon in stop_codons:
        combinations_list.append(combinations_list[-1].copy())
        nucleotides = ['A', 'C', 'G', 'T']
        site_index = random.randrange(61)
        position_index = random.randrange(3)
        substitution_index = random.randrange(3)
        old_codon = list(combinations_list[-1][site_index])
        nucleotides.remove(old_codon[position_index])
        old_codon[position_index] = nucleotides[substitution_index]
        new_codon = ''.join(old_codon)
    combinations_list[-1][site_index] = new_codon

for i,seq in enumerate(combinations_list):
    print(f'>seq_{i}')
    print(''.join(seq))

