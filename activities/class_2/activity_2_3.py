"""Calculations for the Third Class 2 Activity in REB, The Course"""
# import libraries
import numpy as np

# transpose of the stoichiometric coefficient matrix
nu_transpose = np.array([[-1, -2, 1, 0, 0, 0],
                         [0, -3, 1, -1, 1, 0],
                         [-1, 1, 0, 1, -1, 0],
                         [-1, -3, 0, 0, 1, 1],
                         [0, -4, 0, -1, 2, 1],
                         [0, 1, 1, 0, -1, -1]])
N_independent = np.linalg.matrix_rank(nu_transpose)
print()
print(f"Number of independent reactions: {N_independent}")
print()

# working matrix with the first two rows of nu_transpose
W1 = nu_transpose[0:2, :]
rank_W1 = np.linalg.matrix_rank(W1)
print(f"Rank of W1: {rank_W1}")
print()

# working matrix with the first three rows of nu_transpose
W2 = nu_transpose[0:3, :]
rank_W2 = np.linalg.matrix_rank(W2)
print(f"Rank of W2: {rank_W2}")
print()

# working matrix with the first, second, and fourth rows of nu_transpose
W3 = nu_transpose[[0, 1, 3], :]
rank_W3 = np.linalg.matrix_rank(W3)
print(f"Rank of W3: {rank_W3}")
print()
