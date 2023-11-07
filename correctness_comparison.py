# -*- coding: utf-8 -*-
"""correctness_comparison.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1CgadXAfr3kr0i6FVh67RbJFCTjXoFtYu
"""

import matplotlib.pyplot as plt
import numpy as np

exact_results = [3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 8, 7, 7, 7, 7, 7, 9, 8, 9, 9, 12, 10, 10, 11, 10, 11, 10, 13]
algorithm_results = [3, 3, 4, 5, 5, 6, 6, 7, 7, 6, 8, 5, 7, 7, 7, 7, 9, 7, 9, 9, 12, 10, 10, 11, 8, 11, 9, 10]
n_values = range(2, len(algorithm_results) + 2)

bar_width = 0.35
index = np.arange(len(n_values))

# Create a bar plot
plt.figure(figsize=(10, 6))
plt.bar(index, algorithm_results, bar_width, color='red', label='Algorithm Results')
plt.bar(index + bar_width, exact_results, bar_width, color='green', label='Exact Results')
plt.xlabel('n')
plt.ylabel('Results')
plt.title('Algorithm Results vs Exact Results')
plt.xticks(index + bar_width / 2, n_values)
plt.legend()
plt.grid(True)