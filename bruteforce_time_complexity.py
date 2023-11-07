# -*- coding: utf-8 -*-
"""bruteforce_time_complexity.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1CgadXAfr3kr0i6FVh67RbJFCTjXoFtYu
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

data = {
    2: [1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    3: [4, 2, 3, 3, 2, 2, 3, 3, 3, 4, 2, 2, 4, 2, 2, 4, 2, 3, 3, 3],
    4: [6, 8, 7, 8, 6, 7, 7, 6, 7, 7, 10, 4, 7, 6, 5, 7, 8, 5, 7, 7],
    5: [21, 14, 13, 11, 17, 12, 15, 20, 9, 12, 12, 22, 11, 11, 12, 15, 12, 18, 12, 12],
    6: [25, 25, 23, 26, 36, 129, 31, 51, 111, 52, 35, 31, 24, 23, 24, 23, 23, 24, 24],
    7: [50, 49, 47, 39, 48, 42, 50, 46, 56, 53, 40, 41, 78, 48, 38, 38, 39, 39, 39],
    8: [43, 43, 43, 42, 42, 45, 42, 41, 42, 44, 46, 46, 45, 45, 31, 30, 31, 44, 40],
    9: [168, 220, 185, 174, 164, 165, 170, 163, 161, 174, 163, 163, 182, 165, 164, 184, 164, 164, 169],
    10: [368, 286, 288, 311, 287, 280, 271, 292, 282, 282, 313, 278, 276, 274, 291, 280, 284, 283, 274],
    11: [521, 502, 42334, 862, 722, 720, 3774, 720, 760, 689, 683, 683, 710, 754, 738, 725, 661, 714, 704],
    12: [1510, 1267, 1234, 1331, 1221, 1307, 1397, 40884, 1860, 1879, 5930, 1826, 1816, 1586, 1266, 1272, 1255, 1610, 1825],
    13: [1308, 1158, 1171, 1240, 1145, 1171, 1160, 1165, 1125, 1141, 1178, 1279, 1170, 1165, 1211, 1202, 1149, 1181, 1245],
    14: [1345, 1200, 1162, 1541, 1238, 1184, 1195, 1214, 1177, 1175, 1172, 1180, 1171, 1195, 1162, 5585, 1752, 1757, 1709],
    15: [2062, 1988, 1974, 1953, 1985, 1996, 2456, 7037, 3053, 2956, 2663, 2011, 2005, 2745, 2035, 2007, 1978, 1985, 1990],
    16: [2167, 2198, 2141, 2066, 2053, 2114, 2176, 2047, 2158, 2200, 2251, 2143, 2185, 2085, 2063, 2205, 2025, 1994, 2160],
    17: [3876, 3930, 3966, 3956, 4062, 4243, 3950, 4076, 4691, 4050, 4152, 4172, 3979, 4258, 3991, 4002, 4144, 3976, 5973],
    18: [5902, 3597, 3587, 3662, 3761, 3589, 3812, 3739, 3704, 3657, 3570, 3592, 3656, 3553, 3590, 3583, 3560, 3613, 3628],
    19: [3957, 3868, 5264, 3842, 4021, 3796, 3777, 3953, 3871, 4042, 3774, 3830, 3827, 3828, 3892, 3802, 3604, 3792, 3831],
    20: [8125, 7870, 7929, 8437, 7926, 7934, 7601, 7728, 7527, 7683, 7462, 7583, 7362, 8983, 7898, 7207, 7444, 7424, 8677]
}

results = {}

for n, values in data.items():
    mean = sum(values) / len(values)
    sd = math.sqrt(sum((x - mean)**2 for x in values) / len(values))
    sm = sd / math.sqrt(n)
    M1 = mean + 2.093 * sm
    M2 = mean + 2.093 * sm
    results[n] = (M1, M2)

x_values = list(data.keys())
y_values_M1 = [result[0] for result in results.values()]

plt.scatter(x_values, y_values_M1, label='M1')
plt.xlabel('N values')
plt.ylabel('M1 values')
plt.title('Scatter Plot of M1 values for different N values')
plt.legend()

regression_x = np.array(x_values).reshape(-1, 1)
regression_y_M1 = np.array(y_values_M1)

regressor_M1 = LinearRegression()
regressor_M1.fit(regression_x, regression_y_M1)
regression_line_M1 = regressor_M1.predict(regression_x)

slope_M1 = regressor_M1.coef_[0]
print("Slope of the fitted line for M1:", slope_M1)

plt.plot(x_values, regression_line_M1, color='red')

plt.xticks(x_values)
plt.legend()
plt.show()