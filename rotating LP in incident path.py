## I = |(M1/4(sc)*e_s)^T *R* ((M1/4(inc)*e_i)|^2
## https://magnetism.eu/esm/2018/slides/vavassori-practical.pdf
## https://onlinelibrary.wiley.com/doi/book/10.1002/9780470060193
## https://pypolar.readthedocs.io/en/latest/06-Mueller-Matrices.html

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

angle = np.radians(180)
e_i = np.array([1, np.cos(2*angle), np.sin(2*angle), 0])  # LP in Stokes representation

# general retarder
t = np.radians(20) #angle of slow axis, can be rotated
delta = 90 #phase difference between the fast and slow axis, 90 for L/4, 45 for L/2
cos2theta = np.cos(2 * t)
sin2theta = np.sin(2 * t)
cos_delta = np.cos(delta)
sin_delta = np.sin(delta)
    
M_R = np.array([
        [1, 0, 0, 0],
        [0, cos2theta**2 + sin2theta**2 * cos_delta, cos2theta * sin2theta * (1 - cos_delta), sin2theta * sin_delta],
        [0, cos2theta * sin2theta * (1 - cos_delta), cos2theta**2 * cos_delta + sin2theta**2, -cos2theta * sin_delta],
        [0, -sin2theta * sin_delta, cos2theta * sin_delta, cos_delta]])

# iterate over every orientation of LINEAR POLARIZER
T_values = np.radians(np.arange(0, 361, 1))  # angles from 0 to 360 in steps of 1
stokes_vectors = []

for T in T_values:
    cos2Theta = np.cos(2 * T)
    sin2Theta = np.sin(2 * T)

    M_LP = (1/2) * np.array([
        [1, cos2Theta, sin2Theta, 0],
        [cos2Theta, cos2Theta**2, cos2Theta*sin2Theta, 0],
        [sin2Theta, cos2Theta*sin2Theta, sin2Theta**2, 0],
        [0, 0, 0, 0]
    ])

    output = np.matmul(M_LP, e_i)
    stokes_vectors.append(output)

df = pd.DataFrame(stokes_vectors, columns=['S0', 'S1', 'S2', 'S3'])
df['T (degrees)'] = np.degrees(T_values) 

# Plot S0 vs T
plt.figure(figsize=(10, 6))
plt.plot(df['T (degrees)'], df['S0'], label="S0", color='b')
plt.xlabel("Linear Polarizer Angle (degrees)")
plt.ylabel("S0 (Intensity)")
