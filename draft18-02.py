## I = |(M1/4(sc)*e_s)^T *R* ((M1/4(inc)*e_i)|^2
## https://magnetism.eu/esm/2018/slides/vavassori-practical.pdf
## https://onlinelibrary.wiley.com/doi/book/10.1002/9780470060193
## https://pypolar.readthedocs.io/en/latest/06-Mueller-Matrices.html

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

theta = np.radians(180)
e_i = np.array([1, np.cos(2*theta), np.sin(2*theta), 0])  # LP in Stokes representation

# General retarder
t = np.radians(20) #angle of slow axis
delta = 45 #phase difference between the fast and slow axis, 90 for L/4, 45 for L/2
cos2theta = np.cos(2 * t)
sin2theta = np.sin(2 * t)
cos_delta = np.cos(delta)
sin_delta = np.sin(delta)
    
M_L4 = np.array([
        [1, 0, 0, 0],
        [0, cos2theta**2 + sin2theta**2 * cos_delta, cos2theta * sin2theta * (1 - cos_delta), sin2theta * sin_delta],
        [0, cos2theta * sin2theta * (1 - cos_delta), cos2theta**2 * cos_delta + sin2theta**2, -cos2theta * sin_delta],
        [0, -sin2theta * sin_delta, cos2theta * sin_delta, cos_delta]])

# linear polarizer (horizontal transmission)
M_LP = (1/2) * np.array([
    [1, 1, 0, 0],
    [1, 1, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]])

# firstly l/4
e_i2 = np.matmul(M_LP, np.matmul(M_L4, e_i))
print("Output Stokes vector after passing through linear polarizer (horizontal) and L/4:", e_i2)
