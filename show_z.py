import numpy as np
import matplotlib.pyplot as plt
from utils_node import read_stroke


img = np.zeros((128, 128))+255
file = read_stroke("Sample_points.txt")
zz = []
for i in range(len(file)):
    stroke = file[i]
    for (x, y, z) in stroke:
        zz.append(z)

plt.plot(zz)
plt.show()
