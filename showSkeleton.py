import numpy as np
import cv2
from utils_node import read_stroke



scale = 5

img = np.zeros((128*scale, 128*scale))+255
file = read_stroke("Sample_points.txt")
for i in range(len(file)):
    stroke = file[i]
    for (x, y, z) in stroke:
        # img[round(y), round(x)] = 0
        cv2.circle(img, (round(x)*scale, round(y)*scale), 3, 0, thickness=3)

cv2.imshow("a", img)
cv2.waitKey(0)
