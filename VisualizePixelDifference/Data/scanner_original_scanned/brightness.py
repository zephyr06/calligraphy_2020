import cv2
import numpy as np

img=cv2.imread("Kong-opt-dynamic.png",0)
minV=np.min(img)
maxV=np.max(img)
print(minV, maxV)
minV2=40
