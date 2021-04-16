import cv2

img=cv2.imread("ori.png",0)
cv2.imwrite("ori.png", cv2.resize(img, (128*5,128*5), interpolation = cv2.INTER_CUBIC))

