from PIL import Image
import cv2

from os import listdir
from os.path import isfile, join
onlyfiles = [join("../Data/written_results_scanned_using_scanner", f) for f in listdir("../Data/written_results_scanned_using_scanner") if isfile(join("../Data/written_results_scanned_using_scanner", f))]

dim = (640, 640)

for f in onlyfiles:
    if f.endswith(".png"):
        print(f)
        resized = cv2.imread(f)
        resized = cv2.resize(resized, dim, interpolation=cv2.INTER_NEAREST)
        cv2.imwrite(f, resized)