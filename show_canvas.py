""" Show a file created by C++ file_io::WriteCanvas
"""
import argparse

import cv2
import numpy as np

"""
picture show
"""


def run(path):
    # Open file
    f = open(path, 'r')

    # Create numpy array
    mattemp = []
    for line in f:
        linetemp = line.split(' ')
        chatemp = []
        for cha in linetemp:
            if(cha != '' and cha != "\n"):
                chatemp.append(float(cha))
        mattemp.append(chatemp)
    mattemp = np.array(mattemp)
    mattemp = np.uint8(mattemp)

    # Can close file now
    f.close()

    # resize and show
    scale = 5
    mattemp = cv2.resize(mattemp, (128*scale, 128*scale))
    cv2.imshow(path+" Picture by Chebyshev", mattemp)
    cv2.waitKey(0)


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    args = parser.parse_args()
    run(args.path)
