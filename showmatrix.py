"""
Copyright (C) 2019 The Borg Lab - All Rights Reserved

Python file for visualizing the image output from GTSAM optimization
"""
import argparse
import sys

import show_canvas

parser = argparse.ArgumentParser()
parser.add_argument(
    '--file_name', default="res", type=str, metavar="file_name", help="the name of the image to show, could be res, ori, skeleton")

if __name__ == "__main__":
    """
    Main function that takes in a file name and displays an image
    """
    if len(sys.argv) > 3:
        print("The argument is not recognized!")

    args = parser.parse_args()
    file_name = args.file_name

    if file_name == "ori":
        img = cv2.imread("Generated_data_file/ori.png")
        cv2.imshow("Original image", img)
        cv2.waitKey(0)
    elif (file_name != "res" and file_name != "skeleton"):
        print("The argument is not recognized!")
    else:
        path = 'Generated_data_file/'+file_name+"_picture.txt"
        show_canvas.run(path)
