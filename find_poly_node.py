"""
Copyright (C) 2019 The Borg Lab - All Rights Reserved

This file takes in Unicode representations and will generate the necessary data files for GTSAM to optimize.
Calligraphy strokes have fixed width.
"""
import argparse

# from utils import (clear, clearsvg, convertBW_PNG, cv2, find_radius, np,
#                    resizePNG, sample, separate_strokes, record_coefficients, merge_strokes)
from utils_node import *

parser = argparse.ArgumentParser()
parser.add_argument(
    '--unicodes',  default="4e7f", type=str, metavar="unicodes", help="The unicodes of characters to analyze.")
parser.add_argument(
    '--adaptive_z',  default=False,   action="store_true", help="The unicodes of characters to analyze.")
parser.add_argument(
    '--deg',  default=6, type=int, help="The order of characters to analyze.")


if __name__ == "__main__":
    """
    Takes in a unicode of a character and generates the Chebyshev coefficients for the
    x, y, and z coordinates we need to draw.
    """
    # Takes in command line arguments, unicode in hex separated by spaces, one per character.
    args = parser.parse_args()

    num_points = 50  # Number of points to sample per stroke
    # Degree of Chebyshev polynomial that we will generate
    deg = int(args.deg)
    clear = False

    unicode = args.unicodes

    skeleton_index, paths = separate_strokes(unicode)
    convertBW_PNG(unicode)
    record_coefficients(unicode, paths, skeleton_index,
                        num_points, deg, initial_z=1.8*2.5, delta_x=10, adaptive_z=args.adaptive_z)


clearsvg()
# resizePNG()
scalePNG_pixels()
merge_strokes(unicode)
if clear:
    clear(unicode)
# show_initial_direction()
