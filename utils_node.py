"""
Copyright (C) 2019 The Borg Lab - All Rights Reserved

Utils file for unicode2chebyshev
"""
import argparse
import os
import sys
import math

import cairosvg
import cv2
import numpy as np
import svgpathtools
from svgpathtools import svg2paths2, wsvg
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def generate_sample_index(N=1000):
    # Generate index by sampling according to the same PDF
    # sample_even = np.random.rand(N)*2-1 # Range [-1, 1]
    # sample_uneven = np.zeros(N)
    # for i in range(N):
    #     sample_uneven[i] = np.sin(np.pi * sample_even[i] - np.pi/2)
    #
    # sample_uneven = np.sort(sample_uneven)

    cheb2=1
    # Generate index including the beginning and ending points
    if (cheb2):
        sample_uneven = np.zeros(N)  # Range [-1, 1]
        for i in range(N-1):
            sample_uneven[i] = np.cos((N-i-1)/(N-1)*np.pi)
        sample_uneven[N-1]=1
        return sample_uneven
    else:
    # Generate index excluding the beginning and ending points
        sample_uneven = np.zeros(N+1)  # Range [-1, 1]
        for i in range(N+1):
            sample_uneven[i] = np.cos((N-i)/(N+1)*np.pi)
        return sample_uneven[:N]

    # Generate index evenlynum_points
    # return np.linspace(-1,1,N)
def find_index(i,sample_index, stroke_len_in_index):
    for stroke_index in range(int(2/stroke_len_in_index)):
        if sample_index[i] > (-1 + stroke_len_in_index * stroke_index ) and sample_index[i] < (-1 + stroke_len_in_index * (stroke_index+1) ):
            return stroke_index
    print(bcolors.FAIL+"Input parameter for find_index is wrong, please recheck!"+bcolors.ENDC)
    return 0

def sample(stroke, sample_index, expand_part=0.5):
    """
    Given a stroke represented by curves, sample 100 points from the stroke, where t ranges from 0 to 1.

    stroke: list of Bezier curves for one stroke
    num_points: the number of points to be sampled from one stroke
    """
    if len(stroke) < 1:
        print("Input Error of the Sample Function!")

    result_x = []
    result_y = []

    stroke_len_in_big_time = (1-(-1))/len(stroke)
    if(sample_index[-1]!=1):
        for i in range(len(sample_index)):
            t = (sample_index[i] - (-1)) % stroke_len_in_big_time
            stroke_index= int((sample_index[i] - (-1) - t)/ stroke_len_in_big_time)
            t = (sample_index[i] - (-1)) % stroke_len_in_big_time / stroke_len_in_big_time

            curve = stroke[stroke_index]
            if not (t>=0 and t<=1):
                print(1)
            result_x.append(curve.point(t).real)
            result_y.append(curve.point(t).imag)
    else:
        for i in range(len(sample_index[:-1])):
            t = (sample_index[i] - (-1)) % stroke_len_in_big_time
            stroke_index= int((sample_index[i] - (-1) - t)/ stroke_len_in_big_time)
            t = (sample_index[i] - (-1)) % stroke_len_in_big_time / stroke_len_in_big_time


            curve = stroke[stroke_index]
            if not (t>=0 and t<=1):
                print(1)
            result_x.append(curve.point(t).real)
            result_y.append(curve.point(t).imag)
        result_x.append(stroke[len(stroke)-1].point(1).real)
        result_y.append(stroke[len(stroke)-1].point(1).imag)

    result_x = np.array(result_x)
    # Change from SVG axis to Matrix Image representation axis
    result_y = 1024-128-(np.array(result_y))

    return result_x, result_y


def clear(unicode):
    """
    Clears Generated_data_file of all files generated during calculation according to unicode

    unicode: used to find associated files to clear
    """
    files = os.listdir("Generated_data_file")
    for file in files:
        if file[:4] == unicode:
            os.remove("Generated_data_file/"+file)


def clearsvg():
    """
    Clears Generated_data_file of all SVG files generated during calculation
    """
    files = os.listdir("Generated_data_file")
    for file in files:
        if file[-4:-1] == ".sv":
            os.remove("Generated_data_file/"+file)
    print("svg files" + " are removed!")


def resizePNG(outputsize=(128, 128)):
    """
    Resizes output images to 128 x 128

    outputsize: desired outputsize
    """
    files = os.listdir("Generated_data_file")
    for file in files:
        if file[-4:-1] == ".pn":
            file = "Generated_data_file/"+file
            imgtemp = cv2.imread(file, 0)
            # if(imgtemp.shape!=(128,128)):
            #     res = cv2.resize(imgtemp, outputsize,
            #                  interpolation=cv2.INTER_AREA)
            res = cv2.resize(imgtemp, outputsize,
                         interpolation=cv2.INTER_CUBIC)
            # else:
            #     res=cv2.GaussianBlur(imgtemp,(3,3),0)
            cv2.imwrite(file, res)
    print("DownSample finished!")

def scalePNG_pixels():
    files = os.listdir("Generated_data_file")
    for file in files:
        if file[-4:-1] == ".pn":
            file = "Generated_data_file/"+file
            imgtemp = cv2.imread(file, 0)
            minV=np.min(imgtemp)
            if(minV < 44):
                minV=np.round((minV+44)/2.0)
            maxV=255
            # use linear method to sacle the pixel range to 0~255
            aa=(imgtemp-minV)/(maxV-minV)*255
            cv2.imwrite(file, np.uint8(aa))
    print("DownSample finished!")


def convertBW_PNG(unicode):
    """
    Converts images to Black and White

    unicode: unicode of Chinese character
    """
    files = os.listdir("Generated_data_file")
    for file in files:
        if file[-4:-1] == ".pn" and file[:4] == unicode:
            file = "Generated_data_file/"+file
            imgtemp = cv2.imread(file, 0)
            cv2.imwrite(file, 255-imgtemp)
    print("Convert Black & White finished!")

def merge_strokes(unicode):
    """
    merge stroke images into one character image
    """
    files = os.listdir("Generated_data_file")
    mat=np.zeros((128,128))+255
    for file in files:
        if file[-4:-1] == ".pn" and file[:4] == unicode:
            file = "Generated_data_file/"+file
            imgtemp = cv2.imread(file, 0)
            for i in range(128):
                for j in range(128):
                    if(mat[i,j]>imgtemp[i,j]):                    
                        mat[i,j] = imgtemp[i,j]
    cv2.imwrite("Generated_data_file/"+"ori.png", mat)
    #cv2.imwrite("Generated_data_file/"+unicode +"_ori.png", mat)
    print("Merging finished!")


def separate_strokes(unicode):
    """
    Takes in a unicode and extracts individual strokes and stores them as separate images.

    unicode: unicode of Chinese character
    """
    # Find skeleton SVG file for character "unicode" from hanzivg folder.
    original_svg_filename = "../Data/svgs/" + \
        str(int(unicode, 16)) + ".svg"
    cairosvg.svg2png(url=original_svg_filename,
                     write_to="Generated_data_file/"+"ori.png")

    # Extract all the individual strokes in the original file
    paths, attributes, svg_attributes = svg2paths2(original_svg_filename)
    svg_attributes['style'] = "fill:none;stroke:#000000;stroke-width:3;" \
        "stroke-linecap:round;stroke-linejoin:round;"
    for att in attributes:  # Fix the "upside-down" problem
        att['transform'] = "scale(1, -1) translate(0, -900)"

    # Distinguish between skeleton and normal strokes
    skeleton_index = np.linspace(
        int((len(paths) - 4) / 3 + 1), int(len(paths) - 4 + 1 - 2), int((len(paths) - 4) / 3))
    stroke_index = np.linspace(
        0, int((len(paths) - 4) / 3) - 1, int((len(paths) - 4) / 3))
    skeleton_index = list(np.uint(skeleton_index))
    stroke_index = list(np.uint(stroke_index))

    for index in stroke_index:
        filename_temp = "Generated_data_file/"+unicode + \
            '_' + 'stroke' + '{:0>2}'.format(str(index)) + '.svg'
        png_filename = "Generated_data_file/"+unicode + \
            '_' + 'stroke' + '{:0>2}'.format(str(index)) + '.png'
        wsvg(paths[index], attributes=attributes,
             svg_attributes=svg_attributes, filename=filename_temp)
        # Convert svg to png
        cairosvg.svg2png(url=filename_temp, write_to=png_filename, output_width=128, output_height=128)
        # cairosvg.svg

    return skeleton_index, paths

def show_fitting_result(x_samples, x_coeff, index_big, if_show=False, if_save=False):
    if not if_show:
        return None
    else:
        size0=len(x_samples)
        test_index=np.random.uniform(-1,1,size0)
        test_index=sorted(test_index)
        test_res = np.polynomial.chebyshev.chebval(test_index,x_coeff)
        ref_index=range(size0)
        plt.plot(ref_index, test_res,'r--', size0/2*index_big+size0/2,x_samples,'b--')
        if(if_save):
            plt.savefig("uuneven-fitting.png")
        plt.show()


def record_coefficients(unicode, paths, skeleton_index, num_points, deg,  initial_z=3, delta_x=10, adaptive_z=False):
    # change_z=False will generate chebyshev coefficients with fixed z, True will generate one file with adaptive z
    # Open file to store Chebyshev coefficients
    f = open("Generated_data_file/node_cheby_sample.txt", 'w')
    idx = 0
    # Going through all paths in svg and only analyzing the paths that feature skeleton of stroke

    for index in skeleton_index:  # Only analyze for the skeleton paths

        cur_stroke = paths[index]

        png_filename_temp = "Generated_data_file/" + \
                            unicode + '_' + 'stroke' + \
                            '{:0>2}'.format(str(idx)) + '.png'

        original_png_stroke = cv2.imread(png_filename_temp, 0)

        # Sample discrete points for all the component curves of one stroke
        #index_big = generate_sample_index(len(cur_stroke) * 200)
        index_big = generate_sample_index(deg)

        x_samples, y_samples = sample(cur_stroke, index_big)
        z_samples = np.zeros(len(x_samples))+initial_z*100

        #x_samples=x_samples+delta_x

        poly = lagrange(index_big, x_samples)
        x_coeff =Polynomial(poly).coef
        test_index=np.random.uniform(-1,1,200)
        test_index=np.sort(test_index)
        plt.plot(test_index,np.polyval(x_coeff,test_index))
        plt.plot(index_big,x_samples)
        # plt.show()



        # Transform sampled points to Chebyshev coefficients of given degree
        # t_for_cheby = np.linspace(0, 1, len(x_samples))

        np.savetxt(f, np.expand_dims(x_samples, 0))
        np.savetxt(f, np.expand_dims(y_samples, 0))
        np.savetxt(f, np.expand_dims(z_samples, 0))
        # Save the extracted coefficients of each stroke in the file
        idx = idx + 1
    f.close()
    print("The initial direction: ", math.atan2(cur_stroke[0][0].real-cur_stroke[0][1].real,cur_stroke[0][0].imag-cur_stroke[0][1].imag)*180/3.14)


def record_init_trajectory(unicode, paths, skeleton_index, num_points, deg, initial_z=3, adaptive_z=False):
    # change_z=False will generate chebyshev coefficients with fixed z, True will generate one file with adaptive z
    # Open file to store Chebyshev coefficients
    f = open("Generated_data_file/coeff_cheby_sample.txt", 'w')
    idx = 0
    # Going through all paths in svg and only analyzing the paths that feature skeleton of stroke
    color_trajectory=np.zeros([1024,1024,3])+255


    for index in skeleton_index:  # Only analyze for the skeleton paths

        cur_stroke = paths[index]

        # Sample discrete points for all the component curves of one stroke
        index_big = generate_sample_index(len(cur_stroke) * 200)

        x_samples, y_samples = sample(cur_stroke, index_big)
        for i in range(len(x_samples)):
            cv2.circle(color_trajectory, (np.int(y_samples[i]), np.int(x_samples[i])), 5, (0,0,255))
        idx = idx + 1
    f.close()
    # cv2.imshow("a",color_trajectory)
    # cv2.waitKey(0)
    return color_trajectory

def visualize_dataset(unicode, paths, skeleton_index, num_points, deg, adaptive_z=False):
    # change_z=False will generate chebyshev coefficients with fixed z, True will generate one file with adaptive z
    # Open file to store Chebyshev coefficients
    f = open("Generated_data_file/coeff_cheby_sample.txt", 'w')
    idx = 0
    # Going through all paths in svg and only analyzing the paths that feature skeleton of stroke
    samples_big_x=[]
    samples_big_y=[]
    for i in range(len(skeleton_index)):
        skeleton_index.append(int(skeleton_index[i]-1))
    for index in skeleton_index:  # Only analyze for the skeleton paths

        cur_stroke = paths[index]

        png_filename_temp = "Generated_data_file/" + \
                            unicode + '_' + 'stroke' + \
                            '{:0>2}'.format(str(idx)) + '.png'

        original_png_stroke = cv2.imread(png_filename_temp, 0)

        # Sample discrete points for all the component curves of one stroke
        index_big = generate_sample_index(len(cur_stroke) * 200)

        x_samples, y_samples = sample(cur_stroke, index_big)

        samples_big_x.append(x_samples)
        samples_big_y.append(y_samples)
        idx = idx + 1
    f.close()
    big_x = np.array(samples_big_x)
    big_y = np.array(samples_big_y)
    skeleton_img=np.zeros((1024,1024))+255
    for i in range(len(big_x)):
        for j in range(len(big_x[i])):
            cv2.circle(skeleton_img, (np.int(big_y[i][j]), np.int(big_x[i][j])), 3, 0, -1)
            skeleton_img[np.int(big_x[i][j]), np.int(big_y[i][j])]=0
    cv2.imshow("a",skeleton_img)
    cv2.waitKey(0)
    cv2.imwrite("initial_Niao.jpg")
    print("The initial direction: ", math.atan2(cur_stroke[0][0].real-cur_stroke[0][1].real,cur_stroke[0][0].imag-cur_stroke[0][1].imag)*180/3.14)


def merge_pic(img1,img2):
    shape=img1.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                if (img2[i,j,k]<img1[i,j,k]):
                    img1[i, j, k]=img2[i,j,k]
    # cv2.imshow("a", img1)
    # cv2.waitKey(0)
    cv2.imwrite("merge.png",img1)


def read_stroke(fileName, num_points=50):
    '''
    Read txt file to generate a list of strokes

    Args:
        fileName(str): The name of the txt file, including .txt

    Returns:
        list_of_strokes: The list of strokes of the input character
        size: len(stroke) * num_points * 3
    '''

    with open('Generated_data_file/'+fileName) as f:
        lines = f.read().splitlines()

    index=0
    coordinates = []
    for line in lines:
        line = line.strip("\'[]").split(',')
        num = float(line[0])
        if(index%3==0):
            point3 = [num]
        else:
            point3.append(num)
        if(len(point3)==3):
            coordinates.append(point3)
        index = index +1

    num_of_strokes = int(len(coordinates)/num_points)
    stroke = []
    list_of_strokes = []

    for i in range(0, num_of_strokes):
        for j in range(0, num_points):
            stroke.append(coordinates[i*num_points+j])
        list_of_strokes.append(stroke)
        stroke = []
    return list_of_strokes
