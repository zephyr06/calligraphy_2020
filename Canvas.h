/**
 * @file Canvas.h
 * @brief Definition and initialization of Canvas class
 * @author Frank Dellaert
 * @date September 2019
 */

#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <dirent.h>
#include <math.h>
#include <sys/types.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

const int kImageHeight = 128, kImageWidth = 128;
typedef Eigen::Matrix<double, kImageHeight, kImageWidth> Canvas;

/**
 * @brief Initializes an empty white canvas
 * 
 * @return Canvas of size kImageHeight x kImageWidth as specified in Canvas.h
 */
Canvas InitCanvas();

/**
 * @brief Merge two images, only darkening pixel values
 *
 * @param imageToAdd
 * @param inOutImage
 */
void MergeCanvas(
    const Canvas &imageToAdd,
    Canvas *inOutImage);
