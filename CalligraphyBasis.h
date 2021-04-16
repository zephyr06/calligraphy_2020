
/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * Header file that includes the basic functions that do not depend on any type of brush or state
 */
#pragma once

#include <CppUnitLite/TestHarness.h>
#include <dirent.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/inference/Key.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
// #include <gtsam/3rdparty/Eigen/Eigen/E
// #include <eigen/
// #include <gtsam/3rdparty/
#include <gtsam/nonlinear/Values.h>
#include <math.h>
#include <sys/types.h>
#include <opencv2/highgui.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <tuple>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>

#include "Chebyshev2.h"
#include "VirtualBrush.h"
#include "file_io.h"
#include "Parameters.h"

using namespace std;
using namespace gtsam;
using namespace cv;
using symbol_shorthand::A;

#define PI 3.14159265
// const int NUM_SAMPLES = 200;
typedef Eigen::Matrix<double, 1, ORDER> Basis_another;

// Coefficients are the representation of one stroke in the pseudo-spectral way, i.e.,
// the values of the trajectory at Chyebshev points are stored there
typedef Eigen::Matrix<double, ORDER * 3, 1> Coefficients;
typedef Eigen::Matrix<double, 1, kImageHeight * kImageWidth> PictError;

void ReadBrushParameters()
{
}

/**
 * @brief Create Chebyshev points in [-1,1]
 *
 * @param numPoints number of Chebyshev points
 * @return vector<double> vector of Chebyshev points
 */
vector<double> GenerateChebyshevPoints(int numPoints);

/**
 * @brief Create uniform points in [-1,1]
 *
 * @param numPoints number of uniform points to sample
 * @return vector<double> vector of uniformly sampled points
 */
vector<double> GenerateUniform(int numPoints);

/**
 * @brief Sample 3D points using Chebyshev coefficients
 *
 * @param coefficients vector of length ORDER*3 where the first third of the
 * vector is the Chebyshev coefficients for the x's, second third for y's,
 * and last third is for z's
 * @param uniform boolean flag for uniform or Chebyshev sampling
 * @return vector<Point3> vector of Point3 samples
 */
vector<Point3> SamplePoints(Coefficients coefficients, int sampleNum = NUM_SAMPLES, bool uniform = true);

/**
 * @brief Calculates the pixel differences between the original and generated
 * pictures
 * 
 * @param generatedPicture Canvas representing the simulated picture
 * @param targetPicture Canvas representing goal picture
 */
PictError PictureDifference(Canvas generatedPicture, Canvas targetPicture);

/**
 * Used by testAnalyzeCharacter.cpp to read the initial Chebyshev
 * coefficients for optimization
 *
 * @param path: path from which to read Chebyshev coefficients
 * @param ratio: transform our images by this ratio, from 1024 to 128
 * @return vector<Coefficients>: vector of coefficients of the x, y, z
 * Chebyshev polynomials
 */
vector<Coefficients> ReadCoefficients(string path, double ratio = 8);

/**
 * Saving generated sample points (x, y, z) into file
 * "Sample_points.txt"
 *
 * @param resultCoefficients: optimized Chebyshev polynomials
 * @param unicode: unicode of Chinese character
 */
void SaveSamplePoints(vector<Coefficients> resultCoefficients,
                      string unicode);

/**
 * Samples a set of points to use for drawing given the characterCoefficients
 *
 * @param characterCoefficients: Coefficients of current character
 * @return vector<vector<Point3>>: set of points per stroke
 */
vector<vector<Point3>> SampleCharacterPoints(
    const vector<Coefficients> &characterCoefficients, int localNUM_SAMPLES_Draw = NUM_SAMPLES);
