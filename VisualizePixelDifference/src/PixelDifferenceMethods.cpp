#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "PixelDifferenceMethods.h"

using namespace ::std;
typedef Eigen::Matrix<double, 128, 128> Canvas;

const int processImages(std::string originalImagePath, std::string resultImagePath, std::string resultFilename, int dimensions)
{
    // Read in images, resizes to 640 x 640
    Eigen::MatrixXd originalMatrix = getImage(originalImagePath, dimensions);
    Eigen::MatrixXd resultMatrix = getImage(resultImagePath, dimensions);

    // Get the value of differences between original and result
    Eigen::MatrixXd differencesMatrix = subtractMatrices(originalMatrix, resultMatrix);

    // !!!!!!!!!!!!!!!!!!!!!!!!!
    blendMatrices(originalMatrix, resultMatrix, resultFilename);
    // Only take differences above a certain threshold
    return calculatePixelDifference(differencesMatrix);
}

const int calculatePixelDifference(Eigen::MatrixXd diffMatrix)
{
    int pixDiff = 0;
    for (size_t i = 0; i < diffMatrix.rows(); i++)
    {
        for (size_t j = 0; j < diffMatrix.cols(); j++)
        {
            if (diffMatrix(i, j) != 0)
            {
                pixDiff++;
            }
        }
    }
    return pixDiff;
}

/**
 * Gets an image and also resizes it
 */
const Eigen::MatrixXd getImage(std::string path, int dimensions)
{
    cv::Mat image;
    image = cv::imread(path, cv::IMREAD_GRAYSCALE);
    if (image.size == 0)
    {
        cout << "File " << path << " not read" << endl;
    }
    cv::resize(image, image, cv::Size(dimensions, dimensions), cv::INTER_NEAREST);

    // cv::bitwise_not(image, image);
    // cout << "bitwise not" << endl;

    Eigen::MatrixXd eigen_mat;
    cv::cv2eigen(image, eigen_mat);

    // If pixel is above threshold of black, make it all the way black, otherwise white
    for (size_t i = 0; i < eigen_mat.rows(); i++)
    {
        for (size_t j = 0; j < eigen_mat.cols(); j++)
        {
            if (eigen_mat(i, j) < 200)
            {
                eigen_mat(i, j) = 0;
            }
            else
            {
                eigen_mat(i, j) = 255;
            }
        }
    }

    return eigen_mat;
}

const Eigen::MatrixXd subtractMatrices(Eigen::MatrixXd original, Eigen::MatrixXd result)
{
    return original - result;
}

void blendMatrices(Eigen::MatrixXd original, Eigen::MatrixXd result, std::string resultFilename)
{
    Eigen::MatrixXd differencesMatrix = subtractMatrices(original, result);
    for (size_t i = 0; i < differencesMatrix.rows(); i++)
    {
        for (size_t j = 0; j < differencesMatrix.cols(); j++)
        {
            if (differencesMatrix(i, j) == 0)
            {
                // Final will be white
                original(i, j) = 255;
                result(i, j) = 255;
                differencesMatrix(i, j) = 255;
            }
            else if (differencesMatrix(i, j) < 0)
            {
                original(i, j) = 255;
                result(i, j) = 0;
                differencesMatrix(i, j) = 0;
            }
            else
            {
                original(i, j) = 0;
                result(i, j) = 255;
                differencesMatrix(i, j) = 0;
            }
        }
    }

    cv::Mat originalImage;
    cv::eigen2cv(original, originalImage);
    cv::Mat resultImage;
    cv::eigen2cv(result, resultImage);
    cv::Mat diff;
    cv::eigen2cv(differencesMatrix, diff);

    vector<cv::Mat> channels;
    channels.push_back(diff);
    channels.push_back(resultImage);
    channels.push_back(originalImage);

    cv::Mat fin_img(originalImage.rows, originalImage.cols, originalImage.type());
    cv::merge(channels, fin_img);
    cv::imwrite("../ImageOverlaying/blended_SI_scanner/" + resultFilename, fin_img);
}