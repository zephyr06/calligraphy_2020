/*
 * @file SimpleVirtualBrush.h
 * @brief File for the simple virtual brush and the brush state for simple brush
 * 
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 */

#pragma once
#include "VirtualBrush.h"

// TOOD: move these to cpp
#include <opencv2/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/**
 * @brief Simple Virtual Brush Control
 */
class SimpleBrushControl
{
public:
    double x;
    double y;
    double z;
    SimpleBrushControl(const Point3 &point) : x(point.x()), y(point.y()), z(point.z()) {}
};

/**
 * @brief Simple Virtual Brush State
 */
class SimpleBrushState
{
public:
    double x;
    double y;
    double r;
    double kScaleShape;

    /**
     * Constructor
     * 
     * @param kScaleInput the scale of the Canvas
     */
    SimpleBrushState(double kScaleInput = 1)
    {
        x = 0;
        y = 0;
        r = 0;
        kScaleShape = kScaleInput;
    }

    /**
     * Constructor
     * 
     * @param sth points to sample
     * @param kScaleInput
     */
    SimpleBrushState(vector<Point3> sth, double kScaleInput = 1)
    {
        x = 0;
        y = 0;
        r = 0;
        kScaleShape = kScaleInput;
    }

    /**
     * Constructor
     * 
     * @param sth Coefficients of Chebyshev polynomials
     * @param kScaleInput
     */
    SimpleBrushState(Coefficients sth, double kScaleInput = 1)
    {
        x = 0;
        y = 0;
        r = 0;
        kScaleShape = kScaleInput;
    }

    /**
     * @brief Gets the width of the mark
     * 
     * @param z the height of the brush location
     * @return double width of mark
     */
    double Width(double z) // The length of the tip that have contact with the paper
    {
        double res = 0;

        vector<double> coeff{1.19034619e+00, -3.20286161e-04};
        if (z <= 0)
            return 0;
        else
            res = kScaleReal2Image * (coeff[0] * z + coeff[1]) * kScaleShape;
        if (res < 0.01)
        {
            return 0;
        }
        return res;
    }

    /**
     * @brief Updates the parameters of the brush given a control
     * 
     * @param &control the next Point3 coordinate of the brush location
     */
    void Update(const Point3 &control)
    {
        x = control.x();
        y = control.y();
        r = Width(control.z());
    }

    void ShowState()
    {
        cout << x << ", " << y << ", " << r << endl;
    }
};

class SimpleVirtualBrush : public VirtualBrush<SimpleBrushState>
{
public:
    using VirtualBrush::VirtualBrush;

    using HiResImage = Eigen::Matrix<double, 100, 100>;

    /**
     * @brief Draw a mark in a high-resolution image
     *
     * @return HiResImage high resolution image with drawn mark
     */
    HiResImage DrawMarkInHiRes() const
    {
        HiResImage hiResImage = HiResImage::Ones() * 255;

        // picture coordinate: left-down corner is the origin, and x-axis is
        // just the horizontal one The input point uses the picture coordinate
        static const double kScale = 5;
        double R = kScale * currentState_.r / 2.0;

        Point2 imageCenter(50, 50);
        // Point2 center(imageCenter.x() + offset_len, imageCenter.y());

        for (int v = round(imageCenter.y() - R); v < round(imageCenter.y() + R);
             v++)
        {
            int temp = round(
                sqrt(R * R - (v - imageCenter.y()) * (v - imageCenter.y())));
            for (int u = imageCenter.x() - temp; u < imageCenter.x() + temp;
                 u++) // floor rounding used, actually
            {
                if (v >= 0 && v <= kImageWidth && u >= 0 && u <= kImageHeight)
                    hiResImage(v, u) = 0;
            }
        }

        return hiResImage;
    }

    /**
     * @brief Move blurred image of mark to the correct place, including calculation of
     * direction and offset
     *
     * @param blurredImage OpenCv image of size 20x20
     * @return Canvas blurred image
     */
    Canvas MoveMarkToRightPlace(const cv::Mat &blurredImage) const
    {
        // define arguments

        // Moving part
        cv::Mat matWithMark;
        cv::Mat_<double> affineTransform(2, 3);
        affineTransform << 1, 0,
            currentState_.x - 10, //
            0, 1,
            currentState_.y - 10;
        cv::Size dsize(kImageHeight, kImageWidth);
        cv::Scalar borderValue = cv::Scalar(255, 255, 255);

        cv::warpAffine(blurredImage, matWithMark, affineTransform, dsize,
                       cv::INTER_LINEAR, cv::BORDER_CONSTANT, borderValue);

        // Convert and return
        Canvas imageWithMark;
        cv::cv2eigen(matWithMark, imageWithMark);
        return imageWithMark;
    }

    /**
     * @brief New enhanced mark drawing
     *
     * @return Canvas drawing with enhanced mark
     */
    Canvas DrawMark() const override
    {
        if (currentState_.r <= 0)
            // return one white picture when the brush lifts upon the paper
            return InitCanvas();

        // 1a. Draw mark in high-res image
        auto highResImage = DrawMarkInHiRes();

        // 1b. Scale down image to get automatic anti-aliasing
        cv::Mat openCvImage;
        cv::eigen2cv(highResImage, openCvImage);
        cv::Mat antiAliasedImage;
        cv::Size dsize(20, 20);
        double fx = 0.2, fy = 0.2;
        cv::resize(openCvImage, antiAliasedImage, dsize, fx, fy,
                   cv::INTER_AREA);

        // 2. Filter with a Gaussian
        cv::Mat blurredImage = antiAliasedImage;

        // 3. Move mark to the right place
        return MoveMarkToRightPlace(blurredImage);
    }

    /**
     * @brief Draws a single stroke into a dynamic Eigen matrix
     *
     * @param strokeTrajectory sampled points on the stroke trajectory
     * polynomial
     * @return Canvas a dynamic Eigen matrix with the stroke drawn
     */
    Canvas DrawStroke(const vector<Point3> &strokeTrajectory)
    {

        Canvas strokeImage = InitCanvas();
        currentState_ = lastState_;

        // For points after the transition point
        for (int i = 0;
             i < static_cast<int>(strokeTrajectory.size()); i++)
        {
            currentState_.Update(strokeTrajectory[i]);
            const Canvas markImage = DrawMark();
            MergeCanvas(markImage, &strokeImage);
        }

        currentState_.Update(strokeTrajectory.back());

        return strokeImage;
    }
    Canvas DrawStroke(const Coefficients &coeff)
    {
        vector<Point3> points = SamplePoints(coeff);
        return DrawStroke(points);
    }
};
