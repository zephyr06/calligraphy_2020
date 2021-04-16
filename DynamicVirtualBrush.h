

/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * The derived class for the dynamic virtual brush class
 */
#pragma once
#include "DynamicBrushState.h"
#include "VirtualBrush.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/**
 * @brief The dynamic virtual brush derived from virtual brush class
 */
class EnhancedDynamicVirtualBrush : public VirtualBrush<DynamicBrushState>
{
public:
    /**
     * Constructor for EnhancedDynamicVirtualBrush
     * 
     * @param &initialState initial state of brush
     * @param character_size the scale of our character
     **/
    EnhancedDynamicVirtualBrush(const DynamicBrushState &initialState,
                                double character_size = kCharacterSize) //
        : VirtualBrush(initialState, character_size)
    {
    }

    const double kScale = 5;
    const int markLen = 20;
    // Maybe there is a way to solve this problem, the size of HiResImage sets the limits for anti-alising problem
    using HiResImage = Eigen::Matrix<double, 100, 100>;

    /**
     * @brief Draw a mark in a high-resolution image, the drawing parameters are gotten from currentState_
     *
     * @param Null
     * @return HiResImage
     */
    HiResImage DrawMarkInHiRes() const
    {
        HiResImage hiResImage = HiResImage::Ones() * 255;

        // picture coordinate: left-down corner is the origin, and x-axis is
        // just the horizontal one The input point uses the picture coordinate

        double drag_len = kScale * currentState_.drag;
        // double offset_len = kScale * currentState_.offset;
        double width = kScale * currentState_.width;

        // Code in rotation part is not elegent, we need to transform the
        // current mark to the center of the picture, and then rotate, and then
        // transform it back;

        Point2 imageCenter(hiResImage.rows() / 2, hiResImage.cols() / 2);
        // Point2 center(imageCenter.x() + offset_len, imageCenter.y());

        VBCrossSection curve_profile(
            0); // use quadratic curve for the profile of the tip mark
        curve_profile.Fitting(width / 2.0, drag_len);

        for (int v = round(imageCenter.y() - width / 2.0);
             v < round(imageCenter.y() + width / 2.0); v++)
        {
            int temp = round(curve_profile.Calculate(imageCenter.y() - v));
            for (int u = imageCenter.x(); u < imageCenter.x() + temp;
                 u++) // floor rounding used, actually
            {
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
     * @return Canvas
     */
    Canvas moveMarkToRightPlace(const cv::Mat &blurredImage) const
    {
        // define arguments

        // Rotation part
        cv::Point2f center_tran(blurredImage.cols / 2, blurredImage.rows / 2);

        cv::Mat mat_rotate = cv::getRotationMatrix2D(
            center_tran, currentState_.currentDirection, 1);
        cv::Scalar borderColor = cv::Scalar(255, 255, 255);

        cv::warpAffine(blurredImage, blurredImage, mat_rotate,
                       blurredImage.size(), cv::INTER_LINEAR,
                       cv::BORDER_CONSTANT, borderColor);

        // Moving part
        cv::Mat matWithMark;
        cv::Mat_<double> affineTransform(2, 3);
        affineTransform << 1, 0,
            currentState_.x - 10 +
                currentState_.offset *
                    cos(currentState_.currentDirection / 180.0 * CV_PI), //
            0, 1,
            currentState_.y - 10 -
                currentState_.offset *
                    sin(currentState_.currentDirection / 180.0 * CV_PI);
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
     * @brief The function draws a written mark from the current state of the brush
     * 
     * @return Canvas with a drawn written mark using current brush state
     */
    Canvas DrawMark() const override
    {
        if (currentState_.z <= 0)
            // return one white picture when the brush lifts upon the paper
            return InitCanvas();

        // 1a. Draw mark in high-res image
        auto highResImage = DrawMarkInHiRes();

        // 1b. Scale down image to get automatic anti-aliasing
        cv::Mat openCvImage;
        cv::eigen2cv(highResImage, openCvImage);
        cv::Mat antiAliasedImage;
        cv::Size dsize(markLen, markLen);
        double fx = 1.0 / kScale;
        double fy = fx;
        cv::resize(openCvImage, antiAliasedImage, dsize, fx, fy,
                   cv::INTER_AREA);

        // 2. Filter with a Gaussian
        // cv::Mat blurredImage;
        // cv::GaussianBlur(antiAliasedImage, blurredImage, cv::Size(5, 5), 0, 0);

        // 3. Move mark to the right place
        return moveMarkToRightPlace(antiAliasedImage);
    }

    Canvas DrawStroke(const Coefficients &coeff)
    {
        vector<Point3> points = SamplePoints(coeff);
        return DrawStroke(points);
    }
    /**
     * @brief The function draws a stroke given its representation
     *
     * @param strokeTrajectory, the stroke trajectory in vector<Point3> form
     * @return Canvas, the stroke image
     */
    Canvas DrawStroke(const vector<Point3> &strokeTrajectory) override
    {

        Canvas strokeImage = InitCanvas();
        currentState_ = lastState_;
        currentState_.currentDirection = JudgeInitialDirection(strokeTrajectory);

        // warning, changes the `strokeTrajectory` instance variable !!!
        // int strokeTransitionIndex_ = FindTransitionPoint(strokeTrajectory);
        int strokeTransitionIndex_ = 0;

        // For points before the transition point
        for (int i = 0; i < strokeTransitionIndex_; i++)
        {
            currentState_.Update(strokeTrajectory[i], "between");
            const Canvas markImage = DrawMark();
            MergeCanvas(markImage, &strokeImage);
        }
        // For points after the transition point
        for (int i = strokeTransitionIndex_;
             i < static_cast<int>(strokeTrajectory.size()); i++)
        {
            currentState_.Update(strokeTrajectory[i], "in");

            const Canvas markImage = DrawMark();
            MergeCanvas(markImage, &strokeImage);
        }

        currentState_.Update(strokeTrajectory.back(), "last");
        lastState_ = currentState_;

        return strokeImage;
    }

    /**
     * @brief The function finds the transition point in a stroke. 
     * The parameters of the brush changes in a different way before and after the transition point:
     * Before: the inertia is much bigger, and it also changes in a slightly different way
     * After: the inertia is smaller
     * 
     * @param stroke: the stroke to find the transition point
     * @return the index of the transition point
     **/
    int FindTransitionPoint(vector<Point3> stroke)
    {

        double maxDistance = lastState_.offset;

        double dist_sum = 0;

        for (int i = 1; i < static_cast<int>(stroke.size()); i++)
        {

            dist_sum += sqrt((stroke[i][0] - stroke[i - 1][0]) *
                                 (stroke[i][0] - stroke[i - 1][0]) +
                             (stroke[i][1] - stroke[i - 1][1]) *
                                 (stroke[i][1] - stroke[i - 1][1]));

            if (dist_sum > maxDistance)
            {
                return i;
            }
        }
        cout << "Warning: Didn't find the transition point in the stroke, "
                "maybe you need to adjust parameters for FindTransitionPoint"
             << ", maxDistance: " << maxDistance << endl;
        return 1;
    }
};
