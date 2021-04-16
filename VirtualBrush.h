/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * The base class for the virtual brush
 */

#pragma once

#include "Canvas.h"
#include "file_io.h"

#include <gtsam/base/numericalDerivative.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/inference/Key.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>

#include <CppUnitLite/TestHarness.h>

#include "Parameters.h"

using namespace std;
using namespace gtsam;

/**
 * @brief Calculate the sigmoid function given an input
 * @param x, the input value
 * @output the output value
 */
double sigmoid(double x) { return 1 / (1 + exp(-x)); }

/**
 * @brief the base class of Virtual Brush model
 *
 * Template argument BRUSH should be a functional class that has
 *   a default constructor and a constructor that takes stroke input (either vector<Point3> or Coefficients)
 *   a 'render method' that renders something on an image
 *   a Update(Point3) method that takes one point and updates the inner parameters of virtual brush
 */
template <class BRUSH_STATE>
class VirtualBrush
{
protected:
    double character_size; // the width of the character you expect to
                           // write in the real-world coordinate

    BRUSH_STATE lastState_;    // record the direction of the tip after writing
                               // the last stroke, used only by CharacterDraw
    BRUSH_STATE currentState_; // Current state of the brush

public:
    /**
     * @brief Default constructor for use in ChuVB
     */
    explicit VirtualBrush() {}

    /**
     * Constructor
     *
     * @param initialState for the virtual brush
     * @param character_size is the size of the simulated character
     */
    VirtualBrush(const BRUSH_STATE &initialState,
                 double character_size = kCharacterSize) //
        : character_size(character_size)
    {
        lastState_ = initialState;
        currentState_ = initialState;
    }

    // The derived classes have to override this method
    virtual Canvas DrawMark() const = 0;
    // The derived classes have to override this method
    virtual Canvas DrawStroke(const vector<Point3> &strokeTrajectory) = 0;

    /**
     * @brief This function draws a character given its trajectories using this virtual brush
     * 
     * @param strokeTrajectories is the trajectories for all the strokes in this character
     * @return a Canvas that contains the simulated character
     **/
    Canvas DrawCharacter(const vector<vector<Point3>> &strokeTrajectories)
    {
        // Canvas characterImage = draw_stroke_first(strokeTrajectories[0]);
        Canvas characterImage = InitCanvas();

        // cout << "Size : " << strokeTrajectories.size() << endl;
        for (int k = 0; k < int(strokeTrajectories.size()); k++)
        {
            Canvas strokeImage = DrawStroke(strokeTrajectories[k]);
            MergeCanvas(strokeImage, &characterImage);
        }
        return characterImage;
    }

    /**
     * @brief This function is mainly used for debug purpose
     * 
     * @param strokeTrajectories the trajectories of stroke that we want to draw
     * @return Canvas drawn with the strokeTrajectories
     **/
    Canvas DrawCharacterWithState(vector<vector<Point3>> strokeTrajectories)
    {
        // Canvas mergedImage = draw_stroke_first(strokeTrajectories[0]);
        Canvas mergedImage = InitCanvas();

        // showCurrentState();
        for (const auto &strokeTrajectory : strokeTrajectories)
        {
            const Canvas strokeImage = DrawStroke(strokeTrajectory);
            MergeCanvas(strokeImage, &mergedImage);
            // showCurrentState();
            cout << "===================================================="
                 << endl;
        }
        return mergedImage;
    }

    /**
     * @brief it shows the current state, for the debug purpose
     * 
     * @return BRUSH_STATE
     **/
    BRUSH_STATE ReturnCurrentState() { return currentState_; }
    /**
     * @brief it shows the current state, for the debug purpose
     * 
     * @return BRUSH_STATE
     **/
    BRUSH_STATE ReturnLastState() { return lastState_; }
};