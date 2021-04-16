/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * Header for image-based GTSAM optimization for Calligraphy writing trajectory
 */

#pragma once

#include "CalligraphyBasis.h"

/**
 * Calligraphy writing based on pseudo-spectral basis
 *
 * The template argument BRUSH should be class that has
 *   A type BRUSH::State that formalizes the continuous brush state, if any
 *   Constructor BRUSH that takes a BRUSH::State
 *
 */
template <class BRUSH, class BRUSH_STATE>
class CalligraphyPseudoSpectral
{
public:
    /**
     * @brief It draws a stroke given one stroke represented as Coefficients
     * 
     * @param coeff the Coefficients of a stroke
     * @param lastState_ the BrushState to initialize the virtual brush
     * @return The Canvas of the stroke picture
     **/
    static Canvas DrawStroke(Coefficients coeff,
                             const BRUSH_STATE &lastState_)
    {
        // Problem: the new virtual brush is totally initialized every time the
        // function is called
        vector<Point3> points = SamplePoints(coeff);
        // We will search for the direction of the first point later
        BRUSH vb(lastState_);
        Canvas picture = vb.DrawStroke(points);
        return picture;
    }

    /**
     * @brief Given the set of stroke trajectories, draw the character
     *
     * @param characterCoefficients vector of Coefficients representing
     * @param showState If show the state of the virtual brush after it draws each stroke or not
     * Chebyshev polynomials of x, y, z coordinates
     * @return Canvas: character drawn on a Canvas
     */
    static Canvas DrawCharacter(
        const vector<vector<Point3>> &characterTrajectories,
        bool showState = false)
    {
        Canvas picture;

        const BRUSH_STATE lastState_(characterTrajectories[0]);
        // const BRUSH_STATE lastState_;
        BRUSH vb(lastState_);
        if (showState)
            picture = vb.DrawCharacterWithState(characterTrajectories);
        else
            picture = vb.DrawCharacter(characterTrajectories);

        return picture;
    }

    /**
     * @brief Given the set of Chebyshev coefficients, one for each stroke, draw the
     * character
     *
     * @param characterCoefficients: vector of Coefficients representing
     * Chebyshev polynomials of x, y, z coordinates
     * @param showState If show the state of the virtual brush after it draws each stroke or not
     * Chebyshev polynomials of x, y, z coordinates
     * @return Canvas: character drawn on a Canvas
     */
    static Canvas DrawCharacter(
        const vector<Coefficients> &characterCoefficients, int localNUM_SAMPLES_Draw = NUM_SAMPLES_Draw,
        bool showState = false)
    {
        const auto characterTrajectories =
            SampleCharacterPoints(characterCoefficients, localNUM_SAMPLES_Draw);
        return DrawCharacter(characterTrajectories, showState);
    }

    /**
     * @brief Draw the trajectory on the canvas
     * @param characterTrajectories The character trajectories, composed of the points of each stroke
     * @return the Canvas of the character
     **/
    static Canvas DrawTrajectory(
        const vector<vector<Point3>> &characterTrajectories)
    {
        Canvas trajectory_map = InitCanvas();
        for (const auto &strokeTrajectory : characterTrajectories)
        {
            for (const Point3 &point : strokeTrajectory)
            {
                int i = round(point.y());
                int j = round(point.x());
                if (i >= 0 && i < kImageHeight && j >= 0 && j < kImageWidth)
                {
                    trajectory_map(i, j) = 0;
                }
            }
        }
        return trajectory_map;
    }

    /**
     * Factor for Chebyshev coefficients
     */
    class ChebyshevFactor : public NoiseModelFactor1<Coefficients>
    {
    private:
        Canvas targetPicture_;
        BRUSH_STATE lastState_;

    public:
        /**
         * @brief The constructor for the factor
         * @param key the key for a factor
         * @param targetPicture The reference character image
         * @param lastState The BrushState to initialize the virtual brush used in thefactor
         * @param model The model for GTSAM
         **/
        ChebyshevFactor(Key key, const Canvas targetPicture,
                        const BRUSH_STATE lastState, SharedNoiseModel model)
            : NoiseModelFactor1<Coefficients>(model, key),
              targetPicture_(targetPicture),
              lastState_(lastState) {}

        /**
         * @brief Evaluates the error between our current image drawing and desired image
         **/
        Vector evaluateError(const Coefficients &point,
                             boost::optional<Matrix &> H = boost::none) const override
        {
            PictError v;
            // Optimizing per stroke
            Canvas currentPicture = DrawStroke(point, lastState_);
            v = PictureDifference(currentPicture, targetPicture_);
            vector<Point3> sPoints = SamplePoints(point);
            uint num = sPoints.size();
            for (uint i = 0; i < 50; i++)
                v(127 - i, 127) = sPoints[num - i - 1].z() * 10000000 * scaleForRegularizationZ * (50 - i) / 50.0;

            if (H)
            {
                boost::function<Matrix(const Coefficients &)> f =
                    [this](const Coefficients &input) {
                        Canvas temp_character = DrawStroke(input, lastState_);
                        PictError v = PictureDifference(temp_character,
                                                        targetPicture_);
                        vector<Point3> sPoints = SamplePoints(input);
                        uint num = sPoints.size();
                        for (uint i = 0; i < 50; i++)
                            v(127 - i, 127) = sPoints[num - i - 1].z() * 10000000 * scaleForRegularizationZ * (50 - i) / 50.0;
                        return v;
                    };

                *H = numericalDerivative11(f, point,
                                           deltaOptimizer); // 0.5 before works better
            }
            return v;
        }
    };

    /**
     * @brief It finds the last state of the virtual brush after it writes a stroke
     * 
     * @param coeff The Coeffitients for a stroke
     * @param lastState_ the Virtual Brush State to initialize a virtual brush
     * @return the BrushState after the virtual brush writes a stroke
     **/
    static BRUSH_STATE FindLastState(Coefficients coeff,
                                     BRUSH_STATE lastState_)
    {
        BRUSH vb(lastState_);
        vb.DrawStroke(coeff);
        return vb.ReturnLastState();
    }

    /**
     * @brief Run GTSAM optimizer on a specific stroke
     *
     * @param initialCoefficients Initial coefficients for a specific stroke
     * @param lastState_ The last state of the brush from the last previous stroke
     * @param strokeImage Ideal stroke image that we want to draw
     * @return tuple<Coefficients, BRUSH_STATE>
     */
    static tuple<Coefficients, BRUSH_STATE> OptimizeStroke(
        Coefficients initialCoefficients, BRUSH_STATE lastState_,
        Canvas strokeImage)
    {
        Symbol key('a', 0);
        Values initialEstimate;
        auto model =
            noiseModel::Isotropic::Sigma(kImageWidth * kImageHeight, deltaOptimizer);
        cout << "Creating the factor graph......" << endl;

        NonlinearFactorGraph graph;
        graph.emplace_shared<ChebyshevFactor>(key, strokeImage, lastState_,
                                              model);
        initialEstimate.insert(key, initialCoefficients);

        cout << "Begin the optimization......" << endl;
        LevenbergMarquardtParams params;
        params.setlambdaInitial(initialLambda);
        params.setVerbosityLM("SUMMARY");
        // params.setLinearSolverType("SEQUENTIAL_QR");
        params.setLinearSolverType(optimizerType);
        // params.setLinearSolverType("MULTIFRONTAL_CHOLESKY");
        // params.setLinearSolverType("SEQUENTIAL_CHOLESKY");
        LevenbergMarquardtOptimizer optimizer(graph, initialEstimate, params);

        Values result = optimizer.optimize();
        Coefficients resultCoefficients = result.at<Coefficients>(key);
        BRUSH_STATE lastStateAfterOptimize =
            FindLastState(resultCoefficients, lastState_);

        // cout << resultCoefficients << endl;
        return std::make_tuple(resultCoefficients, lastStateAfterOptimize);
    }

    /**
     * @brief Optimizes for strokes in a character
     *
     * @param initialCoefficients the initial estimation of the Coefficients of the strokes in a character
     * @param originalStrokeImages the reference picture set of all the strokes of a character
     * @return vector<Coefficients>
     */
    static vector<Coefficients> OptimizeCharacter(
        vector<Coefficients> initialCoefficients,
        vector<Canvas> originalStrokeImages)
    {
        if (static_cast<int>(initialCoefficients.size()) < 1)
        {
            cout << "Input error, the character must have at least one stroke!";
            return initialCoefficients;
        }
        vector<Coefficients> resultCoefficients;
        Coefficients tempCoefficients;
        // bool dipInk=true;

        // Provide an initialization for the first stroke, this method should
        // be different for different type of virtual brushes
        BRUSH_STATE lastState_(initialCoefficients[0]);

        // Optimize for all the strokes
        for (int i = 0; i < static_cast<int>(initialCoefficients.size()); i++)
        {
            if(dipInk)
            {
                BRUSH_STATE beginState_(initialCoefficients[i]);
                tie(tempCoefficients, lastState_) = OptimizeStroke(
                    initialCoefficients[i], beginState_, originalStrokeImages[i]);
                cout << " Finish optimization of one stroke in dipping ink situation......" << endl;                
            }
            else
            {
                tie(tempCoefficients, lastState_) = OptimizeStroke(
                    initialCoefficients[i], lastState_, originalStrokeImages[i]);
                cout << " Finish optimization of one stroke not in dipping ink situation......" << endl;
            }

            resultCoefficients.push_back(tempCoefficients);
        }
        cout << "******************* Character Optimization finished!" << endl;
        return resultCoefficients;
    }

    /**
     * @brief Calculate the error between the original and optimized image.
     * 
     */
    class CharacterErrorAnalysis
    {
    private:
        Canvas characterTargetImage_;
        vector<Coefficients> initialCoefficients_;
        vector<Coefficients> resultCoefficients_;

    public:
        /**
         * @brief The constructor
         * @param characterTargetImage The reference character image
         * @param initialCoefficients The initial estimation of the Coefficients in the character
         * @param resultCoefficients The Coefficients gotten after optimization
         **/
        CharacterErrorAnalysis(Canvas characterTargetImage,
                               vector<Coefficients> initialCoefficients,
                               vector<Coefficients> resultCoefficients)
            : characterTargetImage_(characterTargetImage),
              initialCoefficients_(initialCoefficients),
              resultCoefficients_(resultCoefficients) {}

        /**
         * @brief Run all error evaluations
         **/
        void run() const
        {
            initialError();
            currentError();
            cout << "The character picture is saved now." << endl;
        }

        /**
         * @brief Calculates error between initial estimation and reference image
         **/
        void initialError() const
        {
            Canvas temp = DrawCharacter(initialCoefficients_, NUM_SAMPLES_Draw);
            cout << "Total initial error: "
                 << PictureDifference(temp, characterTargetImage_).norm()
                 << endl;
            file_io::WriteCanvas(temp,
                                 "../../Generated_data_file/init_picture.txt");
        }

        /**
         * @brief Calculates error between optimized results and reference image
         **/
        void currentError() const
        {
            Canvas temp = DrawCharacter(resultCoefficients_, NUM_SAMPLES_Draw);
            cout << "Saving resulting image into res_picture.txt" << endl;
            file_io::WriteCanvas(temp,
                                 "../../Generated_data_file/res_picture.txt");
            cout << "Total current error: "
                 << PictureDifference(temp, characterTargetImage_).norm()
                 << endl;
        }
    };
};
