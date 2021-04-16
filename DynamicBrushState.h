
/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * The BrushState class for the dynamic virtual brush
 */

#include "CalligraphyBasis.h"
#include "VirtualBrush.h"

/**
 * @brief VBCrossSection represents the half ellipsoid shape that is drawn by the
 * Virtual Brush. This is also the intersection the VB makes with the paper;
 * parametric representation of the mark that is left on paper. Coefficients a,
 * b, c will help fit the shape of our contact surface.
 */
class VBCrossSection
{
private:
    int type; // type == 0: quadratic curve, y=a*x^2
    double a;
    double b;
    double c;

public:
    VBCrossSection(int type, double ai = 0, double bi = 0, double ci = 0)
        : type(type)
    {
        a = ai;
        b = bi;
        c = ci;
    }

    /**
     * @brief Calculating coefficient 'a' from y=a*x^2 given x (input1) and y(input2), if type=0
     * 
     * @param input1 the x in 'y=a*x^2'; one halft of the width
     * @param input2 ; drag
     * @param input3 Not used for now
     */
    void Fitting(double input1, double input2 = 0, double input3 = 0)
    {
        if (type == 0)
        {
            // accept input directly in the fitting frame, generally should be
            // (width/2, drag) if (input1 < 0.000001)

            //    cout << "Warning, poor fitting for VBCrossSection happened!"
            //         << endl;
            if (input1 < 0.01)
            {
                a = 0;
                c = 0;
                return;
            }
            a = input2 / input1 / input1;
            c = input2;
            // cout<<a<<", "<<c<<endl;
            return;
        }
        else
            cout << "Error type, please check available types!" << endl;
    }

    /**
     * @brief Calculate the expected output from the current curve expression
     * 
     * @param input1, the 'x'
     * @param input2 Not used for now
     * @param input3 Not used for now
     **/
    double Calculate(double input1, double input2 = 0, double input3 = 0)
    {
        if (type == 0)
        {
            return -a * input1 * input1 + c;
        }
        else
        {
            cout << "Error type, please check available types!" << endl;
            return 0;
        }
    }

    /**
     * @brief Find the drag based on the current curve, the assumed relationship is 
     * drag = a * (width/2)^2
     * 
     * @param width the width of the current written mark
     * @param input2 Not used for now
     * @param input3 Not used for now
     **/
    double Drag(
        double width, double input2 = 0,
        double input3 = 0)
    { // Accept width, and get corresponding drag
        if (type == 0)
        {
            return a * width * width / 4.0;
        }
        else
        {
            cout << "Error type, please check available types!" << endl;
            return 0;
        }
    }
    double Width(
        double drag, double input2 = 0,
        double input3 = 0)
    { // Accept width, and get corresponding drag
        if (type == 0)
        {
            return 2.0 * sqrt(drag / double(a));
        }
        else
        {
            cout << "Error type, please check available types!" << endl;
            return 0;
        }
    }
};

/**
 * @brief Find initial direction from stroke trajectory
 *
 * @param points sampled trajectory for stroke
 * @param threshold it evaluates if the current points 
 * formulate a set of valid points to judge the initial direction
 * @return the orientation angle in degrees
 */
double JudgeInitialDirection(const vector<Point3> &points,
                             double threshold = initialDirectionThreshold)
{
    if (points.size() < 5)
    {
        cout << "Error: too few points in JudgeInitialDirection" << endl;
        return 0;
    }
    vector<double> seq;
    for (int i = 0; i < static_cast<int>(points.size()) - 1; i++)
    {
        seq.push_back(atan2(-(points[i + 1][1] - points[0][1]),
                            points[i + 1][0] - points[0][0]) *
                      180 / M_PI);
        // cout << seq[i] << endl;
    }
    double judge = 0;
    for (int i = 0; i < static_cast<int>(seq.size()) - 4; i++)
    {
        judge =
            (abs(seq[i + 1] - seq[i]) +
             abs(seq[i + 2] - seq[i + 1]) + abs(seq[i + 3] - seq[i + 2]) +
             abs(seq[i + 4] - seq[i + 3])) /
            4.0;
        if (judge < threshold)
            return seq[i] + 180;
        // cout << i << ", " << seq[i] << endl;
    }
    cout << "Error: Didn't find initial direction, please increase the "
            "threshold!"
         << endl;
    return 0;
}

/**
 * @brief Find initial direction from stroke coefficients, a overloadded function
 *
 * @param strokeCoefficients the Coefficients way of representation of a stroke
 * @param threshold it evaluates if the current points 
 * formulate a set of valid points to judge the initial direction
 * @return the orientation angle in degrees
 */
double JudgeInitialDirection(const Coefficients &strokeCoefficients,
                             double threshold = initialDirectionThreshold)
{
    vector<Point3> points = SamplePoints(strokeCoefficients);
    return JudgeInitialDirection(points, threshold);
}

/**
 * @brief Virtual brush state used for the dynamic virtual brush
 */
class DynamicBrushState
{
public:
    double width;
    double drag;
    double offset;
    double currentDirection;
    double changeInDirection;
    double zLast;
    double x;
    double y;
    double z;
    vector<Point3> trajectoryHistory_; // record previous control command
    vector<double> offset_history_;

    int transitionPoint_;

    /**        
    scale controls the output range of z; higher scale results into
    smaller z; When scale=1, it corresponds to the actual-world
    simulation, though may not generate the best results since the
    accuracy of execution becomes more important
    **/
    // double kScale; // kScale=kImageHeight/character_size
    double inertiaInStroke;
    double inertiaDirection;

public:
    struct DataMember
    {
        double width;
        double drag;
        double offset;
        double currentDirection;
        double changeInDirection;
        double x;
        double y;
        double z;
        double zLast;
    };

    /** 
     * @brief it is used to copy all the state variables
     * 
     * @param Null
     * @return all of the data member of the class cloned
     **/
    DataMember clone()
    {
        DataMember data;
        data.width = width;
        data.drag = drag;
        data.offset = offset;
        data.currentDirection = currentDirection;
        data.changeInDirection = changeInDirection;
        data.x = x;
        data.y = y;
        data.z = z;
        data.zLast = z;
        return data;
    }

    // /**
    //  * @brief Construct a new State from a stroke direction
    //  *
    //  * @param inputState the Virtual Brush State object to copy
    //  */
    // DynamicBrushState &operator=(const DynamicBrushState &inputState)
    // {
    //     return *this;
    // }

    /**
     * @brief The constructor to initialize from another BrushState object
     * 
     * @param &t state with which to initialize
     **/
    DynamicBrushState(const DynamicBrushState &t) : width(t.width), drag(t.drag), offset(t.offset), changeInDirection(t.changeInDirection)
    {

        inertiaInStroke = t.inertiaInStroke;
        inertiaDirection = t.inertiaDirection;

        currentDirection = t.currentDirection;
        x = t.x;
        y = t.y;
        z = t.z;
        zLast = t.zLast;
        transitionPoint_ = t.transitionPoint_;
    }

    /**
     * @brief The default constructor
     **/
    DynamicBrushState() : width(0), drag(0), offset(0), changeInDirection(0)
    {

        inertiaInStroke = dInertiaInStroke;
        inertiaDirection = dInertiaDirection;

        currentDirection = 0;
        x = 0;
        y = 0;
        z = 0;
        zLast = 0;
        transitionPoint_ = 0;
        // cout << "Note! Somewhere called the default constructor!" << endl;
    }

    /**
     * @brief The constructor to initialize from a stroke
     * 
     * @param &firstStrokeCoefficients the Coefficients representation for a stroke
     **/
    DynamicBrushState(const Coefficients &firstStrokeCoefficients) : width(0), drag(0), offset(0), changeInDirection(0)
    {

        inertiaInStroke = dInertiaInStroke;
        inertiaDirection = dInertiaDirection;

        currentDirection = JudgeInitialDirection(firstStrokeCoefficients);
        x = 0;
        y = 0;
        z = 0;
        zLast = 0;
        transitionPoint_ = 0;
    }

    /**
     * @brief The constructor to initialize from a stroke
     * 
     * @param firstStrokePoints the vector representation for a stroke
     **/
    DynamicBrushState(const vector<Point3> &firstStrokePoints) : width(0), drag(0), offset(0), changeInDirection(0)
    {

        inertiaInStroke = dInertiaInStroke;
        inertiaDirection = dInertiaDirection;

        currentDirection = JudgeInitialDirection(firstStrokePoints);
        x = 0;
        y = 0;
        z = 0;
        zLast = 0;
        transitionPoint_ = 0;
    }

    //The following functions are about simple updating methods for the parameters of dynamic virtual brush

    /**
     * @brief it estimates the drag parameter given the 'z' input in the ideal situation
     * 
     * @param z the current z height, z=0 means the tip of the brush just contacts the paper with not deformation, 
     * z>0 means it deforms after contacting the paper
     * @return the length of drag
     **/
    double Drag(
        double z) // The length of the tip that have contact with the paper
    {

        double res = 0;
        vector<double> coeff{0.66453417, 0.00103169};
        if (z <= 0)
            return 0;
        else
            res = kScaleReal2Image * (coeff[0] * z + coeff[1]);
        if (res < 0.01)
        {
            return 0;
        }
        return res;
    }

    /**
     * @brief it estimates the width parameter given the 'z' input in the ideal situation
     * 
     * @param z the current z height, z=0 means the tip of the brush just contacts the paper with not deformation, 
     * z>0 means it deforms after contacting the paper
     * @return the length of width
     **/
    double Width(
        double z) // The length of the tip that have contact with the paper
    {

        double res = 0;
        vector<double> coeff{1.19034619e+00, -3.20286161e-04};
        if (z <= 0)
            return 0;
        else
            res = kScaleReal2Image * (coeff[0] * z + coeff[1]);
        if (res < 0.01)
        {
            return 0;
        }
        return res;
    }

    /**
     * @brief it estimates the offset parameter given the 'z' input in the ideal situation
     * 
     * @param z the current z height, z=0 means the tip of the brush just contacts the paper with not deformation, 
     * z>0 means it deforms after contacting the paper
     * @return the length of offset
     **/
    static double Offset(
        double z)
    {

        double res = 0;
        // Remember to modify the inverse_offset() function, if you want to change there!!!!!!
        vector<double> coeff{-3.35623343e+01, 7.30804484e-01, -2.11326069e-04};
        if (z <= 0)
            return 0;
        else if (z < 0.030)
            res = kScaleReal2Image * (coeff[0] * z * z + coeff[1] * z + coeff[2]);
        if (res < 0.01)
        {
            return 0;
        }
        return res;
    }

    /**
     * @brief Old method about updating the orientation of the brush mark.
     * All calculations are done in the image coordinate frame where x-axis is
     * horizontal and y-axis is vertical.
     *
     * @param nextPoint the control point that updates the orientation parameter
     * @return the updated orientation of the brush, in degrees
     */
    double Orientation(Point3 nextPoint)
    {
        // Angles used in the degree form in this function, sin/cos/tan function
        // requires radius form

        int len = trajectoryHistory_.size();

        // No trajectoryHistory_, beginning of stroke has same orientation as
        // ending of the previous stroke
        if (len == 0)
        {
            return currentDirection;
        }
        else
        {
            Point2 brushMovement;
            // Not sure if nextPoint.x() corresponds to horizontal or vertical
            // change, same with nextPoint.y()
            brushMovement[0] = nextPoint.x() - trajectoryHistory_[len - 1].x();
            brushMovement[1] = nextPoint.y() - trajectoryHistory_[len - 1].y();

            double brushMovementAngle =
                (atan2(-brushMovement[1], brushMovement[0]) * 180 / CV_PI +
                 180);
            if (brushMovementAngle > 180)
                brushMovementAngle -= 360;
            // cout << "Moving angle: " << brushMovementAngle << endl;
            return brushMovementAngle;
        }
    }

    /**
     * @brief Updating the orientation of the brush mark.
     * All calculations are done in the image coordinate frame where x-axis is
     * horizontal and y-axis is vertical.
     *
     * @param nextPoint the control point that updates the orientation parameter
     * @param rootPoint the root point at the last time instance
     */
    void Orientation(Point3 nextPoint, Point3 rootPoint)
    {
        // The brush does not touch the paper, nothing will be updated
        if (nextPoint.z() <= 0)
            return;

        Point2 moveDisplacement(nextPoint.x() - rootPoint.x(),
                                nextPoint.y() - rootPoint.y());
        double moveDistMax = moveDisplacement.norm();
        if (moveDistMax > offset + drag * centerRootRatio)
        {
            currentDirection = atan2(moveDisplacement.y() * -1, moveDisplacement.x()) *
                                   180 / CV_PI +
                               180;
        }
        else
        {
            offset = moveDistMax;
            currentDirection = atan2(moveDisplacement.y() * -1, moveDisplacement.x()) *
                                   180 / CV_PI +
                               180;
        }

        if (currentDirection > 180)
            currentDirection -= 360;
        else if (currentDirection < -180)
            currentDirection += 360;
    }

    /**
     * @brief Calculate the position of root point at last time instance; 
     *  the definition of root point can be found in the paper; 
     * This function is always called before anything is updated
     *
     * @return the root point at the last time instance
     */
    Point3 GetRootPoint()
    {
        return Point3(x + (offset +  drag* centerRootRatio) * cos(currentDirection / 180.0 * CV_PI),
                      y - (offset +  drag* centerRootRatio) * sin(currentDirection / 180.0 * CV_PI),
                      z);
    }

    /**
     * @brief update the dynamic parameters after each control nextPoint;
     * This version of Update assumes that the brush dips ink and restores itself everytime before it writes a stroke
     * 
     * @param nextPoint the control point that updates the parameters of the virtual brush
     * @param type: could only be "between" or "in" or "initial" or "last"
     */
    void Update(const Point3 &nextPoint, string type)
    {

        // Must run this at first
        Point3 rootPoint = GetRootPoint();
        // cout << "The root point: " << rootPoint << endl;
        x = nextPoint[0];
        y = nextPoint[1];
        z = nextPoint[2];

        // DataMember lastState_ = clone();

        int len = trajectoryHistory_.size();
        if (z <= 0 || len == 0)
        {
            trajectoryHistory_.push_back(nextPoint);
            return;
        }

        //Note: the orientation parameter must be updated at last
        if (type == "in")
        {
            width = width * inertiaInStroke +
                    Width(z) * (1 - inertiaInStroke);
            drag = drag * inertiaInStroke +
                   Drag(z) * (1 - inertiaInStroke);
            offset = offset * inertiaInStroke +
                     Offset(z) * (1 - inertiaInStroke) * kScaleOffset;

            Orientation(nextPoint, rootPoint);
        }
        else if (type == "last")
        {
            // lastState_.currentDirection = currentDirection;
            offset = 0;
            // double z_last = inverse_offset(lastState_.offset);
            zLast = 0;
            width = 0;
            drag = 0;
            changeInDirection = 0;

            trajectoryHistory_.clear();
            offset_history_.clear();

            return;
        }
        else
        {
            cout << "Type Input Error, only 'between', 'in' and 'initial' are "
                    "recognizable!"
                 << endl;
            return;
        }

        trajectoryHistory_.push_back(nextPoint);
        offset_history_.push_back(offset);
        // cout << "TEST-current=================" << offset <<
        // endl;

        // showCurrentState();
        // cout << "Current Point: " << nextPoint << endl;
        return;
    }

    //The following two functions are used to find the last state of virtual brush
    /**
     * @brief Estimate the offset of the brush when it finishes one stroke and get lifted upon the paper
     * 
     * @param execlude_rate The points at the beginning/ending part of the stroke are excluded for 
     * consideration since they generally suffer from the errors of high-order polynomials approximation
     * @return The estimated offset
     **/
    double find_offset_last(double execlude_rate = 0.05)
    {
        // Calculate the average offset of the stroke, and use it for the
        // lastState_
        double sum_offset = 0;
        int index_begin = round(offset_history_.size() * execlude_rate);
        int index_end = static_cast<int>(offset_history_.size()) - index_begin;

        for (int i = index_begin; i < index_end; i++)
            sum_offset += offset_history_[i];

        double res = sum_offset / (index_end - index_begin);
        return res;
    }
    /**
     * @brief Find the estimated 'z' from the estimated offset, 
     * the 'z' can be used to calculate the corresponding width and drag parameters
     * 
     * @param offset the estimated offset from find_offset_last
     * @return the estimated 'z'
     **/
    static double inverse_offset(double offset)
    {
        vector<double> coeff{-3.35623343e+01, 7.30804484e-01, -2.11326069e-04};
        if (offset < coeff[2])
        {
            cout << "No solution found in inverse_offset, please check it out!" << endl;
            return 0;
        }
        double a = coeff[0];
        double b = coeff[1];
        double c = coeff[2] - offset / kScaleReal2Image;
        double z1 = (-b + sqrt(b * b - 4 * a * c)) / 2.0 / a;
        // double z2 = (-b - sqrt(b * b - 4 * a * c)) / 2.0 / a;
        // cout << z1 << endl;
        // cout << z2 << endl;
        // cout << Offset(z1) << endl;
        // cout << Offset(z2) << endl;
        if (z1 > 0)
            return z1;
        else
        {
            cout << "No valid solution found in inverse_offset, please check it out!" << endl;
            return 0;
        }
    }

    /**
     * @brief Show the current states of the dynamic virtual brush
     **/
    void showCurrentState() const
    {
        cout << "********************************************" << endl;
        // cout << "Current Virtual Brush state variables: " << endl;
        cout << "width:                   " << width << endl;
        cout << "drag:                    " << drag << endl;
        cout << "offset:                  " << offset << endl;
        cout << "currentDirection:        " << currentDirection
             << endl;
        cout << "changeInDirection:       " << changeInDirection
             << endl;
    }

    /**
     * @brief Edit the current states of the dynamic virtual brush
     * 
     * @param a_width new width
     * @param b_drag new drag
     * @param c_offset new offset
     * @param d_currentDirection new Direction
     * @param e_changeInDirection new changeInDirection
     **/
    void Edit_current_state(double a_width, double b_drag, double c_offset, double d_currentDirection, double e_changeInDirection)
    {
        width = a_width;
        drag = b_drag;
        offset = c_offset;
        currentDirection = d_currentDirection;
        changeInDirection = e_changeInDirection;
    }
};
