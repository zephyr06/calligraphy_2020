
#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <opencv2/core/core.hpp>

cv::FileStorage ConfigParameters("/home/zephyr/Programming/calligraphy_8_1/parameters.yaml", cv::FileStorage::READ);
// cv::FileStorage ConfigParameters("/home/lab/Documents/calligraphy_8_1/parameters.yaml", cv::FileStorage::READ);
using namespace std;
const int kNodeNum = 7;

const double kZ_Tolerance = (double)ConfigParameters["kZ_Tolerance"]; // The tolerance of effective z to generate friction with paper
// Note, if you want to try different number of nodes, you may need to change
// the initialization in the Node constructors
// static const std::vector<double> segmentSeq{7, 5, 5, 3, 2, 1, 0.5}; // It has KNodeNum elements,

// Note: the definition of baseLength is associated with Node. The position of the first node
// is given by the control command, Seven nodes only need 6 segment lengths,
// the same as the baseLength vector. The node is locked with the segment above it.

// static const vector<float> baseLength{7, 5, 3, 3, 2, 1};

vector<float> baseLength;
vector<float> nodeStrainWeights;
vector<float> nodeGeometricWeights;

// each node has an associated segment
const double muExternal = (double)ConfigParameters["muExternal"];

const double weightVerticalInnerFriction = (double)ConfigParameters["weightVerticalInnerFriction"];
const double weightStrain = (double)ConfigParameters["weightStrain"];
const double weightStrainXY = (double)ConfigParameters["weightStrainXY"]; // used for the strain energy in StrainEnergyFactor
const double weightGeometric = (double)ConfigParameters["weightGeometric"];
const double weightZ = (double)ConfigParameters["weightZ"]; //0.5 for exponential mapping
const double weightStrainHorizontal = (double)ConfigParameters["weightStrainHorizontal"];
const double weightTwist = (double)ConfigParameters["weightTwist"];

const double kMax_Gradient = (double)ConfigParameters["kMax_Gradient"];
const double initializationStep = (double)ConfigParameters["initializationStep"];
const double kRatioDrag = (double)ConfigParameters["kRatioDrag"];

const double mergeRatio = (double)ConfigParameters["mergeRatio"];
const double initialLambda = (double)ConfigParameters["initialLambda"];

const double kMinGradientTolerance = (double)ConfigParameters["kMinGradientTolerance"];
const double ZGeometricTol = (double)ConfigParameters["ZGeometricTol"];
const double weightFrictionStrainHorizontal = (double)ConfigParameters["weightFrictionStrainHorizontal"];

const double tipAngle = (double)ConfigParameters["tipAngle"];

const double kStrainLateralForOrientation = (double)ConfigParameters["kStrainLateralForOrientation"];
const double deltaOptimizer = (double)ConfigParameters["deltaOptimizer"];
const double alphaClamp = (double)ConfigParameters["alphaClamp"];
const double weightPorResistance = (double)ConfigParameters["weightPorResistance"];
const int strainOrder = (double)ConfigParameters["strainOrder"];

const double dInertiaInStroke = (double)ConfigParameters["dInertiaInStroke"];
const double dInertiaDirection = (double)ConfigParameters["dInertiaDirection"];
const string Uunicode = (string)ConfigParameters["Uunicode"];
const string optimizerType = (string)ConfigParameters["optimizerType"];

// const int ORDER = (int)ConfigParameters["ORDER"];
const int NUM_SAMPLES = (int)ConfigParameters["NUM_SAMPLES"];
const int NUM_SAMPLES_Draw = (int)ConfigParameters["NUM_SAMPLES_Draw"];
const double initialDirectionThreshold = (double)ConfigParameters["initialDirectionThreshold"];

const double scaleForRegularizationZ = (double)ConfigParameters["scaleForRegularizationZ"];
const double kScaleOffset = (double)ConfigParameters["kScaleOffset"];
const double centerRootRatio = (double)ConfigParameters["centerRootRatio"];
const int dipInk = (int)ConfigParameters["dipInk"];


// const double kCharacterSize = 0.256; // should be 0.256
const double kCharacterSize = 0.12; // should be 0.256
const double kScaleReal2Image = 128 / kCharacterSize;
const int ORDER = 6;
