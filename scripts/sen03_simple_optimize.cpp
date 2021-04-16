/*
 * @file sen03_optimize.cpp
 * @brief Image-based GTSAM optimization for Calligraphy writing trajectory
 * 
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 */

#include "../CalligraphyPseudoSpectral.h"

// Toggle between using SimpleVirtualBrush or EnhancedDynamicVirtualBrush

#include "../SimpleVirtualBrush.h"
using Calligraphy = CalligraphyPseudoSpectral<SimpleVirtualBrush,
                                              SimpleBrushState>;

/**
 * Create and use factor graph for Chinese Calligraphy stroke optimization.
 */
int main(int argc, char *argv[])
{
    // Get correct unicode for character to optimize.
    string unicode;
    if (argc == 1)
        // 77F3-SHI; 4E00-YI; 4E7F-ZHI; 9E1F-NIAO; 4E8C-ER; 7A7A-KONG; 524D-QIAN
        unicode = Uunicode;
    else if (argc == 2)
        unicode = argv[1];
    else
    {
        cout << "Error: Too many input paramters!";
        return 0;
    }

    // ReadBrushSettings();
    // Initialize Chebyshev coefficients by reading from file
    // "node_cheby_sample.txt"
    string folder = "../../Generated_data_file/";
    string coefficientsPath = folder + "node_cheby_sample.txt";

    // Call to coefficients reader from "CalligraphyBasis.h"
    const auto initialCoefficients = ReadCoefficients(coefficientsPath);

    // Reading the original image as well as individual stroke images
    Canvas originalCharacterImage = file_io::ReadCharacterPNG(folder, unicode);
    vector<Canvas> originalStrokeImages =
        file_io::ReadStrokePNGs(folder, unicode);

    // Optimize******************************************
    auto optimizeCoefficients = Calligraphy::OptimizeCharacter(
        initialCoefficients, originalStrokeImages);

    Calligraphy::CharacterErrorAnalysis error_analysis(
        originalCharacterImage, initialCoefficients, optimizeCoefficients);
    error_analysis.run();

    // Saving sampled points based on optimized coefficients
    SaveSamplePoints(optimizeCoefficients, unicode);

    // Save the image of the optimized skeleton
    const string optimizedSkeletonPicturePath =
        "../../Generated_data_file/skeleton_picture.txt";
    file_io::WriteCanvas(Calligraphy::DrawTrajectory(
                             SampleCharacterPoints(optimizeCoefficients)),
                         optimizedSkeletonPicturePath);

    // Save sample points from initial trajectory
    // Calligraphy::SaveSamplePoints(initialCoefficients, unicode);
}
