/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * Image-based GTSAM optimization for Calligraphy writing trajectory
 */

#include "../CalligraphyPseudoSpectral.h"

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
        // 77F3-SHI; 4E00-YI; 4E7F-ZHI; 9E1F-NIAO; 4E8C-ER; 7a7a-KONG; 524d-QIAN
        unicode = Uunicode;
    else if (argc == 2)
        unicode = argv[1];
    else
    {
        cout << "Error: Too many input paramters!";
        return 0;
    }
    // ReadBrushSettings();

    // Initialize Chebyshev coefficients.
    string folder = "../../Generated_data_file/";
    string path_coeff = folder + "node_cheby_sample.txt";

    const auto initial_coeff = ReadCoefficients(path_coeff);
    vector<Canvas> single_strokes = file_io::ReadStrokePNGs(folder, unicode);

    const string path = "../../Generated_data_file/res_picture.txt";
    file_io::WriteCanvas(Calligraphy::DrawCharacter(initial_coeff),
                         path);
    const string path2 = "../../Generated_data_file/init_picture.txt";
    file_io::WriteCanvas(
        Calligraphy::DrawTrajectory(SampleCharacterPoints(initial_coeff)),
        path2);
    SaveSamplePoints(initial_coeff, unicode);
}
