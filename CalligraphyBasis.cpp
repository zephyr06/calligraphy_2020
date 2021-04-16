/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * Includes the basic functions that do not depend on any type of brush or state
 */

#include "CalligraphyBasis.h"

vector<double> GenerateChebyshevPoints(int numPoints)
{
    vector<double> result;
    result.reserve(numPoints);
    const double dt = PI / (numPoints - 1);
    for (int i = 0; i < numPoints - 1; i++)
    {
        result.push_back(-cos(i * dt));
    }
    result.push_back(1.0);
    return result;
}

vector<double> GenerateUniform(int numPoints)
{
    vector<double> result;
    result.reserve(numPoints);
    const double dx = 2.0 / (numPoints - 1);
    for (int i = 0; i < numPoints - 1; i++)
    {
        result.push_back(-1.0 + i * dx);
    }
    result.push_back(1.0);
    return result;
}

vector<Point3> SamplePoints(Coefficients coefficients, int sampleNum, bool uniform)
{
    vector<double> times = uniform ? GenerateUniform(sampleNum)
                                   : GenerateChebyshevPoints(sampleNum);

    auto X = coefficients.block(ORDER * 0, 0, ORDER, 1);
    auto Y = coefficients.block(ORDER * 1, 0, ORDER, 1);
    auto Z = coefficients.block(ORDER * 2, 0, ORDER, 1) / 10000;

    vector<Point3> points;
    points.reserve(sampleNum);
    for (auto &&t : times)
    {
        typename Chebyshev2<ORDER>::EvaluationFunctor fx(t);
        points.emplace_back(Point3(fx(X), fx(Y), fx(Z)));
    }
    return points;
}

PictError PictureDifference(Canvas generatedPicture, Canvas targetPicture)
{
    PictError e;
    for (int i = 0; i < kImageHeight; i++)
    {
        for (int j = 0; j < kImageWidth; j++)
        {
            e(i * kImageWidth + j) =
                abs(generatedPicture(i, j) - targetPicture(i, j));
        }
    }
    return e;
}

vector<Coefficients> ReadCoefficients(string path, double ratio)
{
    ifstream infile(path);
    vector<Coefficients> coeff_vec;

    string line;
    string token;
    string delimiter = " ";
    size_t pos = 0;

    int line_num = 0;
    int inner_index = 0;
    Coefficients coeff_temp;

    if (infile.is_open())
    {
        cout << "The coefficients file opened!" << endl;
        Coefficients coeff_temp;
        while (getline(infile, line))
        {
            while ((pos = line.find(delimiter)) != string::npos)
            {
                token = line.substr(0, pos);
                float temp =
                    atof(token.c_str()) / ratio; // Transform coefficients from
                                                 // 1024-scale to 128-scale
                coeff_temp(line_num % 3 * ORDER + inner_index, 0) = temp;
                line.erase(0, pos + delimiter.length());
                inner_index++;
            }
            coeff_temp(line_num % 3 * ORDER + inner_index, 0) =
                atof(line.c_str()) / ratio; // For the last element

            if (line_num % 3 == 2 && inner_index == ORDER - 1)
            {
                coeff_vec.push_back(coeff_temp);
                Coefficients coeff_temp;
            }
            line_num++;
            inner_index = 0;
        }
    }
    else
    {
        throw std::invalid_argument("Can't open the coefficients file " + path);
    }

    cout << "Finished reading coefficients file; read " << coeff_vec.size()
         << " strokes." << endl;
    return coeff_vec;
}

void SaveSamplePoints(vector<Coefficients> resultCoefficients,
                      string unicode)
{
    vector<Point3> result;
    result.reserve(resultCoefficients.size());

    for (int i = 0; i < static_cast<int>(resultCoefficients.size()); i++)
    {
        vector<Point3> sampledPoints = SamplePoints(resultCoefficients[i]);
        result.insert(result.end(), sampledPoints.begin(), sampledPoints.end());
    }

    // Save under unique filename "<unicode>_Sample_points.txt"
    string uniquePath =
        "../../Generated_data_file/" + unicode + "_Sample_points.txt";
    ofstream file(uniquePath);
    if (file.is_open())
    {
        for (int i = 0; i < static_cast<int>(result.size()); i++)
            file << result[i] << endl;
    }
    else
    {
        cout << "Error opening the Sample_points.txt file, nothing was written."
             << endl;
    }

    // Save under generic filename "Sample_points.txt"
    string genericPath = "../../Generated_data_file/Sample_points.txt";
    ofstream file2(genericPath);
    if (file2.is_open())
    {
        for (int i = 0; i < static_cast<int>(result.size()); i++)
            file2 << result[i] << endl;
    }
    else
    {
        cout << "Error opening the Sample_points.txt file, nothing was written."
             << endl;
    }
}

vector<vector<Point3>> SampleCharacterPoints(
    const vector<Coefficients> &characterCoefficients, int localNUM_SAMPLES_Draw)
{
    vector<vector<Point3>> characterTrajectories;
    characterTrajectories.reserve(characterCoefficients.size());
    for (const auto &strokeCoefficients : characterCoefficients)
    {
        characterTrajectories.push_back(SamplePoints(strokeCoefficients, localNUM_SAMPLES_Draw));
    }
    return characterTrajectories;
}
