#include <iostream>
#include <math.h>

#include "../CppUnitTestFramework.hpp"
#include "PixelDifferenceMethods.cpp"

#define GENERATE_UNIT_TEST_MAIN
using namespace std;

int main()
{

    // vector<std::string> listOfFiles = {"Kong-initial.png",
    //                                    "Kong-opt-dynamic.png"};

    vector<std::string> listOfFiles = {"Si-opt-dynamic1.png",
                                       "Si-opt-dynamic2.png",
                                       "Si-opt-dynamic3.png",
                                       "Si-opt-simpleBrush.png",
                                       "Si-initial.png"};

    // vector<std::string> listOfFiles = {"Wo-initial.png",
    //                                    "Wo-opt-dynamic.png",
    //                                    "Wo-opt-dynamic1.png",
    //                                    "Wo-opt-dynamic2.png"};

    // vector<std::string> listOfFiles = {"si_init_simpleBrush.png",
    //                                    "si_opt_simpleBrush.png",
    //                                    "simu_init_dynamicBrush.png",
    //                                    "simu_opt_dynamicBrush.png"};

    std::string originalImagePath = "../Data/svgtopng/kong_svg2png.png";

    int canvasDim = 640;

    for (string filename : listOfFiles)
    {
        string resultImagePath = "../Data/written_results_scanned_using_scanner/" + filename;
        cout << resultImagePath << endl;
        double pixelDifference = processImages(originalImagePath, resultImagePath, filename, canvasDim);

        cout << "The count of pixels different for file " << resultImagePath << " is: " << pixelDifference << endl;
        cout << "Pixel difference as a percentage of " << canvasDim << " x " << canvasDim << " canvas: " << pixelDifference / (canvasDim * canvasDim) * 100 << "%" << endl;

        cout << endl;
    }

    return 0;
}
