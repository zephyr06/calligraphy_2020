#include "file_io.h"

#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>

using std::string;

namespace file_io
{

Canvas readPNG(string path)
{
    cv::Mat image;
    Canvas pict;
    image = cv::imread(path, 0);
    if (!image.data)
    {
        std::cout << path << std::endl;
        printf("No image data \n");
        return pict;
    }
    cv::cv2eigen(image, pict);
    std::cout << "Finish reading of one picture: " << path << std::endl;
    return pict;
}

Canvas ReadCharacterPNG(string path, string unicode)
{
    cv::Mat image;
    Canvas pict;
    string path_final = path + "ori.png";
    image = cv::imread(path_final, 0);
    if (!image.data)
    {
        printf("No image data \n");
        std::cout << path << std::endl;
        return pict;
    }
    cv::cv2eigen(image, pict);
    std::cout << "The character image: " << unicode << " is read" << std::endl;
    return pict;
}

void read_directory(const std::string &name, Filenames &v)
{
    DIR *dirp = opendir(name.c_str());
    struct dirent *dp;
    while ((dp = readdir(dirp)) != NULL)
    {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

std::vector<Canvas> ReadStrokePNGs(string path, string unicode)
{
    std::vector<Canvas> pngs;
    Filenames files;
    read_directory(path, files);
    sort(files.begin(), files.end());
    for (int i = 0; i < static_cast<int>(files.size()); i++)
    {
        int filename_len = files[i].length();
        // If the filename is too short, continue
        if (filename_len <= 4)
            continue;
        if (files[i].substr(filename_len - 4, filename_len - 1) == ".png" &&
            files[i].substr(0, 4) == unicode)
        {
            pngs.push_back(readPNG(path + files[i]));
        }
    }
    return pngs;
}

void WriteCanvas(const Canvas &picture, const string &path)
{
    std::ofstream file(path);
    if (file.is_open())
        file << picture;
    else
        std::cout << "Error opening the file, could not write to canvas."
                  << std::endl;
}
} // namespace file_io