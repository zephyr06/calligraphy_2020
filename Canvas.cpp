/**
 * @file Canvas.cpp
 * @brief Definition and initialization of Canvas class
 * @author Frank Dellaert
 * @date September 2019
 */

#include "Canvas.h"

Canvas InitCanvas()
{
    Canvas temp(kImageHeight, kImageWidth);
    for (int i = 0; i < kImageHeight; i++)
        for (int j = 0; j < kImageWidth; j++)
            temp(i, j) = 255;
    return temp;
}

void MergeCanvas(
    const Canvas &imageToAdd,
    Canvas *inOutImage)
{ // merge the second into the first
    for (int i = 0; i < kImageHeight; i++)
        for (int j = 0; j < kImageWidth; j++)
            (*inOutImage)(i, j) = std::min((*inOutImage)(i, j), imageToAdd(i, j));
}
