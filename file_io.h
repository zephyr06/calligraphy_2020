/*
 * Copyright (C) 2019 The Borg Lab - All Rights Reserved
 *
 * Header for file Input/Output manipulation
 */

#include <string>
#include "Canvas.h"

namespace file_io
{

/**
 * @brief Used by ReadStrokePNGs() to read a single png file
 *
 * @param path path of PNG to read
 * @return Canvas image read
 */
Canvas readPNG(std::string path);

/**
 * Takes the path and unicode of a character image and
 * reads it as a Canvas, which is an image matrix.
 *
 * @param path path of character image
 * @param unicode unicode of character
 * @return Canvas image matrix of character image
 */
Canvas ReadCharacterPNG(std::string path, std::string unicode);

typedef std::vector<std::string> Filenames;

/**
 * @brief Read from necessary files in directory for GTSAM optimization
 *
 * @param name directory
 * @param v all file names
 */
void read_directory(const std::string &name, Filenames &v);

/**
 * @brief Used by testAnalyzeCharacter.cpp to get all images of individual strokes
 * from a directory
 *
 * @param path path at which to read PNG files
 * @param unicode looking for PNG related to this unicode
 * @return vector<Canvas> vector of images
 */
std::vector<Canvas> ReadStrokePNGs(std::string path, std::string unicode);

/**
 * @brief Write a picture to Generated_data_file folder as a text file.
 *
 * @param picture picture to write
 * @param path textfile name
 */
void WriteCanvas(const Canvas &picture, const std::string &path);
} // namespace file_io
