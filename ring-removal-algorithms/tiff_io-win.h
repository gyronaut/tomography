/* Class to read and write tiff files */

#ifndef TIFF_IO_WIN_H
#define TIFF_IO_WIN_H
 
#pragma once

#include <stdlib.h>
#include <cstring>
#include <string>
#include "tiff.h"
#include "tiffio.h"

 class TiffIO
 {
 private:
	float min_val;
	float max_val;
 public:
	TiffIO();
	~TiffIO();
	float** readFloatImage(std::string input_name, int* w_ptr, int* h_ptr);
	float** read16bitImage(std::string input_name, int* w_ptr, int* h_ptr);
	void writeFloatImage(float** image, std::string output_name, int width, int height);
	void write16bitImage(float** image, std::string input_name, int width, int height); 
 };
 
 #endif
