/* Class to read and write tiff files */

#ifndef TIFF_IO_H
#define TIFF_IO__H
 
#pragma once

#include <string>
#include "tiff.h"
#include "tiffio.h"

 class TiffIO
 {
 private:
	
 public:
	TiffIO();
	~TiffIO();
	float** readFloatImage(string input_name, int* w_ptr, int* h_ptr);
	void writeFloatImage(float** image, string output_name, int width, int height);
 
 };
 
 #endif