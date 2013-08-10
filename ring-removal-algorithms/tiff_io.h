/* Class to read and write tiff files */

#ifndef TIFF_IO_H
#define TIFF_IO_H
 
#pragma once

#include <stdlib.h>
#include <cstring>
#include <string>
#include "tiff.h"
#include "tiffio.h"

 class TiffIO
 {
 private:
	static const int NUM_HISTOGRAM_BINS = 256;
	
	float min_val;
	float max_val;
	int width;
	int height;
	
	void createHistogram(float** image_rows);
 public:	
	int* histogram;
	
	TiffIO();
	~TiffIO();
	float** readFloatImage(std::string input_name, int* w_ptr, int* h_ptr, bool do_histogram);
	void writeFloatImage(float** image, std::string output_name, int width, int height);
 
 };
 
 #endif
