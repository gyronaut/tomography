/* Class to read and write tiff files */

#ifndef IMAGE_FILTERS_H
 #define IMAGE_FILTERS_H
 
 #pragma once
 
 class TiffIO
 {
 private:
	
 public:
	TiffIO();
	~TiffIO();
	float** ReadFloatImage(string input_name, int* w_ptr, int* h_ptr);
	void WriteFloatImage(float** image, string output_name, int width, int height);
 
 };
 #endif