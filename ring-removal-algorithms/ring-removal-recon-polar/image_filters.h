/* Class used to do various image filtering techniques, including 1-D median filters, 1-D mean filters,
 * etc.
 */

 #ifndef IMAGE_FILTERS_H
 #define IMAGE_FILTERS_H
 
 #include <cmath> 
 #include <cstdlib>
 #include <cstdio>
 #pragma once
 
 class ImageFilterClass
 {
 private:
	int round(float x);
	void elemSwap(float* arr, int index1, int index2);
	void elemSwapInteger(int* arr, int index1, int index2);
	int partition(float* median_array, int left, int right, int pivot_index);
	int partition2Arrays(float* median_array, int* position_array, int left, int right, int pivot_index);
	void quickSort(float* median_array, int left, int right);
	void quickSort2Arrays(float* median_array, int* position_array, int left, int right);
	void bubbleIntoPosition2Arrays(float* median_array, int* position_array, int index, int length);
	
 public:
	
	ImageFilterClass();
	~ImageFilterClass();
	void doMedianFilter1D(float *** filtered_image, float*** image, int start_row, int start_col,
						  int end_row, int end_col, char axis, int kernel_rad, int filter_width,
						  int width, int height);
	void doMedianFilterFast1D(float *** filtered_image, float*** image, int start_row, int start_col,
						  int end_row, int end_col, char axis, int kernel_rad, int filter_width,
						  int width, int height);
	void doMeanFilterFast1D(float*** filtered_image, float*** image, int start_row, int start_col,
							int end_row, int end_col, char axis, int kernel_rad, int width,
							int height);
	void doMeanFilter1D(float*** filtered_image, float*** image, int start_row, int start_col,
						int end_row, int end_col, char axis, int kernel_rad, int width, int height);
};

#endif
