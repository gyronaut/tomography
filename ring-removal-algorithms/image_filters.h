/* Class used to do various image filtering techniques, including 1-D median filters, 1-D mean filters,
 * etc.
 */

 #ifndef IMAGE_FILTERS_H
 #define IMAGE_FILTERS_H
 
 #pragma once
 
 class ImageFilterClass
 {
 private:
 
	void ElemSwap(float* arr, int index1, int index2);
	int Partition(float* median_array, int left, int right, int pivot_index);
	void Quicksort(float* median_array, int left, int right);
	
 public:
	
	ImageFilterClass();
	~ImageFilterClass();
	void DoMedianFilter1D(float *** filtered_image, float*** image, int start_row, int start_col,
						  int axis, int kernel_rad, int filter_width, int width, int height);
	void DoMeanFilterFast1D(float*** filtered_image, float*** image, int start_row, int start_col,
							int axis, int kernel_rad, int width, int height);
	void DoMeanFilter1D(float*** filtered_image, float*** image, int start_row, int start_col,
						int axis, int kernel_rad, int width, int height);
};
#endif