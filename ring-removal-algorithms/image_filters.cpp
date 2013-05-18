#include "image_filters.h"

using namespace std;

ImageFilterClass::ImageFiltersClass()
{
}

ImageFilterClass::~ImageFiltersClass()
{
}

void ImageFilterClass::DoMedianFilter1D(float *** filtered_image, float*** image, int start_row,
										int axis, int start_col, int kernel_rad, int filter_width,
										int width, int height)
{
	float* median_array = (float*) calloc(2*kernel_rad+1, sizeof(float));
	for(int n = -kernel_rad; n < kernel_rad + 1; n++){
		int row = start_row;
		int col = start_col + round(float(n)*float(ring_width)/float(2*kernel_rad));
		if(col < 0){
			col = -col;
			if(row < height/2){
				row += height/2;
			}else{
				row -= height/2;
			}
			median_array[n+kernel_rad] = image[0][row][col];
		}else if(col >= pol_width){
			median_array[n+kernel_rad] = 0.0;
		}else{
			median_array[n+kernel_rad] = image[0][row][col];
		}
	}
	Quicksort(median_array, 0, 2*kernel_rad);
	return median_array[kernel_rad];
}