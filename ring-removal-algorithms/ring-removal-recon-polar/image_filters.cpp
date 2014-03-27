#include "image_filters.h"

using namespace std;

ImageFilterClass::ImageFilterClass()
{
}

ImageFilterClass::~ImageFilterClass()
{
}

int ImageFilterClass::round(float x)
{
	if(x != 0.0){
		return floor(x+0.5);
	}else{
		return 0;
	}
}

void ImageFilterClass::elemSwap(float* arr, int index1, int index2)
{
	float store_value = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store_value;
}
void ImageFilterClass::elemSwapInteger(int* arr, int index1, int index2)
{
	int store_value = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store_value;
}

int ImageFilterClass::partition(float* median_array, int left, int right, int pivot_index)
{
	float pivot_value = median_array[pivot_index];
	elemSwap(median_array, pivot_index, right);
	int store_index = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivot_value){
			elemSwap(median_array, i, store_index);
			store_index +=1;
		}
	}
	elemSwap(median_array, store_index, right);
	return store_index;
}

int ImageFilterClass::partition2Arrays(float* median_array, int* position_array, int left, int right, int pivot_index)
{
	float pivot_value = median_array[pivot_index];
	elemSwap(median_array, pivot_index, right);
	elemSwapInteger(position_array, pivot_index, right);
	int store_index = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivot_value){
			elemSwap(median_array, i, store_index);
			elemSwapInteger(position_array, i, store_index);
			store_index +=1;
		}
	}
	elemSwap(median_array, store_index, right);
	elemSwapInteger(position_array, store_index, right);
	return store_index;
}

void ImageFilterClass::quickSort(float* median_array, int left, int right)
{
	if(left < right){
		int pivot_index = int((left + right)/2);
		int new_pivot_index = partition(median_array, left, right, pivot_index);
		quickSort(median_array, left, new_pivot_index - 1);
		quickSort(median_array, new_pivot_index + 1, right);
	}
}

void ImageFilterClass::quickSort2Arrays(float* median_array, int* position_array, int left, int right)
{
	if(left < right){
		int pivot_index = int((left + right)/2);
		int new_pivot_index = partition2Arrays(median_array, position_array, left, right, pivot_index);
		quickSort2Arrays(median_array, position_array, left, new_pivot_index - 1);
		quickSort2Arrays(median_array, position_array, new_pivot_index + 1, right);
	}
}
void ImageFilterClass::bubbleIntoPosition2Arrays(float* median_array, int* position_array, int index, int length)
{
	if(index > 0 && index < length -1){
		if(median_array[index] < median_array[index-1]){
			elemSwap(median_array, index, index-1);
			elemSwapInteger(position_array, index, index-1);
			bubbleIntoPosition2Arrays(median_array, position_array, index-1, length);
		}else if(median_array[index] > median_array[index+1]){
			elemSwap(median_array, index, index+1);
			elemSwapInteger(position_array, index, index+1);
			bubbleIntoPosition2Arrays(median_array, position_array, index+1, length);
		}
	}else if(index == 0){
		if(median_array[index] > median_array[index+1]){
			elemSwap(median_array, index, index+1);
			elemSwapInteger(position_array, index, index+1);
			bubbleIntoPosition2Arrays(median_array, position_array, index+1, length);
		}
	}else if(index == length -1){
		 if(median_array[index] < median_array[index-1]){
			elemSwap(median_array, index, index-1);
			elemSwapInteger(position_array, index, index-1);
			bubbleIntoPosition2Arrays(median_array, position_array, index-1, length);
		}
	}
}
void ImageFilterClass::doMedianFilterFast1D(float *** filtered_image, float*** image, int start_row,
										int start_col, int end_row, int end_col, char axis, 
										int kernel_rad,	int filter_width, int width, int height)
{
	int row, col;
	float* median_array = (float*) calloc(2*kernel_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*kernel_rad+1, sizeof(int));
	if(median_array == NULL || position_array == NULL){
		fprintf(stderr, "Error allocating memory for median filter arrays!");
	}
	if(axis == 'x'){
		for(row = start_row; row <= end_row; row++){
			col = start_col;
			for(int n = -kernel_rad; n < kernel_rad + 1; n++){
				int adjusted_col = col + n;
				int adjusted_row = row;
				if(adjusted_col < 0){
					adjusted_col = -adjusted_col;
					if(row < height/2){
						adjusted_row += height/2;
					}else{
						adjusted_row -= height/2;
					}
					median_array[n+kernel_rad] = image[0][adjusted_row][adjusted_col];
				}else{
					median_array[n+kernel_rad] = image[0][row][adjusted_col];
				}
				position_array[n+kernel_rad] = n+kernel_rad;
			}
			//Sort the array
			quickSort2Arrays(median_array, position_array, 0, 2*kernel_rad);
			filtered_image[0][row][col] = median_array[kernel_rad];
			
			//Roll filter along the rest of the row
			for(col = start_col+1; col <= end_col; col++){
				float next_value = 0.0;
				int next_value_col = col + kernel_rad;
				if (next_value_col < width){
					next_value = image[0][row][next_value_col];
				}
				float last_value;
				int last_value_index; 
				for(int i = 0; i < 2*kernel_rad+1; i++){
					position_array[i] -= 1;
					if(position_array[i] < 0){
						last_value_index = i;
						last_value = median_array[i];
						position_array[i] = 2*kernel_rad;
						median_array[i] = next_value;
					}
				}
				bubbleIntoPosition2Arrays(median_array, position_array, last_value_index, 2*kernel_rad+1);
				filtered_image[0][row][col] = median_array[kernel_rad];
			}
		}
	}else if(axis == 'y'){
		for(col = start_col; col <= end_col; col++){
			row = start_row;
			for(int n = -kernel_rad; n < kernel_rad + 1; n++){
				int adjusted_row = row + n;
				int adjusted_col = col;
				if(adjusted_row < 0){
				// Handle edge cases 
					adjusted_row += height;
					median_array[n+kernel_rad] = image[0][adjusted_row][adjusted_col];
				}else{
					median_array[n+kernel_rad] = image[0][adjusted_row][adjusted_col];
				}
				position_array[n+kernel_rad] = n+kernel_rad;
			}
			//Sort the array
			quickSort2Arrays(median_array, position_array, 0, 2*kernel_rad);
			filtered_image[0][row][col] = median_array[kernel_rad];
			
			//Roll filter along the rest of the col
			for(row = start_row+1; row <= end_row; row++){
				float next_value = 0.0;
				int next_value_row = row + kernel_rad;
				if (next_value_row < height){
					next_value = image[0][next_value_row][col];
				}
				float last_value;
				int last_value_index; 
				for(int i = 0; i < 2*kernel_rad+1; i++){
					position_array[i] -= 1;
					if(position_array[i] < 0){
						last_value_index = i;
						last_value = median_array[i];
						position_array[i] = 2*kernel_rad;
						median_array[i] = next_value;
					}
				}
				bubbleIntoPosition2Arrays(median_array, position_array, last_value_index, 2*kernel_rad+1);
				filtered_image[0][row][col] = median_array[kernel_rad];
			}
		}	
	}
	free(median_array);
	free(position_array);
}


void ImageFilterClass::doMedianFilter1D(float *** filtered_image, float*** image, int start_row,
										int start_col, int end_row, int end_col, char axis, 
										int kernel_rad,	int filter_width, int width, int height)
{
	int row, col;
	float* median_array = (float*) calloc(2*kernel_rad+1, sizeof(float));
	if(axis == 'x'){
		for(row = start_row; row <= end_row; row++){
			for(col = start_col; col <= end_col; col++){
				for(int n = -kernel_rad; n < kernel_rad + 1; n++){
					int subsampl_col = col + round(float(n)*float(filter_width)/float((2.0*kernel_rad)+1.0));
					int adjusted_row = row;
					if(subsampl_col < 0){
						subsampl_col = -subsampl_col;
						if(row < height/2){
							adjusted_row += height/2;
						}else{
							adjusted_row -= height/2;
						}
						median_array[n+kernel_rad] = image[0][adjusted_row][subsampl_col];
					}else if(subsampl_col >= width){
						median_array[n+kernel_rad] = 0.0;
					}else{
						median_array[n+kernel_rad] = image[0][row][subsampl_col];
					}
				}
				//Sort the array - this part is REALLY slow right now... why is that?
				quickSort(median_array, 0, 2*kernel_rad);
				filtered_image[0][row][col] = median_array[kernel_rad];
			}
		}
	}else if(axis == 'y'){
		for(col = start_col; col <= end_col; col++){
			for(row = start_row; row <= end_row; row++){
				for(int n = -kernel_rad; n < kernel_rad + 1; n++){
					int subsampl_row = row + round(float(n)*float(filter_width)/float(2*kernel_rad));
					//Dealing with edge cases - need to make this a bit more elegant...
					if(subsampl_row < 0){
						/*subsampl_col = -subsampl_col;
						if(row < height/2){
							row += height/2;
						}else{
							row -= height/2;
						}
						median_array[n+kernel_rad] = image[0][row][subsampl_col];
						*/
					}else if(subsampl_row= width){
						median_array[n+kernel_rad] = 0.0;
					}else{
						median_array[n+kernel_rad] = image[0][subsampl_row][col];
					}
				}
				//Sort the array - this part is REALLY slow right now... why is that?
				quickSort(median_array, 0, 2*kernel_rad);
				filtered_image[0][row][col] = median_array[kernel_rad];
			}
		}
	}
}

/* Runs slightly faster than the above mean filter, but floating-point rounding causes errors on
 * the order of 1E-10. Should be small enough error to not care about, but be careful...
 */
void ImageFilterClass::doMeanFilterFast1D(float*** filtered_image, float*** image,
 									      int start_row, int start_col, int end_row, int end_col,
										  char axis, int kernel_rad, int width, int height)
{
	float mean = 0, sum = 0, previous_sum = 0, num_elems = float(2*kernel_rad + 1);
	int row, col;
	if(axis == 'x'){
		//iterate over each row of the image subset
		for(row = start_row; row <= end_row; row++){
			sum = 0;
			//calculate average of first element of the column
			for(int n = - kernel_rad; n < (kernel_rad + 1); n++){
				col = n + start_col;
				if(col < 0){
					//col += height;
				}else if(col >= width){
					//col -= height;
				}
				sum += image[0][row][col];
			}
			mean = sum/num_elems;
			filtered_image[0][row][start_col] = mean;
			previous_sum = sum;

			for(col = start_col+1; col <= end_col; col++){
				int last_col = (col - 1) - (kernel_rad);
				int next_col = col + (kernel_rad);
				if(last_col < 0){
					last_col += width;
				}
				if(next_col >= width){
					next_col -= width;
				}
				sum = previous_sum - image[0][row][last_col] + image[0][row][next_col];
				filtered_image[0][row][col] = sum/num_elems;
				previous_sum = sum;
			}
		}
	}else if(axis == 'y'){
		//iterate over each column of the image subset
		for(col = start_col; col <= end_col; col++){
			sum = 0;
			//calculate average of first element of the column
			for(int n = - kernel_rad; n < (kernel_rad + 1); n++){
				row = n + start_row;
				if(row < 0){
					row += height;
				}else if(row >= height){
					row -= height;
				}
				sum += image[0][row][col];
			}
			mean = sum/num_elems;
			filtered_image[0][start_row][col] = mean;
			previous_sum = sum;

			for(row = start_row+1; row <= end_row; row++){
				int last_row = (row - 1) - (kernel_rad);
				int next_row = row + (kernel_rad);
				if(last_row < 0){
					last_row += height;
				}
				if(next_row >= height){
					next_row -= height;
				}
				sum = previous_sum - image[0][last_row][col] + image[0][next_row][col];
				if(image[0][row][col] != 0){
					filtered_image[0][row][col] = sum/num_elems;
				}else{
					filtered_image[0][row][col] = 0.0;
				}
				previous_sum = sum;
			}
		}
	}	
}

void ImageFilterClass::doMeanFilter1D(float*** filtered_image, float*** image, int start_row,
									  int start_col, int end_row, int end_col, char axis, 
									  int kernel_rad, int height, int width)
{
	float mean = 0, sum = 0, num_elems = float(2*kernel_rad +1);
	int col, row;
	if(axis == 'x'){
	
	}else if(axis == 'y'){
		for(col = start_col; col <= end_col; col++){
			for(int n = - kernel_rad; n < (kernel_rad + 1); n++){
				row = start_row + n;
				if(row < 0){
					row += height;
				}else if(row >= height){
					row -= height;
				}
				sum += image[0][row][col];
			}
			mean = sum/num_elems;
			filtered_image[0][row][col] =  mean;
		}
	}	
}
