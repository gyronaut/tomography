
#include <string>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include "time.h"
#include "tiff.h"
#include "tiffio.h"
using namespace std;


#define PI 3.14159265359

int round(float x){
	return (x > 0.0) ? floor(x+0.5) : ceil(x+0.5);
}

float** ReadFloatImage(string input_name, int* w_ptr, int* h_ptr)
{
	int width, height;

	TIFF* tif = TIFFOpen(input_name.c_str(), "r");

	if(!tif){
		return NULL;
	}

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	*w_ptr = width;
	*h_ptr = height;

	float* image = (float*) calloc(width*height, sizeof(float));
	float** image_rows = (float**) calloc(height, sizeof(float *));


	float *buffer;
	tstrip_t strip;
	uint32* bc;
	uint32 stripsize;
	uint16 num_strips= TIFFNumberOfStrips(tif);
	TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc);
	stripsize = bc[0];
	buffer = (float *)_TIFFmalloc(stripsize);
	uint32 location = 0;
	for (strip = 0; strip < num_strips; strip++) {
		if (bc[strip] != stripsize) {
			buffer = (float *)_TIFFrealloc(buffer, bc[strip]);
			stripsize = bc[strip];
		}
		TIFFReadEncodedStrip(tif, strip, buffer, bc[strip]);
		for(int i=0; i< bc[strip]/sizeof(float); i++){
			image[location+i] = buffer[i];
		}
		location += stripsize/sizeof(float);
	}
	_TIFFfree(buffer);
	TIFFClose(tif);

	image_rows[0] = image;
	for (int i=1; i<height; i++) {
		image_rows[i] = image_rows[i-1] + width;
	}
	return image_rows;
}

int FindMinDistanceToEdge(float center_x, float center_y, int width, int height){
	int* dist = new int[4];
	dist[0] = center_x;
	dist[1] = center_y;
	dist[2] = width - center_x;
	dist[3] = height - center_y;
	int min = dist[0];
	for(int i = 1; i < 4; i++){
		if(min > dist[i]){
			min = dist[i];
		}
	}
	return min;
}

float** PolarTransform(float** image, float center_x, float center_y, int width, int height, int* pol_width, int* pol_height, float thresh_max, float thresh_min){
	int max_r = FindMinDistanceToEdge(center_x, center_y, width, height);
	*pol_width = max_r;
	*pol_height = round(4.0*PI*float(max_r));
	
	float** polar_image = (float **) calloc(*pol_height, sizeof(float *));
	for(int i=0; i<*pol_height; i++){
		polar_image[i] = (float *) calloc(*pol_width, sizeof(float));
	}
	for(int row = 0; row<*pol_height; row++){
		for(int r = 0; r<center_x; r++){
			//float theta = float(row)/float(*pol_height)*3.0*PI - PI/2.0; //gives theta in the range [-PI/2, 5PI/2]
			float theta = float(row)/float(2.0*max_r);
			float fl_x = float(r)*cos(theta);
			float fl_y = float(r)*sin(theta);
			int x = round(fl_x + float(center_x));
			int y = round(fl_y + float(center_y));
			
			polar_image[row][r] = image[y][x];
			if(polar_image[row][r] > thresh_max){
				polar_image[row][r] = thresh_max;
			}else if(polar_image[row][r] < thresh_min){
				polar_image[row][r] = thresh_min;
			}
		}
	}

	return polar_image;
}

/*Polar transform that uses simple bilinear interpolation to translate from (x,y) to (r, theta)*/
float** PolarTransformBilinear(float** image, float center_x, float center_y, int width, int height, int* pol_width, int* pol_height, float thresh_max, float thresh_min){
	int max_r = FindMinDistanceToEdge(center_x, center_y, width, height);
	*pol_width = max_r;
	*pol_height = round(4.0*PI*float(max_r));
	
	float** polar_image = (float **) calloc(*pol_height, sizeof(float *));
	for(int i=0; i<*pol_height; i++){
		polar_image[i] = (float *) calloc(*pol_width, sizeof(float));
	}

	float theta, fl_x, fl_y, value;
	int x_1, x_2, y_1, y_2;
	for(int row = 0; row<*pol_height; row++){
		for(int r = 0; r<*pol_width; r++){
			//float theta = float(row)/float(*pol_height)*3.0*PI - PI/2.0; //gives theta in the range [-PI/2, 5PI/2]
			theta = float(row)/float(2.0*max_r);
			fl_x = float(r)*cos(theta) + float(center_x);
			fl_y = float(r)*sin(theta) + float(center_y);
			x_1 = floor(fl_x);
			x_2 = ceil(fl_x);
			y_1 = floor(fl_y);
			y_2 = ceil(fl_y);
			if(x_1 == x_2){
				if(y_1 ==y_2){
					value = image[y_1][x_1];
				}else{
					value = (y_2 - fl_y)*image[y_1][x_1] + (fl_y - y_1)*image[y_2][x_1];
				}
			}else if(y_1 == y_2){
				value = (x_2 - fl_x)*image[y_1][x_1] + (fl_x - x_1)*image[y_1][x_2];
			}else{
				value = (y_2 - fl_y)*((x_2 - fl_x)*image[y_1][x_1] + (fl_x - x_1)*image[y_1][x_2]) + (fl_y - y_1)*((x_2 - fl_x)*image[y_2][x_1] + (fl_x - x_1)*image[y_2][x_2]);
			}
			polar_image[row][r] = value;
			if(polar_image[row][r] > thresh_max){
				polar_image[row][r] = thresh_max;
			}else if(polar_image[row][r] < thresh_min){
				polar_image[row][r] = thresh_min;
			}
		}
	}

	return polar_image;
}

float** InversePolarTransform(float** polar_image, float center_x, float center_y, int pol_width, int  pol_height, int width, int height){
	float** cart_image = (float **) calloc(height, sizeof(float *));
	for(int i = 0; i < height; i++){
		cart_image[i] = (float *) calloc(width, sizeof(float));
	}
	for(int row = 0; row < pol_height; row++){
		for(int r = 0; r < pol_width; r++){
			//float theta = float(row)/float(pol_height)*3.0*PI - PI/2.0;
			float theta = float(row)/float(2.0*pol_width); 
			int x = round(float(r)*cos(theta) + float(center_x));
			int y = round(float(r)*sin(theta) + float(center_y));
			cart_image[y][x] = polar_image[row][r];
		}
	}
	return cart_image;
}


/*Inverse polar transform that uses simple bilinear interpolation to translate from (r, theta) to (x, y)*/
float** InversePolarTransformBilinear(float** polar_image, float center_x, float center_y, int pol_width, int  pol_height, int width, int height){
	float** cart_image = (float **) calloc(height, sizeof(float *));
	for(int i = 0; i < height; i++){
		cart_image[i] = (float *) calloc(width, sizeof(float));
	}
	float r, theta, value, row_float;
	int r_1, r_2, row_1, row_2;
	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++){
			value = 0;
			r = sqrt(float((x - center_x)*(x - center_x)) + float((y - center_y)*(y - center_y)));
			if(r < pol_width - 1){
				r_1 = floor(r);
				r_2 = ceil(r);
				theta = atan2(float(y - center_y), float(x - center_x)) + PI;
				row_float = theta*float(pol_height -1)/(2.0*PI);
				if(row_float > pol_height){
					row_float -= pol_height;
				}
				row_1 = floor(row_float);
				row_2 = ceil(row_float);
				if(r_1 == r_2){
					if(row_1 == row_2){
						value = polar_image[row_1][r_1];
					}else{
						value = (row_2 - row_float)*polar_image[row_1][r_1] + (row_float - row_1)*polar_image[row_2][r_1];
					}
				}else if(row_1 == row_2){
						value = (r_2 - r)*polar_image[row_1][r_1] + (r - r_1)*polar_image[row_1][r_2];
				}else{
					value = (row_2 - row_float)*((r_2 - r)*polar_image[row_1][r_1] + (r - r_1)*polar_image[row_1][r_2]) + (row_float - row_1)*((r_2 - r)*polar_image[row_2][r_1] + (r - r_1)*polar_image[row_2][r_2]);
				}
			}
			cart_image[(height-1)-y][(width-1)-x] = value;
		}
	}
	return cart_image;
}

void ElemSwapInt(int* arr, int index1, int index2){
	int store1 = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store1;
}
void ElemSwap(float* arr, int index1, int index2){
	float store1 = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store1;
}

/*Partition and Quicksort algorithm taken from Wikipedia Quicksort Page:
 *
 * http://en.wikipedia.org/wiki/Quicksort#In-place_version
 *
 */
int Partition(float* median_array, int left, int right, int pivot_index){
	float pivot_value = median_array[pivot_index];
	ElemSwap(median_array, pivot_index, right);
	int store_index = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivot_value){
			ElemSwap(median_array, i, store_index);
			store_index +=1;
		}
	}
	ElemSwap(median_array, store_index, right);
	return store_index;
}

void Quicksort(float* median_array, int left, int right){
	if(left < right){
		int pivot_index = int((left + right)/2);
		int new_pivot_index = Partition(median_array, left, right, pivot_index);
		Quicksort(median_array, left, new_pivot_index - 1);
		Quicksort(median_array, new_pivot_index + 1, right);
	}
}

int partition_2(float* median_array, int* position_array, int left, int right, int pivot_index){
	float pivot_value = median_array[pivot_index];
	ElemSwap(median_array, pivot_index, right);
	ElemSwapInt(position_array, pivot_index, right);
	int store_index = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivot_value){
			ElemSwap(median_array, i, store_index);
			ElemSwapInt(position_array, i, store_index);
			store_index +=1;
		}
	}
	ElemSwap(median_array, store_index, right);
	ElemSwapInt(position_array, store_index, right);
	return store_index;
}

void quicksort_2(float* median_array, int* position_array, int left, int right){
	if(left < right){
		int pivot_index = int((left + right)/2);
		int new_pivot_index = partition_2(median_array, position_array, left, right, pivot_index);
		Quicksort(median_array, left, new_pivot_index - 1);
		Quicksort(median_array, new_pivot_index + 1, right);
	}
}

float DoRadialMedianFilter(float*** polar_image, int start_row, int start_col, int m_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*m_rad+1, sizeof(float));
	for(int n = -m_rad; n < m_rad + 1; n++){
		int row = start_row;
		int col = start_col + round(float(n)*float(ring_width)/float(2*m_rad));
		if(col < 0){
			col = -col;
			if(row < pol_height/2){
				row += pol_height/2;
			}else{
				row -= pol_height/2;
			}
			median_array[n+m_rad] = polar_image[0][row][col];
		}else if(col >= pol_width){
			median_array[n+m_rad] = 0.0;
		}else{
			median_array[n+m_rad] = polar_image[0][row][col];
		}
	}
	Quicksort(median_array, 0, 2*m_rad);
	return median_array[m_rad];
}

void radial_median_filter_faster_inner(float*** polar_image, float*** filtered_image, int start_row, int m_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*m_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*m_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -m_rad; n < m_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*m_rad));
		position_array[n+m_rad] = col + ring_width/2;
		if(col < 0){
			col = -col;
			if(row < pol_height/2){
				row += pol_height/2;
			}else{
				row -= pol_height/2;
			}
			median_array[n+m_rad] = polar_image[0][row][col];
		}else{
			median_array[n+m_rad] = polar_image[0][row][col];
		}
	}
	quicksort_2(median_array, position_array, 0, 2*m_rad);
	filtered_image[0][start_row][0] = median_array[m_rad];
	for(int col = 1; col < pol_width/3; col++){
		for(int i = 0; i < 2*m_rad + 1; i++){
			position_array[i] -= 1;
			if(position_array[i] < 0){
				index_last_elem = i;
				position_array[i] = ring_width/2;
				median_array[i] = polar_image[0][start_row][col + ring_width/2];
			}
		}
		partition_2(median_array, position_array, 0, 2*m_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[m_rad];
	}
}

void radial_median_filter_faster_middle(float*** polar_image, float*** filtered_image, int start_row, int m_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*m_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*m_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -m_rad; n < m_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*m_rad)) + pol_width/3;
		position_array[n+m_rad] = col + ring_width/2 - pol_width/3;
		median_array[n+m_rad] = polar_image[0][row][col];
	}
	quicksort_2(median_array, position_array, 0, 2*m_rad);
	filtered_image[0][start_row][pol_width/3] = median_array[m_rad];
	for(int col = pol_width/3 + 1; col < 2*pol_width/3; col++){
		for(int i = 0; i < 2*m_rad + 1; i++){
			position_array[i] -= 1;
			if(position_array[i] < 0){
				index_last_elem = i;
				position_array[i] = ring_width/2;
				median_array[i] = polar_image[0][start_row][col + ring_width/2];
			}
		}
		partition_2(median_array, position_array, 0, 2*m_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[m_rad];
	}
}

void radial_median_filter_faster_outer(float*** polar_image, float*** filtered_image, int start_row, int m_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*m_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*m_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -m_rad; n < m_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*m_rad)) + 2*pol_width/3;
		position_array[n+m_rad] = col + ring_width/2 - 2*pol_width/3;
		if(col >= pol_width){
			median_array[n+m_rad] = 0.0;
		}else{
			median_array[n+m_rad] = polar_image[0][row][col];
		}
	}
	quicksort_2(median_array, position_array, 0, 2*m_rad);
	filtered_image[0][start_row][2*pol_width/3] = median_array[m_rad];
	for(int col = 2*pol_width/3 + 1; col < pol_width; col++){
		for(int i = 0; i < 2*m_rad + 1; i++){
			position_array[i] -= 1;
			if(position_array[i] < 0){
				index_last_elem = i;
				position_array[i] = ring_width/2;
				int col_2 = col + ring_width/2;
				if( col_2 >= pol_width){
					median_array[i] = 0.0;
				}else{
					median_array[i] = polar_image[0][start_row][col_2];
				}
			}
		}
		partition_2(median_array, position_array, 0, 2*m_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[m_rad];
	}
}

float azi_mean_filter(float*** polar_image, int start_row, int start_col, int m_azi, int pol_height){
	float mean = 0, sum = 0, num_elems = float(2*m_azi +1);
	int col = start_col, row;
	for(int n = - m_azi; n < (m_azi + 1); n++){
		row = start_row + n;
		if(row < 0){
			row += pol_height;
		}else if(row >= pol_height){
			row -= pol_height;
		}
		sum += polar_image[0][row][col];
	}
	mean = sum/num_elems;
	return mean;
}

//Runs slightly faster than the above mean filter, but floating-point rounding causes errors on the order
//of 1E-10. Should be small enough error to not care about, but be careful...
void DoAziMeanFilterFast(float*** mean_filtered_image, float*** polar_image, int col, int m_azi, int pol_height){
	float mean = 0, sum = 0, previous_sum = 0, num_elems = float(2*m_azi + 1);
	int row;
	//calculate average of first element of the column
	for(int n = - m_azi; n < (m_azi + 1); n++){
		row = n;
		if(row < 0){
			row += pol_height;
		}else if(row >= pol_height){
			row -= pol_height;
		}
		sum += polar_image[0][row][col];
	}
	mean = sum/num_elems;
	mean_filtered_image[0][0][col] = mean;
	previous_sum = sum;

	for(row = 1; row < pol_height; row++){
		int last_row = (row - 1) - (m_azi);
		int next_row = row + (m_azi);
		if(last_row < 0){
			last_row += pol_height;
		}
		if(next_row >= pol_height){
			next_row -= pol_height;
		}
		sum = previous_sum - polar_image[0][last_row][col] + polar_image[0][next_row][col];
		mean_filtered_image[0][row][col] = sum/num_elems;
		previous_sum = sum;
	}
}

void DoRingFilter(float*** polar_image, int pol_height, int pol_width, float threshold, int m_rad, int m_azi, int ring_width){
	float** filtered_image = (float **) calloc(pol_height, sizeof(float *));
	float** mean_filtered_image = (float **) calloc(pol_height, sizeof(float *));
	for(int i=0; i<pol_height; i++){
		filtered_image[i] = (float *) calloc(pol_width, sizeof(float));
		mean_filtered_image[i] = (float *) calloc(pol_width, sizeof(float));
	}
	
	//Do radial median filter to get filtered_image
	printf("Performing Radial Filter on polar image... \n");
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col <  pol_width; col++){
			if(col < pol_width/3){
				filtered_image[row][col] = DoRadialMedianFilter(polar_image, row, col, m_rad, ring_width, pol_width, pol_height);
			}else if(col < 2*pol_width/3){
				filtered_image[row][col] = DoRadialMedianFilter(polar_image, row, col, 2*m_rad/3, ring_width, pol_width, pol_height);
			}else{
				filtered_image[row][col] = DoRadialMedianFilter(polar_image, row, col, m_rad/3, ring_width, pol_width, pol_height);
			}
		}
	}
	/*
	printf("Performing Faster Radial Filter on polar_image... \n");
	for(int row = 0; row < pol_height; row++){
		radial_median_filter_faster_inner(polar_image, &filtered_image, row, m_rad, ring_width, pol_width, pol_height);
		radial_median_filter_faster_middle(polar_image, &filtered_image, row, 2*m_rad/3, ring_width, pol_width, pol_height);
		radial_median_filter_faster_outer(polar_image, &filtered_image, row, m_rad/3, ring_width, pol_width, pol_height);
	}
	*/
	//subtract filtered image from polar image to get difference image & do last thresholding
	printf("Calculating Difference Image... \n");
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col <  pol_width; col++){
			polar_image[0][row][col] -= filtered_image[row][col];
			if(polar_image[0][row][col] > threshold || polar_image[0][row][col] < -threshold){
				polar_image[0][row][col] = 0;
			}
			//polar_image[0][row][col] = filtered_image[row][col];
		}
	}

	/*
	//Do Azimuthal mean filter
	printf("Performing Azimuthal Filter on polar image... \n");
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col < pol_width; col++){
			if(col < pol_width/3){
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, m_azi, pol_height);
			}else if(col < 2*pol_width/3){
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, 2*m_azi/3, pol_height);
			}else{
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, m_azi/3, pol_height);
			}
		}
	}
	*/
	
	
	//Do Azimuthal filter #2 (faster mean, does whole column in one call)
	for(int col = 0; col < pol_width; col++){
		if(col < pol_width/3){
			DoAziMeanFilterFast(&mean_filtered_image, polar_image, col, m_azi, pol_height);
		}else if(col < 2*pol_width/3){
			DoAziMeanFilterFast(&mean_filtered_image, polar_image, col, 2*m_azi/3, pol_height);
		}else{
			DoAziMeanFilterFast(&mean_filtered_image, polar_image, col, m_azi/3, pol_height);
		}
	}
	

	//Set "polar_image" to the fully filtered data
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col < pol_width; col++){
			polar_image[0][row][col] = mean_filtered_image[row][col];
		}
	}
}

void WriteFloat(float** corrected_image, string output_name, int width, int height)
{
	float* output = (float *) calloc(width*height, sizeof(float));
	float* buffer = (float *) calloc(width, sizeof(float));
	for(int i = 0; i < height; i++){
		memcpy(buffer, corrected_image[i], width*sizeof(float));
		for(int j = 0; j < width; j++){
			output[j+(i*width)] = buffer[j]; 
		}
	}
	TIFF* resultTif = TIFFOpen(output_name.c_str(), "w");
	if(!resultTif){
		return;
	}
	TIFFSetField(resultTif, TIFFTAG_IMAGEWIDTH, width);  
	TIFFSetField(resultTif, TIFFTAG_IMAGELENGTH, height);    
	TIFFSetField(resultTif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(resultTif, TIFFTAG_SAMPLEFORMAT, 3); //floating point = 3
	TIFFSetField(resultTif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(resultTif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(resultTif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(resultTif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	float *buf = NULL;
	tsize_t linebytes =  width*sizeof(float);
	buf =(float *)_TIFFmalloc(linebytes);
	
	for(uint32 row = 0; row < height; row++){
		for(int col = 0; col < width; col++){
			buf[col] = output[col+row*width];
		}
		if (TIFFWriteScanline(resultTif, buf, row, 0) < 0){
			break;
		}
	}
	TIFFClose(resultTif);
	if(buf){
		_TIFFfree(buf);
	}
	free(buffer);
	free(output);
}

string getName(string name_base, int img_num){
	stringstream stream;
	stream << name_base << img_num << ".tif";
	return stream.str();
}

int main(int argc, char** argv){
	if(!(argc == 12)){
		printf("\nUsage:\n\nring_remover_recon [input path] [output path] [input root] [output root] [first file num] [last file num] [center x y] [max ring width] [ring threshold] [verbose]\n\n");
		printf("      [input path]    Path to the folder containing the input images.\n");
		printf("     [output path]    Path to the folder to hold the filtered images.\n");
		printf("   [sinogram root]    Root name of the image to filter, with no file\n");
		printf("                      number or extension.\n");
		printf("[output file root]    Root of the output file name (number and file extension\n");
		printf("                      will be appended to this).\n");
		printf("  [first file num]    Number of first sinogram to reconstruct.\n");
		printf("   [last file num]    Number of last sinogram to reconstruct.\n");
		printf("    [center: x, y]    X and Y value of the center of rotation\n");
		printf("  [max ring width]    maximum width of the rings to be filtered in pixels.\n");
		printf("  [ring threshold]    Rings are treated as offsets to the corrected image. This\n");
		printf("                      threshold is the maximum value of an offset due to a\n");
		printf("                      ring artifact.\n");
		printf("         [verbose]    0 = No Output messages (default), 1 = Output messages\n");
		return 0;
	}else{
		int first_img_num;
		int last_img_num;
		int verbose;
		int num_files;
		int width=0;
		int height=0;
		int pol_width=0;
		int pol_height=0;
		int m_rad = 30;
		int m_azi = 60;
		int ring_width = 25;
		string input_base, input_name, input_path;
		string output_base, output_name, output_path;
		float center_x=511, center_y=511, thresh_max=0.0015, thresh_min=0.00005, threshold = 0.00034;
		char * filter;
		input_path.assign(argv[1], find(argv[1], argv[1]+255, '\0'));
		output_path.assign(argv[2], find(argv[2], argv[2]+255, '\0'));
		input_base.assign(argv[3], find(argv[3], argv[3]+255, '\0'));
		output_base.assign(argv[4], find(argv[4], argv[4]+255, '\0'));
		if(input_path.substr(input_path.length() - 1) != "/" && input_path.substr(input_path.length() - 1) != "\\"){
			input_path.append("/");
		}
		if(output_path.substr(output_path.length() - 1) != "/" && output_path.substr(output_path.length() - 1) != "\\"){
			output_path.append("/");
		}
		first_img_num = atoi(argv[5]);
		last_img_num = atoi(argv[6]);
		center_x = atof(argv[7]);
		center_y = atof(argv[8]);
		ring_width = atoi(argv[9]);
		threshold = atof(argv[10]);
		verbose = atoi(argv[11]);

		for(int img = first_img_num; img < last_img_num + 1; img++){

			time_t t_start, t_end;
			float **image=0, **polar_image=0, **ring_image=0;
			input_name = getName(input_base, img);
			output_name = getName(output_base, img);

			//read in image, perform thresholding

			image = ReadFloatImage(input_path + input_name, &width, &height);    //pixel data is stored in the form image[row][col]

			if(!image){
				return 1;
			}
			//Translate Image to Polar Coordinates

			polar_image = PolarTransformBilinear(image, center_x, center_y, width, height, &pol_width, &pol_height, thresh_max, thresh_min);

			//Call Ring Algorithm

			time(&t_start);
			DoRingFilter(&polar_image, pol_height, pol_width, threshold, m_rad, m_azi, ring_width);
			time(&t_end);
			double seconds = difftime(t_end, t_start);
			printf("Time to perform Ring Filtering: %.2f sec. \n", seconds);

			//Translate Ring-Image to Cartesian Coordinates

			ring_image = InversePolarTransformBilinear(polar_image, center_x, center_y, pol_width, pol_height, width, height);

			//Subtract Ring-Image from Image
	
			for(int row=0; row < height; row++){
				for(int col=0; col < width; col++){
					image[row][col] -= ring_image[row][col];
				}
			}
	
			//Write out Corrected-Image
			WriteFloat(image, output_path + output_name, width, height);

		}
		return 0;
	}
}
