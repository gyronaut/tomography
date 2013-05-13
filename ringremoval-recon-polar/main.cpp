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

float** read_float_image(string inputName, int* w_ptr, int* h_ptr)
{
	int width, height;

	TIFF* tif = TIFFOpen(inputName.c_str(), "r");

	if(!tif){
		return NULL;
	}

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	*w_ptr = width;
	*h_ptr = height;

	float* img = (float*) calloc(width*height, sizeof(float));
	float** image = (float**) calloc(height, sizeof(float *));


	float *buffer;
	tstrip_t strip;
	uint32* bc;
	uint32 stripsize;
	uint16 numStrips= TIFFNumberOfStrips(tif);
	TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc);
	stripsize = bc[0];
	buffer = (float *)_TIFFmalloc(stripsize);
	uint32 location = 0;
	for (strip = 0; strip < numStrips; strip++) {
		if (bc[strip] != stripsize) {
			buffer = (float *)_TIFFrealloc(buffer, bc[strip]);
			stripsize = bc[strip];
		}
		TIFFReadEncodedStrip(tif, strip, buffer, bc[strip]);
		for(int i=0; i< bc[strip]/sizeof(float); i++){
			img[location+i] = buffer[i];
		}
		location += stripsize/sizeof(float);
	}
	_TIFFfree(buffer);
	TIFFClose(tif);

	image[0] = img;
	for (int i=1; i<height; i++) {
		image[i] = image[i-1] + width;
	}
	return image;
}

int find_min_distance_to_edge(float center_x, float center_y, int width, int height){
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

/*Polar transform that uses simple bilinear interpolation to translate from (x,y) to (r, theta)*/
float** polar_transform_bilinear(float** image, float center_x, float center_y, int width, int height, int* pol_width, int* pol_height, float thresh_max, float thresh_min){
	int max_r = find_min_distance_to_edge(center_x, center_y, width, height);
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


/*Inverse polar transform that uses simple bilinear interpolation to translate from (r, theta) to (x, y)*/
float** inverse_polar_transform_bilinear(float** polar_image, float center_x, float center_y, int pol_width, int  pol_height, int width, int height){
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

void elem_swap_int(int* arr, int index1, int index2){
	int store1 = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store1;
}
void elem_swap(float* arr, int index1, int index2){
	float store1 = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store1;
}

/*Partition and Quicksort algorithm taken from Wikipedia Quicksort Page:
 *
 * http://en.wikipedia.org/wiki/Quicksort#In-place_version
 *
 */
int partition(float* median_array, int left, int right, int pivotIndex){
	float pivotValue = median_array[pivotIndex];
	elem_swap(median_array, pivotIndex, right);
	int storeIndex = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivotValue){
			elem_swap(median_array, i, storeIndex);
			storeIndex +=1;
		}
	}
	elem_swap(median_array, storeIndex, right);
	return storeIndex;
}

void quicksort(float* median_array, int left, int right){
	if(left < right){
		int pivotIndex = int((left + right)/2);
		int newPivotIndex = partition(median_array, left, right, pivotIndex);
		quicksort(median_array, left, newPivotIndex - 1);
		quicksort(median_array, newPivotIndex + 1, right);
	}
}

int partition_2(float* median_array, int* position_array, int left, int right, int pivotIndex){
	float pivotValue = median_array[pivotIndex];
	elem_swap(median_array, pivotIndex, right);
	elem_swap_int(position_array, pivotIndex, right);
	int storeIndex = left;
	for(int i=left; i < right; i++){
		if(median_array[i] <= pivotValue){
			elem_swap(median_array, i, storeIndex);
			elem_swap_int(position_array, i, storeIndex);
			storeIndex +=1;
		}
	}
	elem_swap(median_array, storeIndex, right);
	elem_swap_int(position_array, storeIndex, right);
	return storeIndex;
}

void quicksort_2(float* median_array, int* position_array, int left, int right){
	if(left < right){
		int pivotIndex = int((left + right)/2);
		int newPivotIndex = partition_2(median_array, position_array, left, right, pivotIndex);
		quicksort(median_array, left, newPivotIndex - 1);
		quicksort(median_array, newPivotIndex + 1, right);
	}
}

float radial_median_filter(float*** polar_image, int start_row, int start_col, int M_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*M_rad+1, sizeof(float));
	for(int n = -M_rad; n < M_rad + 1; n++){
		int row = start_row;
		int col = start_col + round(float(n)*float(ring_width)/float(2*M_rad));
		if(col < 0){
			col = -col;
			if(row < pol_height/2){
				row += pol_height/2;
			}else{
				row -= pol_height/2;
			}
			median_array[n+M_rad] = polar_image[0][row][col];
		}else if(col >= pol_width){
			median_array[n+M_rad] = 0.0;
		}else{
			median_array[n+M_rad] = polar_image[0][row][col];
		}
	}
	quicksort(median_array, 0, 2*M_rad);
	return median_array[M_rad];
}

void radial_median_filter_faster_inner(float*** polar_image, float*** filtered_image, int start_row, int M_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*M_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*M_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -M_rad; n < M_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*M_rad));
		position_array[n+M_rad] = col + ring_width/2;
		if(col < 0){
			col = -col;
			if(row < pol_height/2){
				row += pol_height/2;
			}else{
				row -= pol_height/2;
			}
			median_array[n+M_rad] = polar_image[0][row][col];
		}else{
			median_array[n+M_rad] = polar_image[0][row][col];
		}
	}
	quicksort_2(median_array, position_array, 0, 2*M_rad);
	filtered_image[0][start_row][0] = median_array[M_rad];
	for(int col = 1; col < pol_width/3; col++){
		for(int i = 0; i < 2*M_rad + 1; i++){
			position_array[i] -= 1;
			if(position_array[i] < 0){
				index_last_elem = i;
				position_array[i] = ring_width/2;
				median_array[i] = polar_image[0][start_row][col + ring_width/2];
			}
		}
		partition_2(median_array, position_array, 0, 2*M_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[M_rad];
	}
}

void radial_median_filter_faster_middle(float*** polar_image, float*** filtered_image, int start_row, int M_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*M_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*M_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -M_rad; n < M_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*M_rad)) + pol_width/3;
		position_array[n+M_rad] = col + ring_width/2 - pol_width/3;
		median_array[n+M_rad] = polar_image[0][row][col];
	}
	quicksort_2(median_array, position_array, 0, 2*M_rad);
	filtered_image[0][start_row][pol_width/3] = median_array[M_rad];
	for(int col = pol_width/3 + 1; col < 2*pol_width/3; col++){
		for(int i = 0; i < 2*M_rad + 1; i++){
			position_array[i] -= 1;
			if(position_array[i] < 0){
				index_last_elem = i;
				position_array[i] = ring_width/2;
				median_array[i] = polar_image[0][start_row][col + ring_width/2];
			}
		}
		partition_2(median_array, position_array, 0, 2*M_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[M_rad];
	}
}

void radial_median_filter_faster_outer(float*** polar_image, float*** filtered_image, int start_row, int M_rad, int ring_width, int pol_width, int pol_height){
	float* median_array = (float*) calloc(2*M_rad+1, sizeof(float));
	int* position_array = (int*) calloc(2*M_rad+1, sizeof(int));
	int index_last_elem;
	for(int n = -M_rad; n < M_rad + 1; n++){
		int row = start_row;
		int col = round(float(n)*float(ring_width)/float(2*M_rad)) + 2*pol_width/3;
		position_array[n+M_rad] = col + ring_width/2 - 2*pol_width/3;
		if(col >= pol_width){
			median_array[n+M_rad] = 0.0;
		}else{
			median_array[n+M_rad] = polar_image[0][row][col];
		}
	}
	quicksort_2(median_array, position_array, 0, 2*M_rad);
	filtered_image[0][start_row][2*pol_width/3] = median_array[M_rad];
	for(int col = 2*pol_width/3 + 1; col < pol_width; col++){
		for(int i = 0; i < 2*M_rad + 1; i++){
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
		partition_2(median_array, position_array, 0, 2*M_rad, index_last_elem);
		filtered_image[0][start_row][col] = median_array[M_rad];
	}
}

float azi_mean_filter(float*** polar_image, int start_row, int start_col, int M_azi, int pol_height){
	float mean = 0, sum = 0, num_elems = float(2*M_azi +1);
	int col = start_col, row;
	for(int n = - M_azi; n < (M_azi + 1); n++){
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
void azi_mean_filter_fast(float*** mean_filtered_image, float*** polar_image, int col, int M_azi, int pol_height){
	float mean = 0, sum = 0, previous_sum = 0, num_elems = float(2*M_azi + 1);
	int row;
	//calculate average of first element of the column
	for(int n = - M_azi; n < (M_azi + 1); n++){
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
		int last_row = (row - 1) - (M_azi);
		int next_row = row + (M_azi);
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

void ring_filter(float*** polar_image, int pol_height, int pol_width, float threshold, int M_rad, int M_azi, int ring_width){
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
				filtered_image[row][col] = radial_median_filter(polar_image, row, col, M_rad, ring_width, pol_width, pol_height);
			}else if(col < 2*pol_width/3){
				filtered_image[row][col] = radial_median_filter(polar_image, row, col, 2*M_rad/3, ring_width, pol_width, pol_height);
			}else{
				filtered_image[row][col] = radial_median_filter(polar_image, row, col, M_rad/3, ring_width, pol_width, pol_height);
			}
		}
	}
	/*
	printf("Performing Faster Radial Filter on polar_image... \n");
	for(int row = 0; row < pol_height; row++){
		radial_median_filter_faster_inner(polar_image, &filtered_image, row, M_rad, ring_width, pol_width, pol_height);
		radial_median_filter_faster_middle(polar_image, &filtered_image, row, 2*M_rad/3, ring_width, pol_width, pol_height);
		radial_median_filter_faster_outer(polar_image, &filtered_image, row, M_rad/3, ring_width, pol_width, pol_height);
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
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, M_azi, pol_height);
			}else if(col < 2*pol_width/3){
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, 2*M_azi/3, pol_height);
			}else{
				mean_filtered_image[row][col] = azi_mean_filter(polar_image, row, col, M_azi/3, pol_height);
			}
		}
	}
	*/
	
	
	//Do faster Azimuthal filter (faster mean, does whole column in one call)
	for(int col = 0; col < pol_width; col++){
		if(col < pol_width/3){
			azi_mean_filter_fast(&mean_filtered_image, polar_image, col, M_azi, pol_height);
		}else if(col < 2*pol_width/3){
			azi_mean_filter_fast(&mean_filtered_image, polar_image, col, 2*M_azi/3, pol_height);
		}else{
			azi_mean_filter_fast(&mean_filtered_image, polar_image, col, M_azi/3, pol_height);
		}
	}
	

	//Set "polar_image" to the fully filtered data
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col < pol_width; col++){
			polar_image[0][row][col] = mean_filtered_image[row][col];
		}
	}
}

void write_float(float** corrected_image, string outputName, int width, int height)
{
	float* output = (float *) calloc(width*height, sizeof(float));
	float* buffer = (float *) calloc(width, sizeof(float));
	for(int i = 0; i < height; i++){
		memcpy(buffer, corrected_image[i], width*sizeof(float));
		for(int j = 0; j < width; j++){
			output[j+(i*width)] = buffer[j]; 
		}
	}
	TIFF* resultTif = TIFFOpen(outputName.c_str(), "w");
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

string getName(string nameBase, int imgNum){
	stringstream stream;
	string number;
	stream << nameBase << imgNum << ".tif";
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
		printf("                      ring artifact./n");
		printf("         [verbose]    0 = No Output Messages (default), 1 = Output Messages\n");
		return 0;
	}else{
		int firstImgNum, lastImgNum, verbose, numFiles, numZerosEachSide=0, datatype;
		int width=0, height=0, pol_width=0, pol_height=0, M_rad = 30, M_azi = 60, ring_width = 25;
		string inputName, inputBase, outputBase, outputName, inputPath, outputPath;
		float center_x=511, center_y=511, thresh_max=0.0015, thresh_min=0.00005, threshold = 0.00034;
		char * filter;
		inputPath.assign(argv[1], find(argv[1], argv[1]+255, '\0'));
		outputPath.assign(argv[2], find(argv[2], argv[2]+255, '\0'));
		inputBase.assign(argv[3], find(argv[3], argv[3]+255, '\0'));
		outputBase.assign(argv[4], find(argv[4], argv[4]+255, '\0'));
		if(inputPath.substr(inputPath.length() - 1) != "/" && inputPath.substr(inputPath.length() - 1) != "\\"){
			inputPath.append("/");
		}
		if(outputPath.substr(outputPath.length() - 1) != "/" && outputPath.substr(outputPath.length() - 1) != "\\"){
			outputPath.append("/");
		}
		firstImgNum = atoi(argv[5]);
		lastImgNum = atoi(argv[6]);
		center_x = atof(argv[7]);
		center_y = atof(argv[8]);
		ring_width = atoi(argv[9]);
		threshold = atof(argv[10]);
		verbose = atoi(argv[11]);

		for(int img = firstImgNum; img < lastImgNum + 1; img++){

			time_t t_start, t_end;
			float **image=0, **polar_image=0, **ring_image=0;
			inputName = getName(inputBase, img);
			outputName = getName(outputBase, img);

			//read in image, perform thresholding

			image = read_float_image(inputPath + inputName, &width, &height);    //pixel data is stored in the form image[row][col]

			if(!image){
				return 1;
			}
			//Translate Image to Polar Coordinates

			polar_image = polar_transform_bilinear(image, center_x, center_y, width, height, &pol_width, &pol_height, thresh_max, thresh_min);

			//Call Ring Algorithm

			time(&t_start);
			ring_filter(&polar_image, pol_height, pol_width, threshold, M_rad, M_azi, ring_width);
			time(&t_end);
			double seconds = difftime(t_end, t_start);
			printf("Time to perform Ring Filtering: %.2f sec. \n", seconds);

			//Translate Ring-Image to Cartesian Coordinates

			ring_image = inverse_polar_transform_bilinear(polar_image, center_x, center_y, pol_width, pol_height, width, height);

			//Subtract Ring-Image from Image
	
			for(int row=0; row < height; row++){
				for(int col=0; col < width; col++){
					image[row][col] -= ring_image[row][col];
				}
			}
	
			//Write out Corrected-Image
			write_float(image, outputPath + outputName, width, height);

		}
		return 0;
	}
}
