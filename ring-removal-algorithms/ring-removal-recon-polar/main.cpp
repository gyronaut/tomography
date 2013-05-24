
#include <string>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include "time.h"
#include "tiff.h"
#include "tiffio.h"

#include "tiff_io.h"
#include "image_transforms.h"
#include "image_filters.h"

using namespace std;

int round(float x){
	return (x > 0.0) ? floor(x+0.5) : ceil(x+0.5);
}

/* Stuff for trying to get a faster median sort working, didn't really pan out.
 *
void ElemSwapInt(int* arr, int index1, int index2){
	int store1 = arr[index1];
	arr[index1] = arr[index2];
	arr[index2] = store1;
}


int partition_2(float* median_array, int* position_array, int left, int right, int pivot_index){
	float pivot_value = median_array[pivot_index];
	lemSwap(median_array, pivot_index, right);
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
 *
 */

void doRingFilter(float*** polar_image, int pol_height, int pol_width, float threshold, int m_rad, int m_azi, int ring_width, ImageFilterClass* filter_machine){
	float* image_block = (float *) calloc(pol_height*pol_width, sizeof(float ));
	float** filtered_image = (float **) calloc(pol_height, sizeof(float*));
	filtered_image[0] = image_block;
	for(int i=1; i<pol_height; i++){
		filtered_image[i] = filtered_image[i-1]+pol_width;
	}
	printf("Pol_width: %d, pol_height: %d\n", pol_width, pol_height);	
	//Do radial median filter to get filtered_image
	printf("Performing Radial Filter on polar image... \n");
			
	filter_machine->doMedianFilter1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width/3 -1, 'x', m_rad, ring_width, pol_width, pol_height);

	filter_machine->doMedianFilter1D(&filtered_image, polar_image, 0, pol_width/3, pol_height-1, 2*pol_width/3 -1, 'x', 2*m_rad/3, ring_width, pol_width, pol_height);

	filter_machine->doMedianFilter1D(&filtered_image, polar_image, 0, 2*pol_width/3, pol_height-1, pol_width-1, 'x', m_rad/3, ring_width, pol_width, pol_height);

	//subtract filtered image from polar image to get difference image & do last thresholding

	printf("Calculating Difference Image... \n");
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col <  pol_width; col++){
			polar_image[0][row][col] -= filtered_image[row][col];
			if(polar_image[0][row][col] > threshold || polar_image[0][row][col] < -threshold){
				polar_image[0][row][col] = 0;
			}
		}
	}
	
	/* Do Azimuthal filter #2 (faster mean, does whole column in one call)
	 * using different kernel sizes for the different regions of the image (based on radius)
	 */

	printf("Performing Azimuthal mean filter... \n");
	fflush(stdout);
	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width/3-1, 'y', m_azi, pol_width, pol_height);
	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, pol_width/3, pol_height-1, 2*pol_width/3-1, 'y', 2*m_azi/3, pol_width, pol_height);
	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, 2*pol_width/3, pol_height-1, pol_width-1, 'y', m_azi/3, pol_width, pol_height);

	printf("Setting polar image equal to final ring image.. \n");
	fflush(stdout);
	//Set "polar_image" to the fully filtered data
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col < pol_width; col++){
			polar_image[0][row][col] = filtered_image[row][col];
		}
	}
	free(image_block);
	free(filtered_image);
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
		
		ImageFilterClass* filter_machine = new ImageFilterClass();
		ImageTransformClass* transform_machine = new ImageTransformClass();
		TiffIO* tiff_io = new TiffIO();
		
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

			image = tiff_io->readFloatImage(input_path + input_name, &width, &height);    //pixel data is stored in the form image[row][col]

			if(!image){
				return 1;
			}
			//Translate Image to Polar Coordinates

			polar_image = transform_machine->polarTransformBilinear(image, center_x, center_y, width, height, &pol_width, &pol_height, thresh_max, thresh_min);

			//Call Ring Algorithm

			time(&t_start);
			doRingFilter(&polar_image, pol_height, pol_width, threshold, m_rad, m_azi, ring_width, filter_machine);
			time(&t_end);
			double seconds = difftime(t_end, t_start);
			printf("Time to perform Ring Filtering: %.2f sec. \n", seconds);

			//Translate Ring-Image to Cartesian Coordinates
			printf("Doing inverse polar transform...\n");
			ring_image = transform_machine->inversePolarTransformBilinear(polar_image, center_x, center_y, pol_width, pol_height, width, height);

			//Subtract Ring-Image from Image
			printf("Subtracting rings from original image...\n");	
			for(int row=0; row < height; row++){
				for(int col=0; col < width; col++){
					image[row][col] -= ring_image[row][col];
				}
			}
	
			//Write out Corrected-Image
			printf("Writing out corrected image to %s.\n", (output_path+output_name).c_str());
			tiff_io->writeFloatImage(image, output_path + output_name, width, height);
			
			free(ring_image[0]);
			free(ring_image);
			free(polar_image[0]);
			free(polar_image);
			free(image[0]);
			free(image);
		}
		return 0;
	}
}
