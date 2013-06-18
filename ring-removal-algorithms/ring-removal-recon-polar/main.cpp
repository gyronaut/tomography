
#include <string>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include "time.h"

#include "tiff_io-win.h"
#include "image_transforms.h"
#include "image_filters.h"

using namespace std;

int round(float x){
	return (x != 0.0) ? floor(x+0.5) : 0;
}


void doRingFilter(float*** polar_image, int pol_height, int pol_width, float threshold, int m_rad, int m_azi, int ring_width, ImageFilterClass* filter_machine, int verbose){
	float* image_block = (float *) calloc(pol_height*pol_width, sizeof(float ));
	float** filtered_image = (float **) calloc(pol_height, sizeof(float*));
	filtered_image[0] = image_block;
	for(int i=1; i<pol_height; i++){
		filtered_image[i] = filtered_image[i-1]+pol_width;
	}
	if(verbose == 1) printf("Pol_width: %d, pol_height: %d\n", pol_width, pol_height);	
	//Do radial median filter to get filtered_image
	if(verbose == 1) printf("Performing Radial Filter on polar image... \n");
	clock_t start_median = clock();

//	filter_machine->doMedianFilterFast1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width-1, 'x', (ring_width -1)/2, ring_width, pol_width, pol_height);	
//	filter_machine->doMedianFilter1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width-1, 'x', (ring_width - 1)/2, ring_width, pol_width, pol_height);
		
	filter_machine->doMedianFilterFast1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width/3 -1, 'x', m_rad, ring_width, pol_width, pol_height);

	filter_machine->doMedianFilterFast1D(&filtered_image, polar_image, 0, pol_width/3, pol_height-1, 2*pol_width/3 -1, 'x', 2*m_rad/3, ring_width, pol_width, pol_height);

	filter_machine->doMedianFilterFast1D(&filtered_image, polar_image, 0, 2*pol_width/3, pol_height-1, pol_width-1, 'x', m_rad/3, ring_width, pol_width, pol_height);
	
	clock_t end_median = clock();
	if(verbose == 1) printf("Time for median filter: %f sec \n", (float(end_median - start_median)/CLOCKS_PER_SEC));

	//subtract filtered image from polar image to get difference image & do last thresholding

	if(verbose == 1) printf("Calculating Difference Image... \n");
	for(int row = 0; row < pol_height; row++){
		for(int col = 0; col <  pol_width; col++){
			polar_image[0][row][col] -= filtered_image[row][col];
			if(polar_image[0][row][col] > threshold || polar_image[0][row][col] < -threshold){
	//		if(polar_image[0][row][col] < threshold){
				polar_image[0][row][col] = 0;
			}
		}
	}
	
	/* Do Azimuthal filter #2 (faster mean, does whole column in one call)
	 * using different kernel sizes for the different regions of the image (based on radius)
	 */

	if(verbose == 1) printf("Performing Azimuthal mean filter... \n");
	clock_t start_mean = clock();

	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, 0, pol_height-1, pol_width/3-1, 'y', m_azi/3, pol_width, pol_height);
	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, pol_width/3, pol_height-1, 2*pol_width/3-1, 'y', 2*m_azi/3, pol_width, pol_height);
	filter_machine->doMeanFilterFast1D(&filtered_image, polar_image, 0, 2*pol_width/3, pol_height-1, pol_width-1, 'y', m_azi, pol_width, pol_height);

	clock_t end_mean = clock();
	if(verbose == 1) printf("Time for mean filtering: %f sec\n", (float(end_mean-start_mean)/CLOCKS_PER_SEC));

	if(verbose == 1) printf("Setting polar image equal to final ring image.. \n");
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
	if(!(argc == 15)){
		printf("\nUsage:\n\nring_remover_recon [input path] [output path] [input root] [output root] [first file num] [last file num] [center x y] [max ring width] [thresh min max] [ring threshold] [angular min] [verbose]\n\n");
		printf("      [input path]    Path to the folder containing the input images.\n");
		printf("     [output path]    Path to the folder to hold the filtered images.\n");
		printf(" [input file root]    Root name of the image to filter, with no file\n");
		printf("                      number or extension.\n");
		printf("[output file root]    Root of the output file name (number and file extension\n");
		printf("                      will be appended to this).\n");
		printf("  [first file num]    Number of first file to filter.\n");
		printf("   [last file num]    Number of last file to filter.\n");
		printf("    [center: x, y]    X and Y value of the center of rotation (floats)\n");
		printf("  [max ring width]    maximum width of the rings to be filtered in pixels (int).\n");
		printf("  [thresh min max]    min and max values for portion of image to filter (floats).\n");
		printf("  [ring threshold]    Rings are treated as offsets to the corrected image. This\n");
		printf("                      threshold is the maximum value of an offset due to a\n");
		printf("                      ring artifact (float).\n");
		printf("     [angular min]    minimum angle in degrees (int) to be considered ring artifact\n");
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
		int m_azi;
		int angular_min;
		int ring_width = 25;
		int r_scale = 2;
		int ang_scale = 2;
		string input_base, input_name, input_path;
		string output_base, output_name, output_path;
	//	default values for lego only
		float center_x=511, center_y=511, thresh_max=0.0018, thresh_min=0.0006, threshold = 0.00034;
	//	defualt values for large rings only
	//	float center_x=1240.5, center_y=1240.5, thresh_max = 20, thresh_min = -20, threshold = 6;
	//	default values for test case:
	//	float center_x, center_y, thresh_max = 20, thresh_min = -20, threshold = 1.0;
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
		ring_width = atoi(argv[9])*r_scale;
		thresh_min = atof(argv[10]);
		thresh_max = atof(argv[11]);
		threshold = atof(argv[12]);
		angular_min = atoi(argv[13]);
		verbose = atoi(argv[14]);

		for(int img = first_img_num; img < last_img_num + 1; img++){

			clock_t start = clock();
			float **image=0, **polar_image=0, **ring_image=0;
			input_name = getName(input_base, img);
			output_name = getName(output_base, img);

			//read in image, perform thresholding

			image = tiff_io->readFloatImage(input_path + input_name, &width, &height);    //pixel data is stored in the form image[row][col]

			if(!image){
				return 1;
			}
			//Translate Image to Polar Coordinates
			if(verbose == 1) printf("Performing Polar Transformation...\n");
			clock_t start_polar = clock();
//			polar_image = transform_machine->polarTransformBilinear(image, center_x, center_y, width, height, &pol_width, &pol_height, thresh_max, thresh_min, r_scale, ang_scale, ring_width);
			polar_image = transform_machine->polarTransform(image, center_x, center_y, width, height, &pol_width, &pol_height, thresh_max, thresh_min, r_scale, ang_scale, ring_width);
			clock_t end_polar = clock();
			if(verbose == 1) printf("Time for polar Transformation: %f sec\n", (float(end_polar - start_polar)/CLOCKS_PER_SEC));
			m_azi = ceil(float(pol_height)/(360.0))*angular_min;	
			//Call Ring Algorithm

			doRingFilter(&polar_image, pol_height, pol_width, threshold, m_rad, m_azi, ring_width, filter_machine, verbose);
						
			//Translate Ring-Image to Cartesian Coordinates
			if(verbose == 1) printf("Doing inverse polar transform...\n");
			clock_t start_invpol = clock();
//			ring_image = transform_machine->inversePolarTransformBilinear(polar_image, center_x, center_y, pol_width, pol_height, width, height, r_scale, ring_width);
			ring_image = transform_machine->inversePolarTransform(polar_image, center_x, center_y, pol_width, pol_height, width, height, r_scale, ring_width);
			clock_t end_invpol = clock();
			if(verbose == 1) printf("Time for Inverse Polar Transform: %f sec\n", (float(end_invpol - start_invpol)/CLOCKS_PER_SEC));

			//Subtract Ring-Image from Image
			if(verbose == 1) printf("Subtracting rings from original image...\n");	
			for(int row=0; row < height; row++){
				for(int col=0; col < width; col++){
					image[row][col] -= ring_image[row][col];
				}
			}
	
			//Write out Corrected-Image
			if(verbose == 1) printf("Writing out corrected image to %s.\n", (output_path+output_name).c_str());
			tiff_io->writeFloatImage(image, output_path + output_name, width, height);
			clock_t end = clock();
			if(verbose == 1) printf("Total time to perform ring filtering: %f sec\n", (float(end-start))/CLOCKS_PER_SEC);
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
