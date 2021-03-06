#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <cstdio>
//#include "tiff_io.h"
#include "tiff_io-win.h"

using namespace std;

string getName(string name_base, int img_num){
	stringstream stream;
	stream << name_base << img_num << ".tif";
	return stream.str();
}

int main (int argc, char** argv){
	if(argc != 9){
		fprintf(stdout, "Incorrect Usage");
	}else{
		TIFFSetWarningHandler(NULL);
		TIFFSetErrorHandler(NULL);
		int verbose;
		int first_img_num, last_img_num;
		string input_path, input_base, input_name;
		string output_path, output_base, output_name;
		int width, height;
		int brightfield_interval
		int j_start, j_end, j_rad, j_0;
		int step, step_max;

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
		brightfield_interval = atoi(argv[7]);
		verbose = atoi(argv[8]);

		if(verbose > 1 || verbose < 0){
			verbose = 1;
		}

		j_rad = 6;
		step_max = 5;
		
		// Main Iteration over specified images //

		for(int img_num = first_img_num; img_num <= last_img_num; img_num++){
			float **image=0, **filtered_image=0;

			float ratio;
			
			float* derivative_sum = 0;

			input_name = getName(input_base, img_num);
			output_name = getName(output_base, img_num);
			image = tiff_io->read16bitImage(input_path + input_name, &width, &height);
			if(!image){
				if(verbose) fprintf(stderr, "Error: unable to open file %s. Trying next file.\n", (input_path+input_name).c_str());
			}else{
				float* ratio_array = (float*) calloc(width, sizeof(float));

				if(verbose) fprintf(stdout, "Allocating memory for Sub-Sinogram.\n");

				float* sub_sinogram_block = (float*) calloc(height*(2*j_rad + 1), sizeof(float));
				float** sub_sinogram = (float**) calloc(height, sizeof(float*));

				sub_sinogram[0] = sub_sinogram_block;
				for(int i = 1; i < height; i++){
					sub_sinogram[i] = sub_sinogram[i-1] + (2*j_rad + 1);
				}

				if(verbose) fprintf(stdout, "Opened file %s...\n", (input_path+input_name).c_str());
				for(step = 1; step <= step_max; step++){
					j_start = step*j_rad;
					j_end = width - step*j_rad - 1;
					for(j_0 = j_start; j_0 <= j_end; j_0++){

						derivative_sum = (float*) calloc((2*j_rad + 1), sizeof(float));
						
						// set up sub_sinogram(n, j)
						float min = 100000;
						float max = -100000;
						for(int row = 0; row < height; row++){
							for(int sino_col = 0, img_col = j_0 - step*j_rad; sino_col < 2*j_rad+1; sino_col++, img_col += step){
								sub_sinogram[row][sino_col] = image[row][img_col];
								if(sub_sinogram[row][sino_col] > max) max = sub_sinogram[row][sino_col];
								if(sub_sinogram[row][sino_col] < min) min = sub_sinogram[row][sino_col];
							}
						}
						if(j_0 == j_end){
							FILE* subsinogram = fopen("D:\\subsinogram.log", "w");
							fprintf(subsinogram, "J_0: %d\n", j_0);
							for(int row = 0; row < height; row++){
								for(int sino_col = 0; sino_col < 2*j_rad+1; sino_col++){
									fprintf(subsinogram, "%f ", sub_sinogram[row][sino_col]);
								}
								fprintf(subsinogram, "\n");
							}
							fclose(subsinogram);
						}
						

						int total_deriv_sum = 0;
						
						//if(verbose)fprintf(stdout, "Min: %f, Max: %f /n", min, max);
						// normalize sub_sinogram(n, j)
						for(int row = 0; row < height; row++){
							for(int col = 0; col < 2*j_rad+1; col++){
								sub_sinogram[row][col] = (sub_sinogram[row][col] - min)/(max - min);
							}
							// calculate D(n, j) and sum for y(j)
							for(int col = 1; col < 2*j_rad; col++){
								derivative_sum[col] += (2*sub_sinogram[row][col] - sub_sinogram[row][col-1] - sub_sinogram[row][col+1]);
							}
						}

						if(j_0 == j_end){
							FILE* norm_subsinogram = fopen("D:\\norm_subsinogram.log", "w");
							fprintf(norm_subsinogram, "J_0: %d\n", j_0);
							for(int row = 0; row < height; row++){
								for(int sino_col = 0; sino_col < 2*j_rad+1; sino_col++){
									fprintf(norm_subsinogram, "%f ", sub_sinogram[row][sino_col]);
								}
								fprintf(norm_subsinogram, "\n");
							}
							fclose(norm_subsinogram);
						}

						for(int col = 0; col < 2*j_rad+1; col++){
							total_deriv_sum += fabs(derivative_sum[col]);
						}
						if(((derivative_sum[j_rad] <= 0)==(derivative_sum[j_rad-1] >= 0)) && ((derivative_sum[j_rad] <= 0)==(derivative_sum[j_rad+1] >=0)) && (fabs(derivative_sum[j_rad]) > fabs(derivative_sum[j_rad - 1])) && (fabs(derivative_sum[j_rad]) > fabs(derivative_sum[j_rad + 1])) ){
							ratio = (2*j_rad - 4)*(fabs(derivative_sum[j_rad]) + fabs(derivative_sum[j_rad-1]) + fabs(derivative_sum[j_rad+1]))/(3*(total_deriv_sum - (fabs(derivative_sum[j_rad]) + fabs(derivative_sum[j_rad-1]) + fabs(derivative_sum[j_rad+1]))));
							if(ratio > ratio_array[j_0]) ratio_array[j_0] = ratio;
						}
					}
				}

				if(verbose) fprintf(stdout, "Writing out logs.\n");
				FILE* fp = fopen("D:\\ratio_log.txt", "w");
				if(!fp){
					if(verbose) fprintf(stderr, "Error opening the log file, aborting log.\n");
				}else{
					for(int i = 0; i < width; i++){
						fprintf(fp, "%f", ratio_array[i]);
						if(ratio_array[i] > 10){
							fprintf(fp, "    !!!");
						}
						fprintf(fp, "\n");
					}
				}
				FILE* derivlog = fopen("D:\\derivative_log.txt", "w");
				if(!derivlog){
					if(verbose) fprintf(stderr, "Error opening the log file, aborting log.\n");
				}else{
					for(int i = 0; i < 2*j_rad-1; i++){
						fprintf(derivlog, "%f \n", derivative_sum[i]);
					}
				}
				tiff_io->writeFloatImage(image, output_path+output_name, width, height);
			}
		}
	}
	return 0;
}
