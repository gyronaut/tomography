
#include "image_transforms.h"


ImageTransformClass::ImageTransformClass(){};
ImageTransformClass::~ImageTransformClass(){};

int ImageTransformClass::findMinDistanceToEdge(float center_x, float center_y, int width,
											   int height)
{
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

int ImageTransformClass::round(float x){
	return (x != 0.0) ? floor(x+0.5) : 0;
}

float** ImageTransformClass::polarTransform(float** image, float center_x, float center_y, int width,
											int height, int* pol_width, int* pol_height,
											float thresh_max, float thresh_min)
{
	int max_r = findMinDistanceToEdge(center_x, center_y, width, height);
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
float** ImageTransformClass::polarTransformBilinear(float** image, float center_x, float center_y,
													int width, int height, int* p_pol_width,
													int* p_pol_height, float thresh_max,
													float thresh_min)
{
	int max_r = findMinDistanceToEdge(center_x, center_y, width, height);
	int pol_width = max_r;
	int pol_height = round(4.0*PI*float(max_r));
	
	*p_pol_width = pol_width;
	*p_pol_height = pol_height;

	float* image_block = (float *) calloc(pol_height*pol_width, sizeof(float));	
	float** polar_image = (float **) calloc(pol_height, sizeof(float *));
	polar_image[0] = image_block;
	for(int i=1; i<pol_height; i++){
		polar_image[i] = polar_image[i-1] + pol_width;
	}

	float theta, fl_x, fl_y, value;
	int x_1, x_2, y_1, y_2;
	for(int row = 0; row<pol_height; row++){
		for(int r = 0; r<pol_width; r++){
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

float** ImageTransformClass::inversePolarTransform(float** polar_image, float center_x,
												  float center_y, int pol_width, int  pol_height,
												  int width, int height)
{
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
float** ImageTransformClass::inversePolarTransformBilinear(float** polar_image, float center_x,
                                                           float center_y, int pol_width,
                                                           int pol_height, int width, int height)
{
	float* image_block = (float *) calloc(height*width, sizeof(float));
	float** cart_image = (float **) calloc(height, sizeof(float *));
	cart_image[0] = image_block;
	for(int i = 1; i < height; i++){
		cart_image[i] = cart_image[i-1]+width;
	}
	float r, theta, value, row_float;
	int col_1, col_2, row_1, row_2, row_floor, row_ceil;
	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++){
			value = 0;
			r = sqrt(float((x - center_x)*(x - center_x)) + float((y - center_y)*(y - center_y)));
			if(r < pol_width - 1){
				col_1 = floor(r);
				col_2 = ceil(r);
				theta = atan2(float(y - center_y), float(x - center_x));
				if(theta < 0){
					theta += 2*PI;
				}
				row_float = theta*float(pol_height -1)/(2.0*PI);
				if(row_float > pol_height){
					row_float -= pol_height;
				}
				row_floor = floor(row_float);
				row_ceil = ceil(row_float);
				row_1 = row_floor%(pol_height);
				row_2 = row_ceil%(pol_height);
				if(col_1 == col_2){
					if(row_1 == row_2){
						value = polar_image[row_1][col_1];
					}else{
						value = (row_ceil - row_float)*polar_image[row_1][col_1] + (row_float - row_floor)*polar_image[row_2][col_1];
					}
				}else if(row_1 == row_2){
						value = (col_2 - r)*polar_image[row_1][col_1] + (r - col_1)*polar_image[row_1][col_2];
				}else{
					value = (row_ceil - row_float)*((col_2 - r)*polar_image[row_1][col_1] + (r - col_1)*polar_image[row_1][col_2]) + (row_float - row_floor)*((col_2 - r)*polar_image[row_2][col_1] + (r - col_1)*polar_image[row_2][col_2]);
				}
			}
			cart_image[y][x] = value;	
		}
	}
	return cart_image;
}
