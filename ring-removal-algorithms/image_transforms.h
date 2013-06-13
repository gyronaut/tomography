/* Class for doin transformations on an image (polar & inverse polar, with both nearest neighbor
 * and bilinear interpolation)
 */

 #ifndef IMAGE_TRANSFORMS_H
 #define IMAGE_TRANSFORMS_H

 #include <cmath>
 #include <cstdlib>
 #include <cstdio>
 
 #define PI 3.14159265359
 
 #pragma once
 
 class ImageTransformClass
 {
 private:
	int round(float x);
	int findMinDistanceToEdge(float center_x, float center_y, int width, int height);
 public:
	ImageTransformClass();
	~ImageTransformClass();
	float** polarTransform(float** image, float center_x, float center_y, int width, int height,
                               int* pol_width, int* pol_height, float thresh_max, float thresh_min,
                               int r_scale, int ang_scale, int overhang);
	float** polarTransformBilinear(float** image, float center_x, float center_y, int width,
                                       int height, int* pol_width, int* pol_height, float thresh_max,
                                       float thresh_min, int r_scale, int ang_scale, int overhang);
	float** inversePolarTransform(float** polar_image, float center_x, float center_y,
                                      int pol_width, int  pol_height, int width, int height,
                                      int r_scale, int over_hang);
	float** inversePolarTransformBilinear(float** polar_image, float center_x, float center_y,
                                              int pol_width, int  pol_height, int width, int height,
                                              int r_scale, int overhang);
	
 };
 
 #endif
