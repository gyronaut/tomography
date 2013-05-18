/* Class for doin transformations on an image (polar & inverse polar, with both nearest neighbor
 * and bilinear interpolation)
 */

 #ifndef IMAGE_TRANSFORMS_H
 #define IMAGE_TRANSFORMS_H
 
 #pragma once
 
 class ImageTransformClass
 {
 private:
	int FindMinDistanceToEdge(float center_x, float center_y, int width, int height);
 public:
	ImageTransformClass();
	~ImageTransformClass();
	float** PolarTransform(float** image, float center_x, float center_y, int width, int height,
						   int* pol_width, int* pol_height, float thresh_max, float thresh_min);
	float** PolarTransformBilinear(float** image, float center_x, float center_y, int width,
								   int height, int* pol_width, int* pol_height, float thresh_max,
								   float thresh_min);
	float** InversePolarTransform(float** polar_image, float center_x, float center_y,
								  int pol_width, int  pol_height, int width, int height);
	float** InversePolarTransformBilinear(float** polar_image, float center_x, float center_y,
										  int pol_width, int  pol_height, int width, int height);
	
 };
 #endif