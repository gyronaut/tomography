
#include "tiff_io.h"

using namespace std;

TiffIO::TiffIO()
{
}
TiffIO::~TiffIO()
{
}

float** TiffIO::readFloatImage(string image_name, int* w_ptr, int* h_ptr)
{
	int width, height;

	TIFF* tif = TIFFOpen(image_name.c_str(), "r");

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
	tsize_t* bc;
	tsize_t stripsize;
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

void TiffIO::writeFloatImage(float** image_rows, string output_name, int width, int height)
{
	float* output = (float *) calloc(width*height, sizeof(float));
	float* buffer = (float *) calloc(width, sizeof(float));
	for(int i = 0; i < height; i++){
		memcpy(buffer, image_rows[i], width*sizeof(float));
		for(int j = 0; j < width; j++){
			output[j+(i*width)] = buffer[j]; 
		}
	}
	TIFF* tif = TIFFOpen(output_name.c_str(), "w");
	if(!tif){
		return;
	}
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);  
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);    
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, 3); //floating point = 3
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	float *tiff_buffer = NULL;
	tsize_t linebytes =  width*sizeof(float);
	tiff_buffer =(float *)_TIFFmalloc(linebytes);
	
	for(uint32 row = 0; row < height; row++){
		for(int col = 0; col < width; col++){
			tiff_buffer[col] = output[col+row*width];
		}
		if (TIFFWriteScanline(tif, tiff_buffer, row, 0) < 0){
			break;
		}
	}
	TIFFClose(tif);
	if(tiff_buffer){
		_TIFFfree(tiff_buffer);
	}
	free(buffer);
	free(output);
}
