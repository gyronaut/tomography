
#include tiff_io.h

using namespace std;

TiffIO::TiffIO()
{
}
TiffIO::~TiffIO()
{
}

float** TiffIO::ReadFloatImage(string input_name, int* w_ptr, int* h_ptr)
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

void TiffIO::WriteFloat(float** corrected_image, string output_name, int width, int height)
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