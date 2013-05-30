#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "image_filters.h"

int main()
{
	ImageFilterClass* filter_machine = new ImageFilterClass();
	float* number_array = (float*) calloc(100, sizeof(float));
	float** number_row = (float**)calloc(1, sizeof(float *));
	number_row[0] = number_array;
	float* output_array = (float*)calloc(100, sizeof(float));
	float** output_row = (float**)calloc(1, sizeof(float *));
	output_row[0] = output_array;
	float* fast_output_array = (float*)calloc(100, sizeof(float));
	float** fast_output_row = (float**)calloc(1, sizeof(float *));
	fast_output_row[0] = fast_output_array;
	for(int i = 0; i<100; i++){
		number_array[i] = float(i%10);
	}
	filter_machine->doMedianFilter1D(&output_row, &number_row, 0, 5, 0, 94, 'x', 2, 5, 100, 1);
	filter_machine->doMedianFilterFast1D(&fast_output_row, &number_row, 0, 5, 0, 94, 'x', 2, 5, 100, 1);

	FILE* out = fopen("filter_test_output", "w");
	fprintf(out, "Original, Median, Fast Median \n");
	for(int i = 0; i<100; i++){
		fprintf(out, "   %f     %f     %f \n", number_row[0][i], output_row[0][i], fast_output_row[0][i]);
	}
	fclose(out);
}
